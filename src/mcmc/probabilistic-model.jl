##############################
# Probabilistic model



## Model parameters to fit
####################




"""
$(TYPEDEF)
$(FIELDS)

Type to hold information about the fitting parameters. Should be
flexible to easily change what parameters to fit.

A sample of the parameters is just a vector but can be interpreted
together with this type.  This is what is done in `make_fields1d`,
`make_loglikelihood` and `make_fwdm_fn`.

Contains:
- th0: initial parameter values
  - for the 1D fields btilde, fsl and temp, th0 corresponds to the
    deviation from the prior expectation, which is in gb.fsl, etc.
  - for other parameters it is just the straight value.
- names: name of parameters with their indices into th0 (to allow
  groups of parameters with the same name)
- th_ele: elevation of parameters. If just one parameter, then it is
  constant, otherwise linearly interpolated to get values at all
  elevations.
- gb::Bands: reference to the glacier Bands which this corresponds to

Convenience constructors (all slow)
- `Theta0(names, counts, gb, pp, pm, pn[, th])` constructs a theta0 with given names.
   Defaults th to zero-vectors.
- `Theta0(dict, gb, pp, pm, pn)` dict maps names to para-values
"""
struct Theta0{T}
    th0::Vector{T}  # the initial parameters
    names::OrderedDict{Symbol,UnitRange{Int64}} # name and index into th
    th_ele::OrderedDict{Symbol,LinSpace{Float64}} # elevation of parameters
                                           # if they are spatially
                                           # distributed (in 1D).

    ## Work arrays: these are updated
    # 1d fields defined on elevation bands gb.bands.  These are fed to
    # the forward model.
    gb::Bands
    pp::Phys
    pm::MPara
    pn::Num
    function Theta0{T}(th0, names, th_ele, gb, pp, pm, pn) where T
        @assert length(unique(keys(names)))==length(names)
        @assert length(th0)==mapreduce(length, +, 0, values(names))
        new{T}(th0, names, th_ele, gb, pp, pm, pn)
    end
end
# this constructor is slow, but the most convenient constructor.
function Theta0(names, counts, gb::Bands, pp::Phys, pm::MPara, pn::Num, th0::AbstractVector{T}=zeros(F,sum(count))) where T
    d = OrderedDict{Symbol,UnitRange{Int64}}()
    ii = 1
    for (name,count) in zip(names, counts)
        d[Symbol(name)] = ii:ii+count-1
        ii = ii+count
    end
    # make ranges of 1D fields thetas:
    th_ele = OrderedDict{Symbol,LinSpace{Float64}}()
    for name in names
        inds = d[name]
        li = length(inds)
        if li>1
            th_ele[name] = linspace(first(gb.ele), last(gb.ele), li)
        elseif li==1
            mid = (first(gb.ele) + last(gb.ele))/2
            th_ele[name] = linspace(mid,mid,1)
        else
            error("This should not happen...")
        end
    end

    return Theta0{T}(th0, d, th_ele, gb, pp, pm, pn)
end
function Theta0(d::Associative, gb::Bands, pp, pm, pn)
    names = collect(keys(d))
    counts = [length(v) for v in values(d)]
    th0 = F[]
    for vv in values(d)
        append!(th0, vv)
    end
    return Theta0(names, counts, gb, pp, pm, pn, th0)
end
Base.show(io::IO, theta0::Theta0) = println(io, "thetas with n=$(length(theta0.th0)) and keys:\n $((keys(theta0.names)...,))")
function get_varnames(theta0::Theta0; pad=true, addnumbers=true)
    varnames = String[]
    for name in keys(theta0.names)
        if length(theta0[name])==1
            push!(varnames, "$(string(name))")
        else
            for i=1:length(theta0[name])
                push!(varnames, addnumbers ? "$(string(name))$i" : "$(string(name))")
            end
        end
    end
    # pad with spaces:
    if pad
        maxlen = mapreduce(length, max, 0, varnames)
        return map(s->rpad(s, maxlen, " "), varnames)
    else
        return varnames
    end
end
"""
Get elevation for index i
"""
function get_ele(t0::Theta0, i::Int)
    for (n,inds) in t0.names
        if i in inds
            ele = t0.th_ele[n][findfirst(i.==inds)]
            return ele
        end
    end
    error("Index $i not in $t0")
end
# Array interface to directly index into thetas:
Base.IndexStyle(::Type{Theta0}) = IndexLinear()
Base.size(a::Theta0) = size(a.th0)
# Base.getindex(a::Theta0, i::Int) = a.th0[i]
# Base.setindex!(a::Theta0, v, i::Int) = a.th0[i]=v
Base.getindex(a::Theta0, s::Symbol) = @view a.th0[a.names[s]]

"""
Enumeration for the different types of model fitting

- h : fit to thickness only
- iv : fit to IV only
- h_iv : fit to h & IV
- length : fit to glacier length only (also included in h, iv and h_iv)
           Note: using this only necessitates the 1D model evaluation which is about
           100x faster than the 2d model.
- prior : fit to prior only
"""
baremodule FitTarget
using Base
@enum _FitTarget h iv h_iv length prior
end


#####
# Forward model wrapping
#####
"""
    make_fwdm_fn(theta0::Theta0, pm, pp, pn)

Create a function which runs the forward model for a theta as
argument.  One for the 2D model and one for the 1D, the latter only
runs the 1D part and does not produce a blob (i.e. ice thickness and
velocity maps).
"""
function make_fwdm_fn(theta0::Theta0, pp, pm, pn)
    gb = theta0.gb
    gl = gb.gl
    pn = Num(pn, verbose=false)
    pm = deepcopy(pm) # this gets updated in-place

    # check that we know what to set
    for (name,_) in theta0.names
        if !(name in [:btilde, :fsl, :temp]) && !isdefined(pm,name)
            error("Not implemented setting non-MPara parameters. Offending name: $name")
        end
    end

    # pre-calculate an update structure for pm
    args_to_MPara = []
    inds = Tuple{Int,Int}[]
    for (ii,name) in enumerate(fieldnames(pm))
        if haskey(theta0.names,name)
            its = theta0.names[name]
            push!(inds, (ii, its[1]))
            push!(args_to_MPara, nothing)
        else
            push!(args_to_MPara, getfield(pm,name))
        end
    end

    # https://github.com/JuliaLang/julia/issues/15276
    fwd = let inds = inds, args_to_MPara = args_to_MPara, pm = pm,
        theta0 = theta0, gb = gb, pp = pp, pn = pn
        function (theta::Vector)
            # modify pm to allow fitting of its parameters
            for (ias, its) in inds
                args_to_MPara[ias] = theta[its]
            end
            pm_ = MPara(args_to_MPara...)::MPara

            # calculate 1d fields (TODO: maybe use some work-arrays here)
            btilde, temp, fsl = make_fields1d(theta, theta0)

            # TODO: update other parameters?

            # run forward model:
            fwdsol = fwdm(gb, pp, pm_, pn, btilde, fsl, temp)
            return fwdsol
        end
    end

    fwd1d = let inds = inds, args_to_MPara = args_to_MPara, pm = pm,
        theta0 = theta0, gb = gb, pp = pp, pn = pn
        function (theta::Vector)
            # modify pm to allow fitting of its parameters
            for (ias, its) in inds
                args_to_MPara[ias] = theta[its]
            end
            pm_ = MPara(args_to_MPara...)::MPara

            # calculate 1d fields (TODO: maybe use some work-arrays here)
            btilde, temp, fsl = make_fields1d(theta, theta0)

            # TODO: update other parameters?

            # run 1D forward model:
            hs_n, tau_local, tau_mean_n, iv_n, qtot_n = fwdm1d(gb, btilde, fsl, temp, pp, pm_, pn)
            return hs_n, tau_local, tau_mean_n, iv_n, qtot_n, pm_
        end
    end
    return fwd, fwd1d
end

"""
    btilde, fsl, temp = make_fields1d(theta, theta0, names=[:btilde, :temp, :fsl])

Make the 1D fields according to the values in theta0.gb (their
expectation) and theta (their error/correction).

Note that the elevation decreases with parameter index.

TODO: make this more general to allow others to be 1D, say `dist_exp`.
"""
function make_fields1d(theta, theta0::Theta0, names=[:btilde, :temp, :fsl])
    gb = theta0.gb
    out = typeof(gb.bdot)[]
    for name in names
        if haskey(theta0.names, name)
            inds = theta0.names[name]
            li = length(inds)
            ele_theta0 = theta0.th_ele[name]
            f = piecewiselinear(ele_theta0, theta[inds]) # TODO use something other than linear interpolation?
            mean_val = getfield1d(gb, name) # this is the expectation value of the field
            push!(out, f.(gb.ele) .+ mean_val) # add error to expectation
        else # not fitting this parameter
            push!(out, getfield1d(gb, name)) # just return (expectation) value of the field
        end
    end
    return out
end
make_fields1d(theta, theta0::Theta0, name::Symbol) = make_fields1d(theta, theta0, [name])[1]

##################
## Probabilistic model
##################

"""
    get_sigmas(gl::Glacier, pm::MPara) -> sigma_h, sigma_iv

Return total sigmas, i.e. both observational (in gl.iv of gl.h) and
model-sigma (in pm.sigma_*).
"""
function get_sigmas(gl::Glacier, pm::MPara)
    @unpack sigma_h_model, sigma_iv_model = pm
    sigma_h_obs, sigma_iv_obs = gl.h.sigma, gl.iv.sigma
    sigma_h = sqrt(sigma_h_model^2 + sigma_h_obs^2)
    sigma_iv = sqrt(sigma_iv_model^2 + sigma_iv_obs^2)

    return sigma_h, sigma_iv
    # inf = convert(typeof(sigma_iv), Inf)

    # return if fit_target==FitTarget.h_iv
    #     sigma_h, sigma_iv
    # elseif fit_target==FitTarget.h
    #     sigma_h, inf
    # elseif fit_target==FitTarget.iv
    #     inf, sigma_iv
    # elseif fit_target==FitTarget.prior || fit_target==FitTarget.length
    #     inf, inf
    # else
    #     error("Unknown `fit_target`==$fit_target")
    # end
end

############
# Posterior
############
"""
$(SIGNATURES)

Take loglikelihood and logprior and return a function evaluating the
log-posterior and a blob (additional information produced during a
forward model run): `logposterior(theta)`

Note: if the log-prior is -Inf then then log-likelihood is not
evaluated and a empty blob is returned.
"""
function make_logposterior(loglikelihood, loglikelihood1d, logprior)
    f = function (theta)
        lp = logprior(theta)
        # return now if prior==-Inf, this way the forward model is not run
        # with invalid parameters.
        if lp==-Inf
            # return empty blob but needs to be type-stable:
            return -Inf, (NaN, NaN, Array{F}(0), Array{F}(0), Array{F}(0,0), Array{F}(0,0))
        end
        ll, blob = loglikelihood(theta)
        return  ll + lp, blob
    end
    f1d = function (theta)
        lp = logprior(theta)
        # return now if prior==-Inf, this way the forward model is not run
        # with invalid parameters.
        if lp==-Inf
            return -Inf
        end
        ll = loglikelihood1d(theta)
        return  ll + lp
    end
    return f, f1d
end

###
# Likelihood
###

"""
$(SIGNATURES)

Constructs a log-likelihood function: `loglikelihood(theta)` and also
returns the effective sigma_h and sigma_iv

Kwargs:
- iv_comp: how to do the comparison to IV.  Defaults to `:at_bands` for
           2D IV data, and `:at_measurement_points` otherwise.

Does this by making a closure over theta0 (containing gb), pm, pp,
pn.  Do this in a function to avoid global variables.

NOTE: this runs into https://github.com/JuliaLang/julia/issues/15276.
However on Julia 0.6 this particular incarnation of 15276 seems to be
fixed.
"""
function make_loglikelihood(theta0::Theta0, pp::Phys, pm::MPara, pn::Num,
                            fit_target::FitTarget._FitTarget;
                            iv_comp = if theta0.gb.gl.iv.vals isa VAWTools.Gridded
                            :at_bands # this is what is usually used
                            elseif  theta0.gb.gl.iv.vals isa VAWTools.Traj
                            :at_measurement_points
                            elseif theta0.gb.gl.iv.vals isa Void
                            :no_iv
                            else
                            error("Not anticipated")
                            end
                            )
    gb = theta0.gb
    gl = gb.gl
    pn = Num(pn, verbose=false)
    pm = deepcopy(pm) # use a copy as this gets updated in-place

    # # some tests
    # !hasdata(gl.iv) && sigma_iv!=0 &&
    #     warn("$gl has no IV data but sigma_iv is non-zero.")
    # !hasdata(gl.h) && sigma_h!=0 &&
    #     warn("$gl has no thickness data but sigma_h is non-zero.")

    #  forward model function
    fwd_fn, fwd1d_fn = make_fwdm_fn(theta0, pp, pm, pn)

    # https://github.com/JuliaLang/julia/issues/15276
    # "Performance of captured variables in closures"
    f = let gl=gl, gb=gb, fwd_fn=fwd_fn, in_glacier_=in_glacier(gl), fit_target=fit_target, iv_comp=iv_comp,
            gl_h=gl.h, gl_iv=gl.iv  # these two are needed because those fields are not strongly-typed in Glacier

        function (theta::Vector)
            # run forward model:
            fwdsol = fwd_fn(theta)
            pm_ = fwdsol.pm # something might change in pm due to theta
            @unpack vol, vol_above, hs1d, ivs1d, hs2d, ivs2d, ivs2d_band, ivs2d_band_mask = fwdsol
            # make std-dev of h and iv errors:
            sigma_h, sigma_iv = get_sigmas(gl, pm_)

            blob = vol, vol_above, hs1d, ivs1d, hs2d, ivs2d
            # return log-likelihood and blob
            if iv_comp==:at_bands
                return _loglikelihood(gl.dem, gl.iscalving, in_glacier_, gl_h, gl_iv, gl.ivmask, hs2d,
                                      ivs2d_band, ivs2d_band_mask, sigma_h, sigma_iv, gb.x, hs1d, fit_target,
                                      pm_.error_model, pm_.nthin_h, pm_.nthin_iv), blob
            else
                return _loglikelihood(gl.dem, gl.iscalving, in_glacier_, gl_h, gl_iv, gl.ivmask, hs2d,
                                      ivs2d, gl.glaciermask, sigma_h, sigma_iv, gb.x, hs1d, fit_target,
                                      pm_.error_model, pm_.nthin_h, pm_.nthin_iv), blob
            end
        end
    end

    # 1D loglikelihood
    f1d = let gb=gb, fwd1d_fn=fwd1d_fn, fit_target=fit_target, iscalving=gl.iscalving,
              verbose=pn.verbose

        function (theta::Vector)
            # run 1D forward model:
            hs1d = fwd1d_fn(theta)[1]

            # no blob

            # likelihood
            x = gb.x
            if verbose && fit_target!=FitTarget.length
                warn("1D log-likelihood only works with FitTarget.length")
            end
            # TODO add a comparison to Vol-Area scaling?

            # @show length1d(x, hs1d), x[end]
            # @show length_logpdf(length1d(x, hs1d), x, iscalving)
            return length_logpdf(length1d(x, hs1d), x, iscalving)
        end
    end
    return f, f1d, fwd_fn, fwd1d_fn
end

"Return interpolation to test whether an arbitrary point is inside glacier"
in_glacier(gl) = Interp.interpolate((gl.dem.x, gl.dem.y), gl.glaciermask, Interp.Gridded(Interp.Constant()) )

"""
$(SIGNATURES)

Log-Likelihood for given measured-data and given forward model result.
Can be overloaded for specific glaciers.

This assumes that the errors have a normal distribution with std-dev
sigma.

Input
 - gl, ivs1d_measured, ivs1d, hs2d, ivs2d, sigma_h, sigma_iv

If sigmas==Inf then that error is not calculated.  If both sigma_iv and sigma_h
are Inf, then the log-likelihood==0.  This is equivalent of sampling
just the prior (but the blob will still be evaluated and included, i.e. this still runs
the forward model and thus is way slower than just sampling the prior).

Output:
- loglikelihood

Note the log of the product of all the independent Gaussians is:

    -n ( ln(σ) + 1/2*ln(2π) ) - 1/(2σ²) * Σ (x_i - μ_i)²
"""
function _loglikelihood(dem, iscalving, in_glacier, gl_h, gl_iv, gl_ivmask,
                        hs2d, ivs2d, ivs2d_mask, sigma_h, sigma_iv, x, hs1d, fit_target,
                        error_model,
                        nthin_h::Int, nthin_iv::Int) # bug in Julia: the two ::Int are needed for type-inference
    out = 0.0
    if fit_target==FitTarget.h || fit_target==FitTarget.h_iv
        if gl_h isa SomeData{F}
            # For some GlaThiDa glacier there is just a mean thickness value.  This uses this.
            Vmod = sum(hs2d)*step(dem.x)^2
            Vmeas = gl_h.vals
            out += -(Vmod-Vmeas)^2/(2*sigma_h^2) - log(sigma_h*sqrt(2pi))
        elseif gl_h isa SomeData{Traj{F}}
            if error_model==ErrorModel.sqrerr
                ## Normal dist
                te, ne = thick_err(dem, in_glacier, gl_h, hs2d, 1:length(gl_h.vals.splits), nthin_h, sqrerr, +)
                if ne>0
                    out += -te/(2*sigma_h^2) - ne*log(sigma_h*sqrt(2pi))
                end
            elseif error_model==ErrorModel.abserr
                ## Laplace dist (abs-errors), such that sigma is still the stdev
                te, ne = thick_err(dem, in_glacier, gl_h, hs2d, 1:length(gl_h.vals.splits), nthin_h, abserr, +)
                if ne>0
                    out += -te/(sigma_h/sqrt(2)) - ne*log(sigma_h*sqrt(2))
                end
            else
                error("ErrorModel $error_model not implemented")
            end
        else
            error("Not supported gl_h type: $(typeof(gl_h))")
        end
    end
    if fit_target==FitTarget.iv || fit_target==FitTarget.h_iv
        if error_model==ErrorModel.sqrerr
            ## Normal dist
            tiv2d, niv2d = iv2d_err(dem, in_glacier, gl_iv, gl_ivmask, ivs2d, ivs2d_mask, nthin_iv, sqrerr, +)
            if niv2d>0
                out += -tiv2d/(2*sigma_iv^2) - niv2d*log(sigma_iv*sqrt(2pi))
            end
        elseif error_model==ErrorModel.abserr
            ## Laplace dist (abs-errors)
            tiv2d, niv2d = iv2d_err(dem, in_glacier, gl_iv, gl_ivmask, ivs2d, ivs2d_mask, nthin_iv, abserr, +)
            if niv2d>0
                out += -tiv2d/(sigma_iv/sqrt(2)) - niv2d*log(sigma_iv*sqrt(2))
            end
        else
            error("ErrorModel $error_model not implemented")
        end
    end
    ## Make sure the length is correct
    # @show length1d(x, hs1d), x[end]
    # @show length_logpdf(length1d(x, hs1d), x, iscalving)
    out += length_logpdf(length1d(x, hs1d), x, iscalving)

    # TODO add a comparison to V-A scaling?

    return out
end

"""
    length_logpdf(theta::Vector, theta0::Theta0)

This encodes that a glacier should have the length of the domain.

TODO: make parameters adjustable.
"""
function length_logpdf(length, x, iscalving)
    if iscalving
        # then just make sure the glacier reaches the end
        return length<x[end] ? -Inf : 0.0
    else # land terminating, make sure glacier length is close to what it should be
        expectation = x[end]
        sigma = expectation*0.005
        percent = 0.03
        range_abs = [(1-percent) * expectation, (1+percent) * expectation]
        # using the subnormal make the probability drop much faster than Normal beyond range_abs:
        return logprior_scalar_subnormal(length, expectation, sigma, range_abs)
    end
end



# Define a function to quantify a model vs measurements error. The
# standard is to assume a normal distributed mean square error.
"Square error of one observation."
sqrerr(y_measured, y_modelled) = (y_measured-y_modelled)^2

## Note:
# other error functions could be used, as e.g. below.
# Also note the term "heteroscedasticity", which means that the error-sigma
# is dependent on the variable, e.g. relative error.

"Absolute error of one observation."
abserr(y_measured, y_modelled) = abs(y_measured-y_modelled)

"""
Relative square errors [1] using the trick of [2].

[1] https://ieeexplore.ieee.org/abstract/document/159804
[2] https://www.jstor.org/stable/27920135
"""
rel_sqrerr(y_measured, y_modelled) =
    ((y_measured-y_modelled)/y_measured)^2 +
    ((y_measured-y_modelled)/y_modelled)^2

"""
Relative square errors. Not sure this makes sense.
"""
rel_sqrerr_min(y_measured, y_modelled) = min(
    ((y_measured-y_modelled)/y_measured)^2,
    ((y_measured-y_modelled)/y_modelled)^2)


"""
Relative absolute error (LARE)
Chen et al 2010, Eq3
https://www.jstor.org/stable/27920135
"""
rel_abserr(y_measured, y_modelled) =
    abs((y_measured-y_modelled)/y_measured) +
    abs((y_measured-y_modelled)/y_modelled)

rel_abserr_min(y_measured, y_modelled) = min(
    abs((y_measured-y_modelled)/y_measured),
    abs((y_measured-y_modelled)/y_modelled))


"""
$(SIGNATURES)

Reduction (sum) of ice thickness (square) errors. Can be overloaded
for a particular glacier.

Input either:
- dem: gl.dem
- in_glacier : interpolation, true if in glacier
- gl_h: gl.h
- hs2d: 2D thickness model output
or:
- gl
- hs2d

Optional positional args:
- lines=1:length(gl_h.vals.splits) : which radar lines to process
- thin=10 : thin the radar lines.  Needed because of spatial correlation.
- errfn=sqrerr : the error function to use `(value, measurement) -> comparison`
- reduction=+ : how to reduce the errors

Return:
- reduction of errors (sum of square errors)
- number of points used

TODO
- think about weights: maybe two very closely spaced points should
  count less than two far ones?  -> there is some correlation
  between them, what's the length?
"""
function thick_err(dem, in_glacier, gl_h::SomeData{Traj{F}}, hs2d, lines=1:length(gl_h.vals.splits),
                   thin=10, errfn=sqrerr, reduction=+, neutral_element=0.0)
    ihs  = Interp.interpolate((dem.x, dem.y), hs2d, Interp.Gridded(Interp.Linear()) )
    out = neutral_element
    n = 0 # counter for number of points at which an error is evaluated
    start = 3
    step = thin # TODO hack...
    tr = gl_h.vals
    for line in lines
        for i = tr.splits[line][1]:step:tr.splits[line][end]
            x,y,h = tr.x[i],tr.y[i],tr.v[i]
            if in_glacier[x,y]==1 # only inside glacier
                err = errfn(ihs[x,y], h)
                isnan(err) && continue
                out = reduction(out, errfn(ihs[x,y], h))
                n += 1
            end
        end
    end
    return out, n
end
thick_err(gl::Glacier, hs2d::Matrix, lines=1:length(gl.h.vals.splits), thin=10, errfn=sqrerr, reduction=+, neutral_element=0.0) =
    thick_err(gl.dem, in_glacier(gl), gl.h, hs2d, lines, thin, errfn, reduction, neutral_element)
thick_err(gl::Glacier, hs2d::Gridded, lines=1:length(gl.h.vals.splits), thin=10, errfn=sqrerr, reduction=+, neutral_element=0.0) =
    thick_err(gl.dem, in_glacier(gl), gl.h, hs2d.v, lines, thin, errfn, reduction, neutral_element)

"""
$(SIGNATURES)

Sum of 2D ice velocity (square) errors. Can be overloaded for a particular
glacier.

Input:
- gl
- ivs2d: 2D velocity model output
- ivs2d_mask: where 2D velocity model output should be compared to observations
- gl_iv: IV observations
- gl_ivmask: and mask

Optional:
- thin=10: thin the velocity measurements.  Needed because of correlation.

"""
function iv2d_err(dem, in_glacier, gl_iv::SomeData{Gridded{F}}, gl_ivmask, ivs2d, ivs2d_mask,
                  thin=10, errfn=sqrerr, reduction=+)
    iv_ = gl_iv.vals
    iv_obs, flag  = size(iv_.v)==size(ivs2d) ? (iv_.v,true) :
        (Interp.interpolate((iv_.x, iv_.y), iv_.v, Interp.Gridded(Interp.Linear()) ), false)
    iv_obs_mask = flag ? gl_ivmask :
        Interp.interpolate((iv_.x, iv_.y), gl_ivmask, Interp.Gridded(Interp.Constant()) )

    out = 0.0
    n = 0
    step = thin # TODO hack...
    for iy=1:step:length(dem.y)
        y = dem.y[iy]
        for ix=1:step:length(dem.x)
            # skip points not in this mask
            !(ivs2d_mask[ix,iy]) && continue
            x = dem.x[ix]

            # check that observation-mask is true:
            if (flag ? iv_obs_mask[ix,iy] : iv_obs_mask[x,y]==1)
                obs = (flag ? iv_obs[ix,iy] : iv_obs[x,y])::Float64
                out = reduction(out, errfn(obs, ivs2d[ix,iy]))
                n += 1
            end
        end
    end
    return out, n
end
# for trajectory data
function iv2d_err(dem, in_glacier, gl_iv::SomeData{Traj{F}}, gl_ivmask, ivs2d, ivs2d_mask,
                  thin=1, errfn=sqrerr, reduction=+)
    iiv  = Interp.interpolate((dem.x, dem.y), ivs2d, Interp.Gridded(Interp.Linear()) )
    out = 0.0
    n = 0
    step = thin # TODO hack...
    for i=1:step:length(gl_iv.vals.v)
        x = gl_iv.vals.x[i]
        y = gl_iv.vals.y[i]
        iv = gl_iv.vals.v[i]
        if in_glacier[x,y]==1 # only inside glacier
            out = reduction(out, errfn(iiv[x,y], iv))
            n += 1
        end
    end
    return out, n
end
iv2d_err(dem, in_glacier, gl_iv::SomeData{Void}, gl_ivmask, ivs2d, ivs2d_mask,
         thin=10, errfn=sqrerr, reduction=+) = 0.0, 0

iv2d_err(gl, ivs2d::Matrix, ivmask=gl.ivmask, thin=10, errfn=sqrerr, reduction=+) =
    iv2d_err(gl.dem, in_glacier(gl), gl.iv, ivmask, ivs2d, thin, errfn, reduction)
iv2d_err(gl, ivs2d::Gridded, ivmask=gl.ivmask, thin=10, errfn=sqrerr, reduction=+) =
    iv2d_err(gl.dem, in_glacier(gl), gl.iv, ivmask, ivs2d.v, thin, errfn, reduction)


## Not used anymore:
# "Sum of 1D ice velocity errors.  This does not work well."
# function iv1d_err(ivs1d_measured, ivs1d, thin=3)
#     out = 0.0
#     n = 0
#     step = thin
#     for i=1:step:length(ivs1d_measured)
#         iv1 = ivs1d_measured[i]
#         iv2 = ivs1d[i]
#         if iv1>0 # mask
#             out += sqrerr(iv1,iv2)
#             n += 1
#         end
#     end
#     return out, n
# end
