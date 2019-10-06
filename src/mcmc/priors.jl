############
# Priors
############
#
# need logprior(theta) = p(.)*p(.)
#
# assume indpendence for now

# priors = Dict{Symbol, Union{Function, Vector{Function}}()

"""
$(SIGNATURES)

Combines the independent component-wise specified priors in
`component_wise_priors` with any free-form priors in
`full_theta_priors`.

`component_wise_priors` is a `Dict{Symbol,Any}` which maps the name of
the parameter to a probability function `p`.  This `p` takes the value
of the parameter and returns the log-prior.

`full_theta_priors` is a tuple of functions.  Each function takes a
`theta` vector and returns the log-prior for it.

If only theta is passed in then it uses `def_component_wise_priors()`
and `def_full_theta_priors()` as defaults.

Return:
- the function `logprior(theta)`
- a function `logprior_debug(theta)` which will print the values of the individual
  components of the prior and return them
"""
function make_logprior(theta0::Theta0,
                       component_wise_priors::Dict{Symbol}=def_component_wise_priors(),
                       full_theta_priors::Dict{Symbol}=def_full_theta_priors(theta0))
    cwps = Dict{Symbol,Any}()
    for (i,name) in enumerate(keys(theta0.names))
        if haskey(component_wise_priors, name)
            cwp = component_wise_priors[name]
            inds = theta0.names[name]
            f = _make_one_logprior(inds, cwp)
            cwps[name] = f
        end
    end
    len_theta0 = length(theta0.th0)

    # return function which adds all prior-functions
    logprior = function (theta::Vector)
        # make sure the length of theta is correct
        @assert length(theta)==len_theta0 "$(length(theta)) != $len_theta0"
        out = 0.0
        for fn in values(cwps)
            out += (fn(theta))::Float64
        end
        for fn in values(full_theta_priors)
            out += (fn(theta))::Float64
        end
        out
    end

    logprior_debug = function (theta::Vector, verbose=true)
        # make sure the length of theta is correct
        @assert length(theta)==len_theta0 "$(length(theta)) != $len_theta0"

        out = 0.0
        out2 = Float64[]
        verbose && println("** Component-wise priors:")
        for (k,fn) in cwps
            out += (fn(theta))::Float64
            push!(out2, fn(theta))
            verbose && println("$k: $(fn(theta))")
        end
        verbose && println("** Full priors:")
        for (k,fn) in full_theta_priors
            out += (fn(theta))::Float64
            push!(out2, fn(theta))
            verbose && println("$k: $(fn(theta))")
        end
        out, out2
    end

    return logprior, logprior_debug
end

function _make_one_logprior(theta_inds, prior::Function)
    function (theta::Vector)
        out = 0.0
        for i in theta_inds
            out += prior(theta[i])
        end
        out
    end
end

## Useful prior functions

"""
    logprior_scalar(val, expectation, sigma, range_abs, range_rel, range_rel_cutoff)

For a given value `val` and its `expectation`, return a
log-probability according to the sigma, range_abs, range_rel,
range_rel_cutoff (the fields of SomeData).  It assumes a truncated
Normal distribution.
"""
function logprior_scalar(val, expectation, sigma, range_abs, range_rel, range_rel_cutoff)
    # check that val is in abs-range and rel-range:
    range_abs[1]<=val<=range_abs[2] || return -Inf
    if abs(expectation)>range_rel_cutoff
        range_rel[1]<=val/expectation<=range_rel[2] || return -Inf
    end

    # Gaussian
    return -(val-expectation)^2/(2*sigma^2) - log(sigma)
end

"""
    logprior_scalar_subnormal(val, expectation, sigma, range_abs, range_rel, range_rel_cutoff)

For a given value `val` and its `expectation`, return a
log-probability according to the sigma, range_abs, range_rel,
range_rel_cutoff (the fields of SomeData).  If val is within the range,
then assume a Normal distribution, outside a super-Normal with power 10.

I.e. once outside the range the probability falls sharply.
"""
function logprior_scalar_subnormal(val, expectation, sigma, range_abs, range_rel=(-Inf,Inf), range_rel_cutoff=Inf)
    # Normal dist
    out = -(val-expectation)^2/(2*sigma^2) - log(sigma)

    # if val is outside abs-range or rel-range use sub-Normal, aka sub-Gaussian
    exponent = 10
    if val<range_abs[1]
        out += -(val-range_abs[1])^exponent
    elseif range_abs[2]<val
        out += -(val-range_abs[2])^exponent
    end

    @assert  range_rel_cutoff==Inf # not implemented

    return out
end

"""
    logprior_1dfields(theta, theta0::Theta0, name)

Priors for the parameters encoding 1d fields (aside-reminder: those
parameters are offsets).  Note that we do the priors on the parameters
and not the field in order to avoid correlations.

The values for the prior distribution come from the `SomeData`
structure of the gl.bdot, gl.dhdt, gl.fsl, gl.temp fields.

Somewhat inconsistently:
- the SomeData-ranges are checked for the value at all elevation bands
- the theta-values themselves are then assumed to have a Gaussian PDF
  with SomeData.sigma and mu=0 (as they are offsets).

This is a full_theta_priors.
"""
function logprior_1dfield(theta, theta0::Theta0, name)
    gb = theta0.gb
    gl = gb.gl
    field1d = make_fields1d(theta, theta0, name)
    field1d_expectation = getfield1d(gb, name)
    if name==:btilde # combo of bdot and dhdt
        # here we need to get
        sigma = sqrt(gl.bdot.sigma^2 +
                     gl.dhdt.sigma^2)
        @assert all(gl.bdot.range_rel==gl.dhdt.range_rel) "range_rel need to be equal for bdot and dhdt"
        @assert all(gl.bdot.range_rel_cutoff==gl.dhdt.range_rel_cutoff) "range_rel_cutoff need to be equal for bdot and dhdt"
        @unpack range_rel, range_rel_cutoff = gl.bdot

        range_abs = (gl.bdot.range_abs[1]+gl.dhdt.range_abs[1],
                     gl.bdot.range_abs[2]+gl.dhdt.range_abs[2])
    else
        @unpack sigma, range_abs, range_rel, range_rel_cutoff = getfield(gl, name)
    end

    # check that all are in abs-range and rel-range:
    @inbounds for (i,fd) in enumerate(field1d)
        range_abs[1]<=fd<=range_abs[2] || return -Inf
        if abs(field1d_expectation[i])>range_rel_cutoff
            range_rel[1]<=fd/field1d_expectation[i]<=range_rel[2] || return -Inf
        end
    end

    # Gaussian distribution: run only over actual parameters as their
    # value on the elevation bands will be correlated.
    inds = theta0.names[name]
    out = 0.0
    @inbounds @fastmath for i in inds
        thi = theta[i]
        out -= thi*thi
    end
    return out/(2*sigma^2)
end


"""
    calc_expected_terminus_flux(terminus_flux::Tuple, gb::Bands)
    calc_expected_terminus_flux(fkind::FluxKind._FluxKind, vals, gb::Bands)

Return expected flux at terminus in m^3/s.  NOTE, this value will still be scaled in
terminus_flux_logprior
"""
calc_expected_terminus_flux(terminus_flux::Tuple, gb::Bands) =
    calc_expected_terminus_flux(terminus_flux[1], terminus_flux[2].vals, gb::Bands)
function calc_expected_terminus_flux(fkind::FluxKind._FluxKind, vals, gb::Bands)
    tflux_expectation = if fkind==FluxKind.abs
        # absolute flux
        return vals
    elseif fkind==FluxKind.rel_bdot_teminus
        flux_n(gb, gb.bdot)[end] # flux at terminus due to expectation bdot
    elseif fkind==FluxKind.rel_bdot_max
        out = maximum(flux_n(gb, gb.bdot)) # maximum flux due to expectation bdot (== total accumulation rate)
        if out<=0
            # if all bdot is negative, i.e. no accumulation area, use the maximum apparent mass balance:
            maximum(flux_n(gb, gb.bdot-gb.dhdt))
        else
            out
        end
    elseif fkind==FluxKind.rel_btilde_teminus
        (flux_n(gb, gb.bdot-gb.dhdt)[end]) # expected flux at terminus
    else
        error("Unknown flux-kind: $fkind")
    end
    # this needs to be bigger than zero because sigma is relative... hack it!
    #@assert tflux_expectation>0 "Expectation value of terminus flux zero or negative: $tflux_expectation"
    if !(tflux_expectation>0)
        tflux_expectation = -1.0
    end
    return tflux_expectation
end

function _calc_flux(theta, theta0::Theta0)
    gb = theta0.gb
    # calculate 1d fields
    btilde = make_fields1d(theta, theta0, :btilde) # value at current MCMC step
    # flux
    return flux_n(gb, btilde)
end



"""
    terminus_flux_logprior(theta::Vector, theta0::Theta0)

A prior for the whole of btilde-1D to enforce that the total flux at
the terminus is according to specs.  See `Glacier` type-def.

Notes:
- this is a full_theta_prior
"""
function terminus_flux_logprior(theta::Vector, theta0::Theta0)
    gb = theta0.gb
    gl = gb.gl
    fkind, terminus_flux = gl.terminus_flux
    fkind==FluxKind.none && return 0.0
    @unpack vals, sigma, range_abs, range_rel, range_rel_cutoff = terminus_flux
    tflux = _calc_flux(theta, theta0)[end]
    # expected flux at terminus
    tflux_expectation = calc_expected_terminus_flux(fkind, vals, gb)

    # hack to catch glaciers with no accumulation area
    # TODO make this proper
    if fkind!=FluxKind.abs && tflux_expectation==-1
        tflux_expectation = 0.0
        # some heuristic:
        heuristic_max_flux = area(gl)/1e6/100
        sigma = heuristic_max_flux / 10
        fkind = FluxKind.abs
    end

    if fkind==FluxKind.abs
        return logprior_scalar_subnormal(tflux, tflux_expectation, sigma, range_abs, range_rel, range_rel_cutoff)
    elseif fkind==FluxKind.rel_bdot_teminus || fkind==FluxKind.rel_btilde_teminus || fkind==FluxKind.rel_bdot_max
        # all SomeData-stuff is relative here...
        @assert range_rel == (-Inf,Inf)
        @assert range_rel_cutoff == Inf
        # convert sigma etc to relative values
        return logprior_scalar_subnormal(tflux, tflux_expectation*vals, tflux_expectation.*sigma,
                                         tflux_expectation.*range_abs, range_rel, range_rel_cutoff)
    else
        error("FluxKind $fkind not implemented")
    end
end

"""
    continuous_glacier_logprior(theta::Vector, theta0::Theta0)

Prior to make sure that the glacier is not in two or more pieces,
i.e. that h=0 for some stretch and the h>0 again downstream.
Only takes 0 or -Inf as values.

Note: this will impact the allowed values of btilde only.
"""
function continuous_glacier_logprior(theta::Vector, theta0::Theta0)
    flux = _calc_flux(theta, theta0)
    @assert length(flux)>0
    ifirst = findfirst(flux.>0)
    ilast = findlast(flux.>0)
    if ifirst==0 || ilast==0
        # there is some odd bug I try to catch...
        println("\n========================================\n=== continuous_glacier_logprior bug ===")
        @show flux
        @show theta
        @show theta0
        println("\n========================================\n")
        return -Inf
    end
    return any(flux[ifirst:ilast].<=0) ? -Inf : 0.0
end

"""
    continuous_glacier_from_top_logprior(theta::Vector, theta0::Theta0)

Prior to make sure that the glacier starts at the top and is not in
two or more pieces, i.e. that h=0 for some stretch and the h>0 again
downstream.  Only takes 0 or -Inf as values.

Note: this will impact the allowed values of btilde only.
"""
function continuous_glacier_from_top_logprior(theta::Vector, theta0::Theta0)
    flux = _calc_flux(theta, theta0)
    @assert length(flux)>0
    ilast = findlast(flux.>0)
    ilast==0 && return -Inf
    @assert flux[1]==0
    return any(flux[2:ilast].<=0) ? -Inf : 0.0
end

"""
    def_component_wise_priors()

Default component-wise priors.  Probably mostly for MPara.
"""
function def_component_wise_priors()
    Dict(:dist_exp => x-> 0.1<x<=1 ? 0.0 : -Inf,
         :alpha_min => x -> 0<x ? 0.0 : -Inf,
         :iv_h_exp => x -> 0.1<x<5 ? 0.0 : -Inf,
         :iv_dist_exp1 => x -> 0<=x<=1 ? 0.0 : -Inf, # Huss-type
         :iv_dist_exp2 => x -> 0<x<20 ? -(x-4)^2/(2*0.5^2) : -Inf, # Nye-type; theory says n+1, effects for >20 are minimal, thus cap it.
         :sigma_h_model => x -> 0<x<200 ? 0.0 : -Inf, # -(x-20)^2/(2*10^2)
         :sigma_iv_model => x -> 0<x<200 ? 0.0 : -Inf, # -(x-10)^2/(2*10^2)
         )
end


"""
    def_full_theta_priors(theta0)

Default full-theta priors.  Includes all priors on the 1D fields and
 the prior on the flux at the terminus `terminus_flux_logprior`.
"""
function def_full_theta_priors(theta0)
    # extra priors, namely the terminus_flux and that the glacier is continuous:
    out = Dict{Symbol,Function}(
        :teminus_flux => theta -> terminus_flux_logprior(theta, theta0), # uses values of gl.terminus_flux
        :continuous_glacier => theta -> continuous_glacier_from_top_logprior(theta, theta0))
    # priors on the 1d-fields originating from the fields gl.bot, gl.dhdt, gl.temp, gl.fsl
    for n in [:btilde, :temp, :fsl]
        if n in keys(theta0.names)
            out[n] = theta -> logprior_1dfield(theta, theta0, n)
        end
    end
    return out
end
