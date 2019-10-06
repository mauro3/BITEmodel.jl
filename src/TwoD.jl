# The 2D model:
# - interpolate data onto flow band
# - pass it to 1D model
# - extrapolate data back to 2D


################
# Forward model 2D
################

# """
# The forward model needs to implement the function:

#     fwdm{P<:APara, M<:APara, N<:APara}(dem, b, fq, po::P, mp::M, nu::N) -> (h, us)

# Input:

# - 2D fields (or a scalar):
#   - dem: DEM of the glacier, includes a mask of the glacierized cells
#   - bdot: mass balance
#   - dhdt: dh/dt
#   - fq: fraction of deformational ice flux
#   - A: ice flow parameter
# - pp: misc. physical parameters: physical constants
# - mp: misc. model parameters: somewhat physical things, but quite ad-hoc
# - np: misc. numerical parameters: e.g. discretisation size, etc.

# Output:
# - h: ice thickness
# - us: surface flow speed

# Units: use SI (m,kg,s)
# """
# function fwdm end

# NOTE: do not access gb.fsl, gb.btilde, gb.temp, etc here!  Use
# values passed in as function arguments

"""
    fwdm(gl::Glacier, pp::Phys, pm::MPara, pn::Num)
    fwdm(gb, pp::Phys, pm::MPara, pn::Num,
              btilde = gb.bdot-gb.dhdt,
              fsl = gb.fsl,
              temp = gb.temp,
              )

The 2D forward model function does the following steps:
- interpolate onto 1D (this can be done separately with `make_bands!(gl,pp,pm,pn)`)
- run 1D model
- extrapolate onto 2D

Input:
- `fwdm(gl::Glacier, pp::Phys, pm::MPara, pn::Num)`
- `fwdm(gb, pp::Phys, pm::MPara, pn::Num[, btilde, fsl, temp])` does
  no 2D->1D and uses the 1D-fields if provided.  Otherwise uses
  gb.bdot, etc.


Parameters (to be fitted):
- A: ice flow factor (scalar)
- fq: ice sliding ratio (field)


Return:

    hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_local, vol_ratio, gb

Units: use SI (m,kg,s) except for flow velocity and SMB (.../a)

TODO: add method which takes the map_onto_band variables
"""
function fwdm(gl::Glacier, pp::Phys, pm::MPara, pn::Num)
    # map onto bands
    gb, pm = make_bands(gl, pp, pm, pn)
    return fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
end


function fwdm(gb, pp::Phys, pm::MPara, pn::Num,
              btilde, fsl, temp)

    gl = gb.gl
    # call fwdm1d
    hs1d_n, taus1d_local, taus1d_n, ivs1d_n, qtot1d_n = fwdm1d(gb, btilde, fsl, temp, pp, pm, pn)
    # extrapolate back to 2D
    hs2d, taus2d, ivs2d, ivs2d_band, ivs2d_band_mask, vol, vol_above =
        extrapolate(hs1d_n, taus1d_n, ivs1d_n, btilde, fsl, gl, gb, pp, pm)
    vol_ratio = check_vol_1D_2D(gb, hs1d_n, vol)
    if pn.verbose && !(0.9<vol_ratio<1.1)
        warn("Volume of 1D and 2D model is more than 10% different.  Ratio 1D/2D=$vol_ratio")
    end
    return FWDSol(hs1d_n, ivs1d_n, taus1d_n, taus1d_local,
                  hs2d, ivs2d, ivs2d_band, ivs2d_band_mask, taus2d,
                  vol, vol_above, vol_ratio,
                  gb, pp, pm, pn, btilde, fsl, temp, Dict{Symbol,Any}())
end


################
# 2D->1D
################
#
# The 2D data gets binned into elevation bands.

"""

Holds elevation bands information for a glacier, which are geometric
properties and 1D fields (their prior expectation value in a Bayesian
context).  Band indices start at the top of the glacier.

Also holds entities to allow extrapolation from bands back to 2D.

Note that gb.bands are the bands starting from the top and their location corresponds to the band-start.
gb.ele gives the band midpoints, gb.bands_ is like gb.bands but includes also the end-point (terminus).

"""
@with_kw struct Bands{Gl<:Glacier}
    gl::Gl
    dem::Gridded{F} # possibly smoothed surface DEM [m] (if pm.window_dem_smooth>0) otherwise original
    bands::SRangeL{Float64} # elevation bands starts.  Runs from top (divide) to bottom (terminus)
    bands_::SRangeL{Float64} # elevation bands starts and end
    bandi::Vector{Vector{Int}} # indices into gl.dem.v for each band
    ele::SRangeL{Float64} # mean elevation of bands, i.e. bands mid-points
    ws::Vector{F}       # width
    ls::Vector{F}       # length
    malphas::Vector{F}  # mean slope of a band
    x::Vector{F}        # distance from divide to beginning of bands (i.e. to nodes)
    xmid::Vector{F}     # distance from divide to middle of bands (i.e. cell-centers)

    # Smoothed 1D fields which represent the expectation value of
    # their prior distribution.
    # TODO: figure out how to band the sigmas of gl.bdot, etc.
    bdot::Vector{F} # [m ice/a]
    dhdt::Vector{F} # [m ice/a]
    fsl::Vector{F} # []
    temp::Vector{F} # [C]
    # original fields (not smoothed)
    fields1d_original::Dict{Symbol,Vector{F}}

    # 1D->2D extrapolation:
    # weights used in extrapolation (TODO put elsewhere?)
    weights::Vector{Vector{F}}
    # precualculated boxcar filter matrix
    boxcar_M::Tuple{SparseMatrixCSC{F,Int},F}
    boxcar_M_iv::Tuple{SparseMatrixCSC{F,Int},F,Vector{Dict{Int,Vector{VAWTools.Line}}},
                       Matrix{F}, Matrix{F}}
    # 2D surface slope calculated with a possibly smoothed dem (the dem in above field)
    alpha2d::Matrix{F}
end
# function Base.show(io::IO, gb::Bands)
#     print(io, "Elevation bands discretisation of:\n  ")
#     Base.show(io, gb.gl)
# end
Base.length(gb::Bands) = length(gb.bands)

"Return a 1D field, including the non-field btilde"
function getfield1d(gb, name)
    if name==:btilde
        return apparent_mass_balance(gb.bdot, gb.dhdt)
    else
        return getfield(gb, name)
    end
end


"""
$(SIGNATURES)

Calculate bands.  This is one of the core functions!

In:
- glacier
- pp
- MPara: model parameters

Out:
- Bands
"""
function make_bands(gl::Glacier, pp::Phys, pm::MPara, pn::Num)
    @unpack dem, glaciermask = gl
    @unpack bandsize, window_width_smooth, window_dem_smooth,
            slope_min, slope_max, iv_window_frac, window = pm
    @unpack calc_boxcar_M = pn

    @assert bandsize>0 "Need a pm.bandsize>0"

    bands, bandi, malphas, areas, lengths, widths, x, xmid, dem_smoothed, alpha2d_smoothed =
        make_1Dglacier(dem, -bandsize, glaciermask;
                       window_dem_smooth=window_dem_smooth,
                       window_width_smooth=window_width_smooth,
                       alpha_min=slope_min, alpha_max=slope_max,
                       FILL=FILL,
                       min_bands=3,
                       min_bin_number_ends=5,
                       verbose=pn.verbose
                       )

    ele = bands+step(bands)/2 # mean elevation of each band
    pm = MPara(pm, bandsize=-step(bands))

    # h-extrapolation weights `sin(alpha)^(-n/(n+2))`
    weights = extrapolation_weigths(bands, bandi, widths, lengths,
                                       pm.alpha_min, alpha2d_smoothed, pp.n)

    # thickness smoothing operator
    dx = step(dem.x)
    if calc_boxcar_M
        boxcar_M = boxcar_matrix(F, round(Int,window/dx), gl.glaciermask, (!).(gl.glaciermask))
        window_boxcar_M = window
        # IV smoothing operator for :masscons method
        boxcar_M_iv, boundaries, ux, uy = VAWTools.get_iv_boxcar_M(F, dem_smoothed, glaciermask, bands, bandi, lengths, iv_window_frac)
        window_boxcar_M_iv = iv_window_frac
    else # set to dummies which will trigger the use of the straight
         # boxcar in extrapolate
        boxcar_M = speye(0)
        window_boxcar_M = -9999.0
        boxcar_M_iv = speye(0)
        window_boxcar_M_iv = -9999.0
        boundaries = Dict{Int,Vector{VAWTools.Line}}[]
        ux = zeros(F,0,0)
        uy = zeros(F,0,0)
    end


    # 1d fields: see what is in gl and do the appropriate thing.
    # Steps:
    # - if 2D: interpolate to 1D & smooth
    # - if 1D: interpolate to elevation bands & smooth
    # - if scalar: make constant for all elevation bands

    fields1d = Dict{Symbol,Vector{F}}()
    fields1d_original = Dict{Symbol,Vector{F}}()
    for name in [:iv, :bdot, :dhdt, :fsl, :temp]
        glfld = getfield(gl, name)
        if isa(glfld, SomeData{Gridded{F}}) # 2D
            gr = glfld.vals
            bandi_ = VAWTools.bandi_for_other_grid(bands, bandi, dem_smoothed, gr, gr.v.!=FILL)
            vals = map_onto_bands(bandi_, gr.v, mean, FILL)
            vals = [(isnan(v) || v==FILL) ? glfld.vals_prior : v for v in vals]
            for v in vals
                if name!=:iv && (isnan(v) || v==FILL)
                    @show name, v, glfld.vals_prior
                end
            end
            fields1d_original[name] = vals
            # smooth: This does not conserve mass!  TODO is this good?
            fields1d[name] = smooth_vector(ele, vals, glfld.smooth_window)
        elseif isa(glfld, SomeData{Gridded1d{F}}) # 1D
            # interpolate onto elevation bands
            if glfld.vals.x==ele # no interpolation
                fields1d_original[name] = glfld.vals.v
                # smooth: This does not conserve mass!  TODO is this good?
                fields1d[name] = smooth_vector(ele, glfld.vals.v, glfld.smooth_window)
            else
                # linear interpolation
                f = piecewiselinear(glfld.vals.x, glfld.vals.v)
                fields1d_original[name] = f.(ele)
                # smooth: This does not conserve mass!  TODO is this good?
                fields1d[name] = smooth_vector(glfld.vals.x, glfld.vals.v,
                                               glfld.smooth_window, ele)
            end
        elseif isa(glfld, SomeData{F}) # 0D
            nb = length(bands)
            fields1d_original[name] = zeros(F,nb) + glfld.vals
            fields1d[name] = zeros(F,nb) + glfld.vals
        elseif isa(glfld, SomeData{Traj{F}})
            fields1d_original[name] = glfld.vals.v
            fields1d[name] = F[] # TODO
        else
            if name==:iv && isa(glfld, SomeData{Void})
                fields1d_original[name] = F[]
                fields1d[name] = F[]
            else
                error("Need to provide a value for field $name in Glacier datastruct.")
            end
        end
        if name!=:iv
            any(isnan.(fields1d_original[name])) && error("1D field $name has NaN: $(fields1d_original[name])")
            any(isnan.(fields1d[name])) && error("Smoothed 1D field $name has NaN: $(fields1d_original[name])")
        end
    end

    return Bands{typeof(gl)}(gl, dem_smoothed,
                             bands, first(bands):step(bands):last(bands)+step(bands), bandi,
                             ele, widths, lengths, malphas, x, xmid,
                             fields1d[:bdot], fields1d[:dhdt], fields1d[:fsl], fields1d[:temp],
                             fields1d_original,
                             weights, (boxcar_M, window_boxcar_M),
                             (boxcar_M_iv, window_boxcar_M_iv, boundaries, ux, uy),
                             alpha2d_smoothed), pm::MPara
end

VAWTools.map_onto_bands(gb::Bands, field, fn=mean) = map_onto_bands(gb.bandi, field, fn, FILL)
VAWTools.bandi_for_other_grid(gb::Bands, othergrid::Gridded, othermask=trues(size(othergrid.v)) ) =
    VAWTools.bandi_for_other_grid(gb.bands, gb.bandi, gb.dem_smoothed, othergrid, othermask)


"""
    map_iv_onto_bands(gl::Glacier, gb::Bands)

Maps the IV in `gl.iv` (which is usually on a different grid to
`gl.dem`) onto the elevation bands.

Does not take mass-flux into account.

Not in use.
"""
function map_iv_onto_bands(gl::Glacier, gb::Bands)
    bandi = VAWTools.bandi_for_other_grid(gb, gl.iv.vals, gl.ivmask)
    return map_onto_bands(bandi, gl.iv.vals.v, mean, FILL)
end



"""
    extrapolation_weigths(bands, bandi, ws, ls,
                               alpha_min, alpha2d, n)

Pre-calculates the weights used in the extrapolation due to local
slope.
"""
function extrapolation_weigths(bands, bandi, ws, ls,
                               alpha_min, alpha2d, n)
    weights = Vector{F}[]
    for (i,band) in enumerate(bands)  # @inbounds does not help
        w, l = ws[i], ls[i]
        inds = bandi[i]
        weight = similar(inds, F)
        for (j,jj) in enumerate(inds)
            alpha_ = max(alpha_min, alpha2d[jj])
            weight[j] = sin(alpha_)^( -n/(n+2))
        end
        push!(weights, weight)
    end
    return weights
end


# """
# Values given on bands.
# """
# struct BandValues{N<:Number}
#     gb::Bands
#     v::Vector{N}
# end
# BandValues{T}(gb::Bands, field::Gridded{T}) = BandValues{T}(gb, map_onto_bands(b, field, mean, FILL))

################
# 1D->2D
################
#

"""
$(SIGNATURES)

    -> hs2d, taus2d, ivs2d, ivs2d_band, ivs2d_band_mask

Extrapolation of the output of the 1D model to 2D.

TODO:
- use directional/an-isotropic smoothing
"""
function extrapolate(hs_n, taus_n, ivs_n, btilde, fsl, gl::Glacier, gb::Bands,
                               pp::Phys, pm::MPara)
    @unpack n, rho, rhosw, g, r = pp
    @unpack alpha_min, window, cap_at_floatation,
            dist_land_maxdist, dist_exp,
            iv_extrapolation_scaling, iv_tau_exp, iv_h_exp, iv_dist_exp1,
            iv_dist_exp2, iv_nye,
            iv_flux_dir_window, iv_window_frac, iv_scale_factor = pm
    @unpack dist_land, dist_land_maxdist_calc, glaciermask = gl
    @unpack weights, ws, ls, bands, bandi, boxcar_M, boxcar_M_iv = gb
    @assert dist_land_maxdist_calc>=dist_land_maxdist "Calculated distance in gl.dist_land is smaller than pm.dist_land_maxdist."

    dx = step(gl.dem.x)
    hs2d = zeros(gl.dem.v)
    taus2d = zeros(gl.dem.v)  # effective driving stress
    ivs2d = zeros(gl.dem.v)

    hs = on_cells(hs_n)
    taus = on_cells(taus_n)
    ivs = on_cells(ivs_n)

    ## Local ice thickness weighted proportional to
    # - (distance to land)^dist_exp
    # - (sin α)^(-n/(n+2))
    for (i,band) in enumerate(bands)  # @inbounds does not help
        h, w, l = hs[i], ws[i], ls[i]
        vol = w*l*h
        inds = bandi[i]
        weight = weights[i] # pre-calculated weight due to surface slope
        vol2d = 0.0
        for (j,jj) in enumerate(inds)
            hs2d[jj] = weight[j] * min(dist_land[jj], dist_land_maxdist)^dist_exp
            vol2d += hs2d[jj]*dx^2
        end
        # scale to conserve total volume of band (note this is non-linear in h)
        for jj in inds
            hs2d[jj] *= vol/vol2d
        end
    end

    # Smooth.  Does not quite conserve mass (because of edge effects),
    # also increases thickness at the edge!
    #
    # TODO this takes a lot of time, even with the sparse multiply;
    # in particular for larger smoothing windows.
    # Consider: for the fitting, one could only smooth at the radar points.
    if boxcar_M[2]==window
        hs2d = apply_boxcar_matrix(boxcar_M[1], hs2d)
    else
        # hs2d = boxcar(hs2d, round(Int,window/dx), ones(size(hs2d)), !gl.glaciermask)
        hs2d = boxcar(hs2d, round(Int,window/dx), gl.glaciermask, .!gl.glaciermask)
    end

    # Remove bits above flotation (does not conserve mass) and calculate
    # volume above flotation.
    dA = step(gl.dem.x) * step(gl.dem.y)
    vol, vol_above = 0.0, 0.0
    fac = 1/((rhosw-rho)/rhosw)
    for j=1:length(hs2d)
        !(glaciermask[j]) && continue
        zs = gl.dem.v[j]
        h = hs2d[j]
        hmax = zs*fac
        if cap_at_floatation && h>hmax
            hs2d[j] = hmax
            h = hmax
        end
        vol += h*dA
        # thickness such that just grounded
        hmin = max(0.0, - (zs - h) * rhosw/rho)
        vol_above += (h-hmin)*dA
    end

    ## Tau
    # Clarke et al 2013 assume that tau is the same over an elevation
    # band, we do this too:
    map_back_to_2D!(taus2d, bandi, taus)
    ## Ice velocity
    for (i,band) in enumerate(bands)  # @inbounds does not help
        # First make a trial u
        inds = bandi[i]
        for jj in inds
            d = min(dist_land[jj], dist_land_maxdist)

            len_scale = min(ws[i]/2,hs[i])
            d2 = min(dist_land[jj], min(len_scale, dist_land_maxdist) )
            # u-trial made with extrapolation, which can happen in a few different ways.
            # It is a trial-solution in the sense that it has not the correct scale
            # but the desired cross-valley/closest-to-outcrop shape.
            ivs2d[jj] =
                 !iv_nye * (
                    taus[i]^iv_tau_exp * # surface slope term on band
                    hs2d[jj]^iv_h_exp *  # ice thickness term
                    d^iv_dist_exp1 ) + # distance to land term 1
                  iv_nye * (
                      1 - ( (len_scale-d2)/len_scale )^iv_dist_exp2 ) # Nye 1965 inspired, Eq 6
        end

        # Scale the trial such that it "resembles" 1D field.
        if iv_extrapolation_scaling==:huss
            # The makes the 2D mean equal to the 1D flow.  Not mass conserving.
            mean_ = mean(ivs2d[inds])
            for jj in inds
                ivs2d[jj] = ivs[i]*ivs2d[jj]/mean_
            end
        elseif iv_extrapolation_scaling==:masscons
            # The makes the 2D flow such that the flux from one
            # elevation band to the next is correct.
            #
            # Whole step is performed below
        elseif iv_extrapolation_scaling==:none
            # Just scale by a constant.  Probably bad.
            for jj in inds
                ivs2d[jj] *= iv_scale_factor # a scale factor.  Only needed if no scaling is applied later.
            end
        else
            error("Unsupported value of iv_extrapolation_scaling = $iv_extrapolation_scaling")
        end
    end
    # Smooth IV
    if iv_extrapolation_scaling==:huss || iv_extrapolation_scaling==:none
        if boxcar_M[2]==window
            ivs2d = apply_boxcar_matrix(boxcar_M[1], ivs2d)
        else
            ivs2d = boxcar(ivs2d, round(Int,2*window/dx), gl.glaciermask, .!gl.glaciermask)
        end
    elseif iv_extrapolation_scaling==:masscons
        # This is both scaling and smoothing in one
        scale_surface_iv = 1./((1-r)*fsl + r)
        if boxcar_M_iv[2]==iv_window_frac # use the pre-calculated smoothing matrix
            boundaries = boxcar_M_iv[3]
            bM = boxcar_M_iv[1]
            ux, uy = boxcar_M_iv[4], boxcar_M_iv[5]
            ivs2d, ivs2d_band, ivs2d_band_mask = VAWTools.calc_u(flux_n(gb, btilde)*year2sec,
                                                     boundaries, ivs2d, hs2d,
                                                     ux, uy, dx, glaciermask, bands, bandi, ls,
                                                     iv_flux_dir_window,
                                                     bM, scale_surface_iv)
        else
            ux,uy = VAWTools.get_ux_uy(gb.dem, gl.glaciermask)
            binmat = VAWTools.bins2matrix(gl.dem, bands, bandi)
            boundaries = VAWTools.calc_boundaries(bands,bandi,binmat)
            ivs2d, ivs2d_band, ivs2d_band_mask = VAWTools.calc_u(flux_n(gb, btilde)*year2sec, boundaries, ivs2d, hs2d,
                                                                 ux, uy, dx, glaciermask, bands, bandi, ls,
                                                                 iv_flux_dir_window,
                                                                 iv_window_frac, scale_surface_iv)
        end
    else
        error("Unsupported value of iv_extrapolation_scaling = $iv_extrapolation_scaling")
    end

    return hs2d, taus2d, ivs2d, ivs2d_band, ivs2d_band_mask, vol, vol_above
end

"""
Forward model solution

```
@with_kw struct FWDSol
    # Fields from the 1D calculation
    hs1d::Vector{F}
    ivs1d::Vector{F}
    taus1d::Vector{F}
    taus1d_local::Vector{F}

    # Fields from the 2D extrapolation
    hs2d::Matrix{F}
    ivs2d::Matrix{F} # the smoothed ivs2d_band
    ivs2d_band::Matrix{F} # the calculated, mass conserving IV across bands
    ivs2d_band_mask::BitMatrix
    taus2d::Matrix{F}

    vol::Float64 # total ice volume (calculated from 2D field)
    vol_above::Float64 # ice volume above flotation
    vol_ratio::Float64 # ration between ice volume calculated from 1D and 2D field.
                       # Gives an indication of how much volume as not conserved in the extrapolation.

    # model input
    gb::Bands
    pp::Phys
    pm::MPara
    pn::Num
    btilde::Vector{F}
    fsl::Vector{F}
    temp::Vector{F}

    # to add various other stuff
    misc::Dict{Symbol,Any}
end
```
"""
@with_kw struct FWDSol
    # Fields from the 1D calculation
    hs1d::Vector{F}
    ivs1d::Vector{F}
    taus1d::Vector{F}
    taus1d_local::Vector{F}

    # Fields from the 2D extrapolation
    hs2d::Matrix{F}
    ivs2d::Matrix{F} # the smoothed ivs2d_band
    ivs2d_band::Matrix{F} # the calculated, mass conserving IV across bands
    ivs2d_band_mask::BitMatrix
    taus2d::Matrix{F}

    vol::Float64 # total ice volume (calculated from 2D field)
    vol_above::Float64 # ice volume above flotation
    vol_ratio::Float64 # ration between ice volume calculated from 1D and 2D field.
                       # Gives an indication of how much volume as not conserved in the extrapolation.

    # model input
    gb::Bands
    pp::Phys
    pm::MPara
    pn::Num
    btilde::Vector{F}
    fsl::Vector{F}
    temp::Vector{F}

    # to add various other stuff
    misc::Dict{Symbol,Any}
end

volume(gb::Bands, hs1d_n) = sum(on_cells(hs1d_n).*gb.ws.*gb.ls)
volume(hs2d, dx) = sum(hs2d)*dx^2
volume(fwdsol) = volume(fwdsol.hs2d, step(fwdsol.gb.gl.dem.x))

##########
# Checks
##########
"Check volume of 1D vs 2D.  Returns mass1D/mass2D"
check_vol_1D_2D(gb, hs1d_n, vol) = volume(gb, hs1d_n)/vol
