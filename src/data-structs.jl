using Parameters

##############################
# Parameters
##############################

# conversion factor
const day2sec = 24*60*60
const year2sec = 365*day2sec

# "Convert m^3/s in Gt/a"
# qms2Gta(flux, rho) = flux*rho/1e3/1e9*year2sec

# A ice flow factor lookup table (fixed):
# Temperature dependence of A (Table 3.4, p.75), for n=3 [Pa^-3 s^-1]
const As = (0.0:-5:-50,[ 2.4e-24
                         9.3e-25
                         3.5e-25
                         2.1e-25
                         1.2e-25
                         6.8e-26
                         3.7e-26
                         2.0e-26
                         1.0e-26
                         5.2e-27
                         2.6e-27])
const _A_itr = Interp.interpolate((reverse(As[1]),), reverse(As[2]), Interp.Gridded(Interp.Linear()))
"""
Ice flow rate factor A as function of temperature.  Only valid for n=3

Interpolated from a lookup table (C&P p. 75)

Input:
- temp: ice temp
- pp:: physical para

Out:
- A
"""
ice_flow_factor(temp) = _A_itr[temp]

"""
$(TYPEDEF)
$(FIELDS)

Fairly fixed physical parameters, not used in tuning.

Units:
- kg, m, s, Celsius
"""
@with_kw struct Phys{FN<:Function} <:APara @deftype F
    # general
    g=9.81 # [m/s^2]
    # ice (references are to Cuffey&Patterson 2010)
    rhow=1000      # water density
    rhoi=910       # ice density
    rho=870        # glacier density including firn layer (kg m-3)
    rhosw=1025     # sea water density (kg m-3)
    n=3            # Glen's n
    @assert n==3   # needs to be ==3! for the A-table to be valid
    A::FN = ice_flow_factor
    r=(n+1)/(n+2) # (average deformational speed) / (surface deformational speed) for shallow ice
    @assert r==(n+1)/(n+2)
end

"""
What error is fitted:
- sqrerr: corresponds to Normal distribution of errors
- abserr: corresponds to Laplac distribution of errors
"""
baremodule ErrorModel
using Base: @enum, getindex
@enum(_ErrorModel, sqrerr, abserr)
end


"""
$(TYPEDEF)
$(FIELDS)

Model parameters.  These are not particularly fixed, ranging from
educated to wild guesses.  Some may be used in the inversion, some define the
inversion.

TODO:
"""
@with_kw struct MPara<:APara @deftype F
    # Standard deviation of errors due forward model limitations
    # (assuming perfectly tuned parameters):
    sigma_h_model = 10 # of ice thickness
    sigma_iv_model = 10 # of surface velocity
    @assert sigma_h_model>=0
    @assert sigma_iv_model>=0
    # total sigma is: sqrt(pm.sigma_h^2 + gl.h.sigma^2)

    # whether to assume Normal-distributed errors (:sqrerr)
    # or Laplace-distributed errors (:abserr)
    error_model::ErrorModel._ErrorModel=[ErrorModel.sqrerr, ErrorModel.abserr][1]

    # poor man's avoidance of spatial correlation: only use every n'th datapoint:
    nthin_h::Int=1
    nthin_iv::Int=10

    ## 2D->1D
    #########

    # Smooth dem to get smooth alpha, smoother bands.  This is in
    # particular important when there is cross-flow bumpiness, such as
    # on Unteraar.  However, it can also be bad.  YMMV, check!
    # This is then also used for h and iv extrapolation. The window
    # should probably be chosen around the size of the offending surface
    # features, or zero.
    window_dem_smooth = 100

    # Allowed minimal/maximal slope of bands
    slope_min = deg2rad(0.4)
    slope_max = deg2rad(60.0)

    bandsize = 30 # (m) arguably this should go into Num

    # Smooth the width of the bands over +/- this elevation range (m)
    # This makes the 1D iv much smoother, however it can make the
    # malphas much noisier!  In particular it can introduce bad
    # artifacts at the edges!
    window_width_smooth = 50

    ## 1D
    #####
    ## Scalars
    # Distance d over which to average tau in units of local ice thickness.
    # Takes average over elevation bands within +/- mean_tau_dist*h; no distance weighting is done.
    mean_tau_dist = 3.0  # Kamb and Echelmeyer 1986 (Fig 2 & Conclusions) suggest a triangularly-weighted
                         # window of length ±2*(1.5 to 2.5 ice thicknesses)

    # use Nye's ice-flow shape factor, or not
    shapeF::Bool = true

    # Whether to use the magic flux corrector
    use_flux_corrector::Bool=false

    # The ice thickness at the top of the glacier.  If set to
    # negative, set ice thickness at top equal to ice thickness one
    # below (good for ice-caps).
    # If set to NaN it will use -1 if isicecap==true and 0 otherwise.
    h_at_top = NaN

    ## 1D->2D
    #########

    # smoothing bed-dem with boxcar filter, probably around mean ice thickness:
    window = 500 # [m]

    # threshold on slopes used in 2D extrapolation
    alpha_min = deg2rad(1.5)

    # max distance to land-margin to consider
    dist_land_maxdist = 3e3 # [m] in Matthias' code: crit_distnu
    dist_exp = 0.3 # exponent to use in the distance to land-margin.  <1 gives more U-shape
    @assert dist_exp>0

    # Cap ice thickness at flotation thickness.  (Setting to true will violate mass-cons)
    cap_at_floatation::Bool = false

    ## IV extrapolation
    #
    # Trial u is made with the formula:
    #
    # !iv_nye *(
    # taus[i]^iv_tau_exp * # surface slope term on band
    # hs2d[jj]^iv_h_exp *  # ice thickness term
    # d[jj]^iv_dist_exp1) + # distance to land term 1
    # iv_nye * (1 - ( (len_scale-d2)/len_scale )^iv_dist_exp2 ) # Nye 1965 inspired
    #
    # Shallow ice is: tau^n*h^(n+1)
    iv_tau_exp = [0, 3][2] # no-effect at 0
    @assert iv_tau_exp>=0
    iv_h_exp = [0, 4][2] # no-effect at 0
    @assert iv_h_exp>=0
    iv_dist_exp1 = [0, 0.3][2] # no-effect at 0
    @assert 0<=iv_dist_exp1<=1

    # Nye 1965 extrapolation:
    iv_nye::Bool = true # to switch-on Nye 1965-extrapolation (and switch off Huss-extrapolation)
    iv_dist_exp2 = [3+1, 1][1] # larger makes IV more even flow-field across glacier, n+1 according to theory.
                               # Ceases to have much effect when >10.
    @assert iv_dist_exp2>0

    # How the trial u is scaled.  Best :masscons
    iv_extrapolation_scaling::Symbol = [:huss, :masscons, :none][2]
    # :huss has no parameters

    # :masscons parameters
    # Over how many grid cells to average the flux direction.  Does not
    # seem to make much difference (at least when using a smoothed DEM).
    iv_flux_dir_window = 300 # [m]
    iv_window_frac = 0.7 # smoothing distance as fraction of max band width

    # :none parameters
    iv_scale_factor = 1 # a scale factor.  Only needed if :none is used.

    # Afl_all::Vector{F}=[0.0255,0.0175,0.0405,0.0485,0.05,0.025,0.014,0.0255,0.0255]               # flow rate factor


    # # cs parabelflörmig über den Gletscher verteilt
    # cs=0.              # minimaler Wert der Korrektur-Funktion
    # rancurv::Vector{F}=[4,6]      # Exponent

    # # Options for redistribution of ice thickness within elevation band
    # fact_ip=0.75
    # #fact_ip=0.5           # exponent regulating distribution of thickness values in the cross section
    # # (0<fact<1) low: U-valley, high: V-Valley
    # slp_crit=1.5             # Cut off threshold for slope influence [deg]

    # rho_glacier=870.     # kg m-3 (glacier density including firn layer)

    # # surface elevation changes
    # dh_dtelev_boundary=30    # % lowermost x %
    # dh_dtsurf::Vector{F}=[-1,-7.5]       # m/a (buttressed / not buttressed)
    # #dh_dtsurf::Vector{F}=[-0.25,-3]       # m/a (buttressed / not buttressed)
    # #dh_dtsurf::Vector{F}=[-2,-15]       # m/a (buttressed / not buttressed)

    # meanslope_lower=.4      # [deg] Lower boundary for slope in elevation band ...
    # coldice_temp=-8
    # fact_a=6.
    # crit_distnu=3000 # [m] critical distance from nunataks - no influence when further away!


end

"""
$(TYPEDEF)
$(FIELDS)

Numerical parameters for deterministic model.
"""
@with_kw mutable struct Num<:APara @deftype F
#    bstep=10       # elevation band discretisation
    reltol=0.05    # relative tolerance for iterative solver
    minthick=1.0    # minimal ice thickness, acts as a abstol
    niter::Int=100 # max number of iterations in Picard loop
    relaxfac=0.1   # relax the Picard iteration: 0 -> no relaxation, 0.1 recommended
    relaxmin=0.1  # relax the Picard iteration: 0 -> no relaxation, 0.05  recommended (TODO?)
    @assert relaxmin>=2*reltol
    # performance
    calc_boxcar_M::Bool=true # set to false when only doing a few forward runs.
    # output
    plotyes::Bool = false
    verbose::Bool = false
    # misc flags
    test::Bool=true # activate some tests
    error_on_neg_qtot::Bool = false # if false return h=0 where qtot<0, otherwise error
end


##############################
# Data
##############################

const SomeDataValType = Union{F,AGridded{F},Traj{F},Void}
"""
$(TYPEDEF)

Holds measurements or synthetic guesses (scalar, 1D-trajectory,
1D-elevation band data, 2D or missing) and their errors expressed as
standard-deviation and ranges (once collapsed into 1D elevation
bands).

See function `logprior_scalar` and `logprior_1dfields` for implementations.

Fields:

    vals::Union{F,AGridded{F},Traj{F},Void}}
    sigma::F # Standard deviation of error
    range_rel::Tuple{F,F} # Range of allowed values relative to local value
    range_rel_cutoff::F # only check range_rel for local values `lv` with `abs(lv)>range_rel_cutoff`
    range_abs::Tuple{F,F} # Range of allowed values in absolute values
    smooth_window::F # elevation range to smooth its 1D version over
                       (not applicable to scalars).  This essentially gives
                       a smoothness constraint on the data
    vals_prior::Union{F,AGridded{F},Traj{F},Void}} # use this value if the vals produce FILL or NaN

TODO:
- add a field which give the time/date range of the measurements.
- think about how to transform a, say, sigma which applies to 2D data,
  into a sigma which applies to the corresponding elevation band data.
- allow non-scalar sigma, etc.
"""
@with_kw struct SomeData{V<:SomeDataValType, VV<:SomeDataValType}
    vals::V  # mean/expectation value
    sigma::F=Inf # Standard deviation of vals distribution (Gaussian)
    range_rel::Tuple{F,F}=(-Inf,Inf) # Range of allowed values relative to local expectation value
    range_rel_cutoff::F=Inf # only check range_rel for local values `lv` with `abs(lv)>range_rel_cutoff`
    range_abs::Tuple{F,F}=(-Inf,Inf) # Range of allowed values in absolute values
    smooth_window::F=0 # elevation range to smooth its 1D version over
    vals_prior::VV=F(NaN) # use this value if the vals produce FILL or NaN
end
SomeData(val::V,sigma,range_rel,range_rel_cutoff,range_abs,smooth_window, vals_prior::VV) where {V,VV} =
    SomeData{V,VV}(val,sigma,range_rel,range_rel_cutoff,range_abs,smooth_window, vals_prior)
SomeData(val::V) where {V} = SomeData{V,F}(val, Inf, (-Inf,Inf), Inf, (-Inf,Inf), 0.0, NaN)
SomeData(val::V,sigma,range_rel,range_rel_cutoff,range_abs) where {V<:Union{F,Void}} =
    SomeData{V,F}(val,sigma,range_rel,range_rel_cutoff,range_abs,0.0,NaN)
hasdata(::SomeData{Void}) = false
hasdata(::SomeData{F}) = true
hasdata(sd::SomeData) = length(sd.vals.v)>0

"""
Regions of the world according to RGI and some extras (with numbers >99)

See also: symbol2region, region2symbol
"""
baremodule Region
using Base: @enum, getindex
@enum(_Region,
      Alaska=1,
      WesternCanadaUSA,
      ArcticCanadaNorth,
      ArcticCanadaSouth,
      Greenland,
      Iceland,
      Svalbard, # and Jan Mayen
      Scandinavia,
      RussianArctic,
      NorthAsia,
      CentralEurope,
      CaucasusMiddleEast,
      CentralAsia,
      SouthAsiaWest,
      SouthAsiaEast,
      LowLatitudes,
      SouthernAndes,
      NewZealand,
      SubAntarctic,
      # Misc
      NotKnown=100,
      Synthetic=101,
      # Other regions of the world
      AntarcticPeninsula=110,
      # ITMIX
      ITMIX=201,
      ITMIX2=202,
      )
end
let
    sym = Symbol.(instances(Region._Region))
    int = Int.(instances(Region._Region))
    global symbol2region(sym) = Region._Region(int[findfirst(x->x==sym, Region.sym)])
    global region2symbol(reg) = Symbol(reg)
end


# """
# The different kinds of surfaces we need to consider:

# - `this_glacier` -- the glacier under study
# - `other_glacier` -- where it borders on other glaciers
# - `land_margin` -- where it borders on solid land
# - `calving_margin` -- where it borders on something it can calf into
#                       (water for tide-water glaciers, air for hanging glaciers)
# """
# @enum SurfaceTypes this_glacier other_glacier_margin land_margin calving_margin

"""
$(TYPEDEF)

Glacier ID, type to dispatch on a particular glacier.  As
type-parameter use a symbol of the RGI, GLIMS id or some other ID.
"""
abstract type GlacierID end
Base.:(==)(a::GlacierID, b::GlacierID) = getid(a)==getid(b)

"Return glacier ID, usually a symbol or Int"
getid(gid::GlacierID) = gid.id
"Return glacier RGI number, a string"
getrgi(::GlacierID)::String = error("to be defined")
"Return the name of a glacier, if it has one (string)"
getname(gid::GlacierID)::String = gid.name
"Return the region (an instance of the `Region` enumeration)"
getregion(gid::GlacierID)::Region = gid.region


"""
$(TYPEDEF)

Glacier outline and outcrops.

The outline gives the glacier outline, usually including outcrops.
However, outcrops can also be specified with separate polygons.

The polygons follow the "bigpoly" convention as in
https://github.com/mauro3/VAWTools.jl/blob/67db8fe62daabfc10d313a46e827349644738ca1/src/VAWTools.jl#L607

    gid::GID
    outline::Matrix{Float64} # size==(2,n)
    splits::Vector{Int} # i-th poly has indices splits[i]:splits[i+1]-1
    outcrops::Matrix{Float64}
    outcrops_splits::Vector{Int}

"""
struct Outline{GID<:GlacierID}
    gid::GID
    outline::Matrix{Float64} # size==(2,n)
    splits::Vector{Int}
    outcrops::Matrix{Float64}
    outcrops_splits::Vector{Int}
end
Outline(gid::GID, o, s) where {GID} =
    Outline{GID}(gid, o, s, Matrix{Float64}(0,0), Int[])
Outline(gid::GID, o) where {GID} =
    Outline{GID}(gid, o, [1,size(o,2)+1] , Matrix{Float64}(0,0), Int[])
hasoutcrops(outl::Outline) = length(outl.outcrops)>0

"Approximate area of an outline"
function area(outline::Outline, gridpts=100)
    ol = outline.outline
    dx = (maximum(ol[1,:])-minimum(ol[1,:]))/gridpts
    dy = (maximum(ol[2,:])-minimum(ol[2,:]))/gridpts
    dx = min(dx,dy)
    x = minimum(ol[1,:]):dx:maximum(ol[1,:])
    y = minimum(ol[2,:]):dx:maximum(ol[2,:])
    pts = 0
    for xx in x, yy in y
        if inpoly((xx,yy), ol)
            pts +=1
        end
    end
    return pts*dx^2
    # TODO: remove outcrops
end

# ## refactor Glacier to use this instead of iscalving at some point
# """
# $(TYPEDEF)

# Holds glacier attributes, roughly following RGI.
# """
# @with_kw_noshow struct GlacierAttrs
#     form::GlacierAttrs.Form
#     terminus::GlacierAttrs.Terminus
#     surging::GlacierAttrs.surging = GlacierAttrs.not_surging
# end

# module GlacierAttrs
# @enum(Form,
#       glacier
#       glacier_stystem
#       icecap
#       icecap_catchment
#       )
# @enum(Terminus,
#       land,
#       marine,
#       lake,
#       dry_calving,
#       regenerated,
#       shelf)
# @enum(Surging,
#       not_surging,
#       possible,
#       probable,
#       observed)
# end



"""
Enumeration for the different types of teminus-flux specification available:

- abs: absolute value [m^3/s]
- rel_bdot_teminus: as ratio in terms of flux from integrating
                    expectation of prior of bdot
- rel_btilde_teminus: as ratio in terms of flux from integrating
                      expectation of prior of btilde
- rel_bdot_max: as ratio in terms of maximum flux from integrating
                expectation of prior of bdot
- none: do not look at the flux
"""
baremodule FluxKind
using Base
@enum _FluxKind abs rel_bdot_teminus rel_btilde_teminus rel_bdot_max none
end

"""
$(TYPEDEF)

Holds all input for a particular glacier, i.e. all our prior
knowledge of a glacier. This data stays fixed during fitting.

Note that no fields should contain NaNs
"""
@with_kw_noshow struct Glacier
    gid::GlacierID # RGI 4.0 number except for Peninsula
    iscalving::Bool # true for calving glacier (tide-water, fresh-water, hanging glacier)
    isicecap::Bool # true for ice-cap like glacier with several outlets
    proj4::String # map projection, a Proj4 string

    # DEM (assumed error free)
    dem::Gridded{F}         # surface DEM (cropped) [m].  Note that
                            # the dem.x,dem.y is also the
                            # computational 2D grid.
    outline::Outline # glacier outline and possibly additional outcrops
    cropbox::Tuple{Float64,Float64,Float64,Float64} # crop-box
    # "True on the glacier"
    glaciermask::Matrix{Bool}
    @assert size(dem)==size(glaciermask)
    # "True on solid ground"
    landmask::Matrix{Bool}
    @assert size(dem)==size(landmask)
    # "distance to closest land point (SurfaceTypes: land_margin)"
    dist_land::Matrix{F}
    # "dist_land was calculated up to this distance"
    dist_land_maxdist_calc::F
    # DEM products, all calculated
    alpha::Matrix{F}       # surface slope [radians] (prod)
    elemean::F
    elemax::F
    elemin::F
    elemedian::F

    # Distributed priors
    ###############
    # Measurements/estimates of input fields (1D, 2D or nothing)
    bdot::SomeData        # mass balance  [m ice /a]
    dhdt::SomeData        # surface elevation change [m ice/a]
    temp::SomeData        # ice temperature [C]
    fsl::SomeData         # basal sliding fraction []

    # Scalar priors
    ###############
    # Flux at terminus can be set to a few different settings:
    # - (:abs, SomeData(flux...)) : absolute flux at terminus [m^3/s].
    # - (:rel_bdot_teminus, SomeData(flux...)) : as fraction of the
    #    integrated bdot at the terminus --> calving
    # - (:rel_btilde_teminus, SomeData(flux...)) : as fraction of the
    #    integrated btilde at the terminus --> calving
    # - (:rel_bdot_max, SomeData(flux...)) : as fraction of the
    #    maximum of the integrated bdot --> calving or land terminating
    #
    # Probably (:abs, SomeData(0.0)) for land-terminating glaciers,
    # and (:rel_bdot_terminus, SomeData(1.0)), or similar, for
    # tide-water glaciers.  Sigma is small for land-terminating and
    # big for tide water.  Range: can restrict it to be non-negative
    # (thus h>0).
    #
    # Note that for all rel_* kinds the SomeData fields vals, sigma and range_abs,
    # are in fact relative to the "*".  E.g. sigma for rel_bdot_max, corresponds to
    # an absolute sigma of bdot_max*sigma.
    terminus_flux::Tuple{FluxKind._FluxKind, SomeData} # [m^3/s] or []

    ## Measurements of forward model output fields
    #####
    h::SomeData   # Traj{F}: measured ice thickness at radar lines, abuses .vals.err field to hold surface elevation [m]
    iv::SomeData    # iv velocity: measured surface flow speed: on different grid to DEM [m/a], possibly at points
    ivmask::Matrix{Bool} # true on the glacier, matches grid of .iv

    ###
    misc::Dict{Symbol,Any}=Dict{Symbol,Any}() # for stuff
end
## This does not work with @with_kw:
# function Glacier(args...)
#     mask = args[8]
#     for (i,arg) in enumerate(args)
#         if _hasnan(arg)
#             error("Argument #$i of Glacier constructor call contains NaNs! Not allowed to enable certain performance tricks.")
#         end
#         if _notmasked(arg, mask)
#             error("Argument #$i of Glacier constructor call contains unmasked FILL=$FILL values!")
#         end
#     end
#     Glacier(args...)
# end

getname(gl::Glacier) = getname(gl.gid)
getregion(gl::Glacier) = getregion(gl.gid)

Base.show(io::IO, gl::Glacier) = print(io, """Glacier $(getid(gl.gid)) with name "$(getname(gl))" in region "$(getregion(gl))" """)
area(gl::Glacier) = sum(gl.glaciermask) * step(gl.dem.x)^2

volume_area_mean_h(gl::Glacier) = volume_area_mean_h(area(gl), isicecap=gl.isicecap)

"Make sure there are no NaN lurking"
_hasnan(arg::AbstractArray{T}) where {T<:AbstractFloat} = any(isnan.(arg))
function _hasnan(arg::AGridded{T}) where {T<:AbstractFloat}
    VAWTools.haserror(arg) ? (_hasnan(arg.v) && _hasnan(arg.v)) : _hasnan(arg.v)
end
_hasnan(arg::SomeData) = hasdata(arg) ? _hasnan(arg.vals) : false
_hasnan(arg::AbstractFloat) = isnan(arg)
_hasnan(arg::Tuple) = mapreduce(_hasnan, |, false, arg)
_hasnan(arg::Any) = false

function _notmasked(arg, mask)
    if isa(arg,AbstractArray) && eltype(arg)<:AbstractFloat
        return any(arg[mask].==FILL)
    elseif isa(arg,AGridded) && eltype(arg.v)<:AbstractFloat
        return any(arg.v[mask]==FILL)
    elseif isa(arg,SomeData)
        hasdata(arg) ? _notmasked(arg.vals) : false
    elseif isa(arg, Tuple)
        return mapreduce(x->_notmasked(x,mask), |, false, arg)
    else
        return false
    end
end


# """
# Container which holds all needed parameters for a model run.
# """
# mutable struct Model
#     ph::Phys
#     pm:MPara
#     dirs::Dirs
#     gl::Glacier
#     num::Num
#     # res::Result
#     # Model(ph,pa,dirs,gl,num,res) = new(ph,pa,dirs,gl,num,res)
#     Model(ph,pa,dirs,gl,num)     = new(ph,pa,dirs,gl,num)
# end
