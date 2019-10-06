#####
# Synthetic Bench
#####
#
# For a Bench-like glacier as used in SHMIP
# https://shmip.bitbucket.io/instructions.html#sec-2-2

function LoadPara(gid::SyntheticGlacier; pl_kwargs...)
    pl_kwargs = deepcopy(pl_kwargs)
    pl = LoadPara(; gid=gid, pl_kwargs...)
    merge!((a,b)->a, pl.dataset_opts[ParaData][:pm],
           Dict(
               :bandsize => 30,
               :window_dem_smooth=>0,
               :window_width_smooth=>0, # can make a noisy tau
               :window=>100.0)
           )
    for DK in datakinds()
        pl.dataset_LOADERS[DK] = SyntheticBenchLoader
    end
    return pl
end

# surface functions
_bench_g(y) = 0.5e-6 * abs(y)^3
_bench_ginv(x) = (x/0.5e-6).^(1/3)
_bench_surface(x,y) = x<0 ? zero(x) : 100(x+200)^(1/4) + 1/60*x - 2e10^(1/4) + 1

# bed function from SHMIP:
const para_bench = 0.05  # gives depth of overdeepening
_bench_h(x, para) = (-4.5*x/6e3 + 5) * (_bench_surface(x,0)-_bench_f(x, para)) /
    (_bench_surface(x,0)-_bench_f(x, para_bench)+eps())
_bench_f(x,para) = (_bench_surface(6e3,0) - para*6e3)/6e3^2 * x^2 + para*x
_bench_bed(x,y, para) = _bench_f(x,para) + _bench_g(y) * _bench_h(x,para)

# outline, upper branch
_bench_outline(x) = _bench_ginv( (_bench_surface(x,0)-_bench_f(x,0.05))/(_bench_h(x,0.05)+eps()) )

# bed function from forward model

function load!(::DataSet{OutlineData, SyntheticBenchLoader}, gid::ANY)
    x = 0:10:6e3
    _outline = _bench_outline.(x)
    x = vcat(0, x[1:end-1], reverse(x), 0)
    y = vcat(0, -_outline[1:end-1], reverse(_outline), 0)
    outline = hcat(x, y)'
    iscalving = false
    isicecap = false
    proj = "n/a"
    outline = Outline(gid, outline)
    return outline, iscalving::Bool, isicecap::Bool, proj::String
end
const _dx = 25.0
const _buf = 10*_dx
function load!(ds::DataSet{DEMData, SyntheticBenchLoader}, outline,
               pp::Phys)::Tuple{Gridded,Matrix{Bool},Matrix{Bool},Matrix{F},F,Tuple,Any}
    opts = ds.opts
    grid = get(opts, :grid, (-_buf:_dx:6000+_buf, -500-_buf:_dx:500+_buf))
    cropbox = get(opts, :cropbox, (minimum(grid[1]), maximum(grid[1]), minimum(grid[2]), maximum(grid[1])))
    dist_land_maxdist = get(opts, :dist_land_maxdist, -10)

    dem = Gridded(grid[1], grid[2],_bench_surface.(grid[1], grid[2]'))
    glaciermask, landmask = maskit(dem, outline)
    dist_land = dist_to_mask(dem, landmask, glaciermask, dist_land_maxdist)
    proj = "n/a"
    return dem, glaciermask, landmask, dist_land, dist_land_maxdist, cropbox, proj
end

"""
opts:
- use_forward_model -- if true use fwd to generate IV field
  - if this is set then also pass in the model parameters:
    opts[:fwd_para] = (gl, pp, pm, pn)
  - as a hack it returns the function of the bed elevation inside
    the opts-struct: opts[:fwd_iv_fn] and opts[:fwd_iv_err_fn]
- iv_grid -- grid to return iv on
- err_sigma -- add a random error with this std
- err_cor_length -- and this correlation length
"""
function load!(ds::DataSet{IVData, SyntheticBenchLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::Tuple{SomeData,Matrix{Bool}}
    opts = ds.opts
    @get! use_forward_model = opts | false
    if use_forward_model
        gl, pp, pm, pn = opts[:fwd_para]
        @get! iv_grid = opts | (gl.dem.x[1:3:end], gl.dem.y[1:3:end])
        @get! err_cor_length = opts | 200.0
        @get! err_sigma = opts | 2.0


        # run forward model
        iv1d, iv2d = fwdm(gl, pp, pm, pn)[[3,6]]
        itp = Interp.interpolate((gl.dem.x, gl.dem.y), iv2d,
                                 Interp.Gridded(Interp.Linear()) )
        ivfn = (x,y) -> itp[x,y]
        opts[:fwd_iv_fn] = ivfn # return the "true" thickness function also
        opts[:fwd_true_iv1d] = iv1d

        # add some error (poor-man's spatially correlated errors)

        err_grid = (linspace(gl.dem.x[1], gl.dem.x[end], max(2, (gl.dem.x[end]-gl.dem.x[1])÷err_cor_length+1)),
                    linspace(gl.dem.y[1], gl.dem.y[end], max(2, (gl.dem.y[end]-gl.dem.y[1])÷err_cor_length+1)))


        itp_err = Interp.interpolate(err_grid, randn(length.(err_grid)...)*err_sigma,
                                 Interp.Gridded(Interp.Linear()) )
        err_fn = (x,y) -> itp_err[x,y]
        opts[:fwd_iv_err_fn] = err_fn

        # make it
        ivs = zeros(F, length.(iv_grid)...)
        ivmask = falses(length.(iv_grid)...)
        for (j,y) in enumerate(iv_grid[2])
            for (i,x) in enumerate(iv_grid[1])
                if 0<=x<=6e3 && _bench_outline(x)>=abs(y)
                    ivs[i,j] = ivfn(x,y) + err_fn(x,y)
                    # avoid negative values:
                    ivs[i,j] = max(ivs[i,j], 0)
                    ivmask[i,j] = true
                else
                    ivs[i,j] = 0
                    ivmask[i,j] = false
                end
            end
        end
        iv = SomeData(vals=Gridded(iv_grid..., ivs), sigma=err_sigma)
    else
        # only one point with velocity
        iv = SomeData(vals=Gridded(2000.0:2000, 0.0:0, 30.0+zeros(1,1)), sigma=5)
        ivmask = trues(1,1)
    end
    return iv, ivmask
end

"""
opts:
- use_forward_model=false -- if true, use the forward-model to generate
  the ice thickness
  - if this is set then also pass in the model parameters:
    opts[:fwd_para] = (gl, pp, pm, pn)
  - as a hack it returns the function of the bed elevation inside
    the opts-struct:  opts[:fwd_thick_fn] and opts[:fwd_thick_err_fn]
- para_bench=0.05 -- how much of an overdeepeing there is
- profile_xlocs=[2000.0,4000.0] -- locations of cross-profiles
- profile_ylocs=[50.0] -- locations of along-profiles
- profile_margin_dist=100.0 -- distance from margin at profile start
- profile_step=10.0 -- step between profile points
- err_sigma -- add a random error with this std
- err_cor_length -- and this correlation length
"""
function load!(ds::DataSet{ThicknessData, SyntheticBenchLoader}, dem,
               glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    opts = ds.opts
    @get! use_forward_model = opts | false
    if use_forward_model
        gl, pp, pm, pn = opts[:fwd_para]
        h1d, hs2d = fwdm(gl, pp, pm, pn)[[1,4]]
        itp = Interp.interpolate((gl.dem.x, gl.dem.y), hs2d,
                                 Interp.Gridded(Interp.Linear()) )
        hfn = (x,y) -> itp[x,y]
        opts[:fwd_thick_fn] = hfn # return the "true" thickness function also
        opts[:fwd_true_h1d] = h1d
    else
        @get! para_bench = opts | 0.05
        # TODO add length-profiles, add arbitrary trajectories
        hfn = (x,y) -> _bench_surface(x,y) - _bench_bed(x,y,para_bench)
    end

    # produce data along profiles:
    @get! profile_xlocs = opts | [850, 2000.0, 4000.0]
    @get! profile_ylocs = opts | [-100.0, 200.0]
    @get! profile_margin_dist = opts | 100.0
    @get! profile_step = opts | 10.0
    @get! err_sigma = opts | 10.0
    @get! err_cor_length = opts | 100.0
    @assert err_cor_length<1080 "err_cor_length too big"

    x = F[]
    y = F[]
    thick = F[]
    surf = F[]
    splits = UnitRange{Int}[]
    hfns_err = []
    for xl in profile_xlocs
        w = _bench_outline(xl)-profile_margin_dist
        ys = -w:profile_step:w
        xs = ys*0 + xl

        # make an error
        ys_err = linspace( ys[1], ys[end], max(2, (ys[end]-ys[1])÷err_cor_length+1))
        hfn_err = let
            tmp = Interp.interpolate((ys_err,), randn(length(ys_err))*err_sigma,
                                     Interp.Gridded(Interp.Linear()) )
            (x,y) -> tmp[y]
        end

        push!(hfns_err, hfn_err)


        s = _bench_surface.(xs,ys)
        push!(splits, length(thick)+1:length(thick)+length(xs))
        append!(surf, s)
        append!(thick, hfn.(xs,ys) + hfn_err.(xs,ys))
        append!(x,xs)
        append!(y,ys)
    end

    for yl in profile_ylocs
        # find starting point
        if abs(yl)<75
            # the outline has \inf slope at 0,6e3:
            x0, x1 = 1.0, 6e3-5
        else
            rootfn = x -> abs(yl) - _bench_outline(x)
            x0 = Roots.fzero(rootfn, 1e-3, 2e3)
            x1 = Roots.fzero(rootfn, 4e3, 6e3-1e-3)
        end
        xs = x0+profile_margin_dist:profile_step:x1-profile_margin_dist
        ys = 0*xs + yl

        # make an error
        xs_err = linspace( xs[1], xs[end], max(2, (xs[end]-xs[1])÷err_cor_length+1))
        hfn_err = let
            tmp = Interp.interpolate((xs_err,), randn(length(xs_err))*err_sigma,
                                     Interp.Gridded(Interp.Linear()) )
            (x,y) -> tmp[y]
        end
        push!(hfns_err, hfn_err)

        s = _bench_surface.(xs,ys)
        push!(splits, length(thick)+1:length(thick)+length(xs))
        append!(surf, s)
        append!(thick, hfn.(xs,ys) + hfn_err.(xs,ys))
        append!(x,xs)
        append!(y,ys)
    end
    opts[:fwd_thick_err_fn] = hfns_err

    thickness = SomeData(vals=Traj(x,y, thick, surf, splits),
                         sigma=err_sigma)
    return thickness::SomeData
end

"""
Construct a plausible SMB field:

- ELA set to median elevation (2D) (can be set)
- gradient 0.008 (range 0.002 (continental) to 0.02 (maritime))
- gradient reduction in accumulation area: 0.55 (range 0.1 to 1)
"""
function load!(ds::DataSet{BdotData, SyntheticBenchLoader}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    opts = ds.opts
    ela_elevation = get(opts, :ela_elevation, median(dem.v[glaciermask]))
    gradient = get(opts, :smb_gradient, 0.008)
    gradient_corr_accumulation = get(opts, :gradient_corr_accumulation, 0.55)
    bd = dem.v-ela_elevation
    for i in eachindex(bd)
        b = bd[i]
        bd[i] = b>0 ? b*gradient*gradient_corr_accumulation : b*gradient
    end
    # scale to m ice / a
    scale_m_ice!(bd, pp)

    bdot = SomeData(vals=Gridded(dem.x, dem.y, bd),
                    sigma=get(opts, :sigma, 2.0))
    return bdot::SomeData
end

"""
Makes a 2D dh/dt field such that: first dh/dt= 0 from top to ELA
elevation, then such that flux at terminus is given.

NOTE: this can give a glacier of zero thickness in the middle as the
flux can go negative, as long as it comes back to the given terminus
flux.

- ELA is taken where bdot==0.  Or can be specified in ds.opts[:ela_elevation].
  Note that if the ELA is inconsistent with bdot, then this produces a inconsistent
  dhdt in the sense that flux will go negative for some of the glacier.
- scaled such that terminus flux is ok (or set manually with
  ds.opts[:dhdt_scale].

- constraints: dh/dt in accum<=bdot, dh/dt in abl >=bdot
  - at least for the 1D case
- then just the usual uncertainty on btilde?
"""
function load!(ds::DataSet{DhdtData, SyntheticBenchLoader},
               bdot, terminus_flux,
               dem, glaciermask, target_proj, outline, cropbox,
               iscalving, pp::Phys)
    opts = ds.opts
    dx = step(dem.x)

    twoD_bdot = if isa(bdot.vals, Gridded) # 2D bdot
        bdot_2D = bdot.vals.v
        true
    elseif isa(bdot.vals, Gridded1d) # 1D bdot
        error("Not supported. This introduces mistakes, in particular for small glaciers.  Extrapolate a 1D bdot to 2D when loading.")
        false
    else
        error("Not supported type of bdot: $bdot, typeof = $(typeof(bdot))")
    end

    # get ELA
    if haskey(opts, :ela_elevation)
        # set by hand
        ela_elevation = opts[:ela_elevation]
    else
        if twoD_bdot
            bands, bandi = VAWTools.bin_grid(dem, 30, glaciermask, min_bands=4)
            bdot1d = VAWTools.map_onto_bands(bandi, bdot.vals.v, mean, FILL)
        else
            bdot1d = bdot.vals.v
            bands = bdot.vals.x
        end
        i = findlast(bdot1d.<=0)
        if i==0
            ela_elevation = minimum(dem.v[glaciermask])
        elseif i==length(bdot1d)
            ela_elevation = maximum(dem.v[glaciermask])
        else
            ela_elevation = bands[i]
        end
    end
    has_accumulation_area = ela_elevation < maximum(dem.v[glaciermask])
    has_ablation_area = ela_elevation > minimum(dem.v[glaciermask])

    msk_bdot = glaciermask .& (bdot_2D.!=FILL) # in case there are masked bdot values within glacier
    smb = sum(bdot_2D[msk_bdot]*dx^2)

    if has_accumulation_area
        if smb<0
            # Trial dhdt has linear dependence on elevation in ablation zone, constant in accumulation zone
            dd = dem.v-ela_elevation
            # above ELA dh/dt== small.  Do not use 0 when there is no ablation area.
            small = has_ablation_area ? 0.0 : 0.1
            for i in eachindex(dd)
                d = dd[i]
                dd[i] = d>0 ? small : d
            end
            dd[.!glaciermask] .= FILL
        else # growing glacier
            dd = dem.v-ela_elevation
        end
    else
        # Trial dhdt is flat
        dd = dem.v*0 + 1
        dd[.!glaciermask] .= FILL
    end
    # Now scale such that flux is satisfied
    if haskey(opts, :dhdt_scale)
        scale = opts[:dhdt_scale]
    else
        if twoD_bdot
            bdot_2D = bdot.vals.v
        else
            error("needs updating")
            # extrapolate to 2D
            bands, bandi = VAWTools.bin_grid(dem, -step(bdot.vals.x), glaciermask) #, min_bands=0)
            bdot_2D = VAWTools.map_back_to_2D(size(dem.v), bandi, reverse(bdot.vals.v))
            ## check that collapsing back to 1D is identity
            # @show tmp_ = VAWTools.map_onto_bands(bandi, bdot_2D)
            # @show bdot.vals.v[isnan.(tmp_)]
            # @assert all(tmp_[length.(bandi).>0] .≈ bdot.vals.v[length.(bandi).>0])
        end

        # Scale such that flux at terminus as requested
        fkind = terminus_flux[1]
        vals = terminus_flux[2].vals
        if fkind==FluxKind.abs
            # absolute flux
            target_flux = vals
        elseif fkind==FluxKind.rel_bdot_teminus # target_flux at terminus relative to flux due to expectation bdot
            bdotflux_at_terminus = sum(bdot_2D[glaciermask])*dx^2
            target_flux = bdotflux_at_terminus*vals
        elseif fkind==FluxKind.rel_bdot_max # maximum flux due to expectation bdot (== total accumulation rate)
            inds = (bdot_2D.>0) .& glaciermask
            max_bdotflux = sum(bdot_2D[inds])*dx^2
            target_flux = max_bdotflux*vals
        elseif fkind==FluxKind.rel_btilde_teminus
            error("Flux-kind: $fkind not supported with this loader $ds")
        elseif fkind==FluxKind.none && !iscalving
            # land terminating
            target_flux = 0.0
        else
            error("Unknown flux-kind: $fkind")
        end

        unscaled_dhdt = sum(dd[msk_bdot]*dx^2)
        scale = (smb-target_flux)/unscaled_dhdt
    end
    # TODO: here a offset instead of a scale could also be used.
    #       But how to determine its value?  Maybe minimize sum(abs(dhdt))?
    dd *= scale
    dd[.!glaciermask] .= FILL
    dhdt = SomeData(vals=Gridded(dem.x, dem.y, dd),
                    sigma=get(opts, :sigma, 2.0))
    return dhdt::SomeData
end
function load!(ds::DataSet{FslData, SyntheticBenchLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    opts = ds.opts
    fsl = get(opts, :fsl, [0.5,0.5] ) # sliding at top and bottom
    top = maximum(dem.v[glaciermask])
    bottom = max(0,0, minimum(dem.v[glaciermask]))
    ele = linspace(bottom, top, length(fsl))
    fsl = SomeData(vals=Gridded1d(ele, collect(fsl)),
                   range_abs=(0,0.999),
                   sigma=Inf)
    return fsl::SomeData
end
function load!(ds::DataSet{TempData, SyntheticBenchLoader}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    opts = ds.opts
    temp = get!(opts, :temp, [-0.5, -0.5])
    top = maximum(dem.v[glaciermask])
    bottom = max(0,0, minimum(dem.v[glaciermask]))
    ele = get!(opts, :ele, linspace(bottom, top, length(temp)))
    temp = SomeData(vals=Gridded1d(ele, collect(temp)),
                    sigma=1,
                    range_abs=(-20,0))
    return temp::SomeData
end
function load!(ds::DataSet{TerminusFluxData, SyntheticBenchLoader}, dem, glaciermask, iscalving, pp::Phys)
    fluxkind = get(ds.opts, :fluxkind, FluxKind.rel_bdot_max)
    if iscalving
        # assume that half the balance is lost through flux
        terminus_flux = get(ds.opts, :terminus_flux, SomeData(vals=1/2,
                                                              sigma=0.1,
                                                              range_abs=(0.005,10.0)
                                                              )
                            )
    else
        # # approximately zero flux at terminus +/- 1% of max bdot flux accumulation
        # terminus_flux = get(ds.opts, :terminus_flux, SomeData(vals=0.0,
        #                                                       sigma=0.001
        #                                                       )
        #                     )
        fluxkind = FluxKind.none
        terminus_flux = SomeData(nothing)
    end

    return fluxkind, terminus_flux
end
