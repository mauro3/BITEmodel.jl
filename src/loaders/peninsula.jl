# Loading of Matthias' Peninsula data
#
# big ones: 49, 197, 281, 402, 422, 428

function LoadPara(gid::PeninsulaGlacier; pl_kwargs...)
    glacier = getid(gid)
    pl_kwargs = deepcopy(pl_kwargs)
    pl = LoadPara(; gid=gid, data_root=joinpath(LoadPara().data_root,
                                          "Peninsula"), pl_kwargs...)
    for (k,v) in (ThicknessData=>10.0,
                  BdotData=>0.2,
                  DhdtData=>0.2,
                  IVData=>1.0)
        pl.dataset_opts[k][:sigma]=v
    end
    updated_pm_defaults = Dict(
        :bandsize => 30,
        :window_width_smooth=>100, # can make a noisy tau
        :cap_at_floatation=>false,
        :dist_exp => 0.5,
        :h_at_top => -1,
    )

    # update pm but keeping any explicitly set (via kwargs) options
    merge!((a,b)->a, pl.dataset_opts[ParaData][:pm], updated_pm_defaults)

    # Set loaders & files
    id = getid(gid)
    id5 = VAWTools.int2str5(id)
    pl = deepcopy(pl)

    dl = pl.dataset_LOADERS
    fl = pl.dataset_files
    # Set the LOADERS
    dl[OutlineData] = PeninsulaLoader
    if isfile(joinpath(pl.data_root, "outline/outcrops/$id5.xyn"))
        fl[OutlineData] = [joinpath(pl.data_root, "outline/$id5.xyn"),
                           joinpath(pl.data_root, "outline/outcrops/$id5.xyn")]
    else
        fl[OutlineData] = [joinpath(pl.data_root, "outline/$id5.xyn")]
    end
    dl[ParaData] = PeninsulaLoader

    dl[DEMData] = PeninsulaLoader
    fl[DEMData] = [joinpath(pl.data_root, "dem/APDEM100m.bin")]

    dl[IVData] = PeninsulaLoader
    fl[IVData] = [joinpath(pl.data_root, "velocity/speed_measured.grid")]

    dl[BdotData] = PeninsulaLoader
    fl[BdotData] = [joinpath(pl.data_root, "RACMO/downscale_RACMO_nn.bin")]

    dl[DhdtData] = PeninsulaLoader

    dl[ThicknessData] = PeninsulaLoader
    f = joinpath(pl.data_root, "icebridge/single_glacier/$id5.dat")
    if isfile(f)
        fl[ThicknessData] = [f]
    end

    dl[FslData] = SyntheticBenchLoader
    dl[TempData] = SyntheticBenchLoader
    dl[TerminusFluxData] = SyntheticBenchLoader


    return pl
end

function load!(ds::DataSet{OutlineData, PeninsulaLoader}, gid)
    iscalving = true
    isicecap = false
    proj = ""

    # outlines & outcrops
    outline, splits = VAWTools.concat_poly(read_xyn(ds.files[1]))

    if length(ds.files)==2
        outcrops, splits_oc = VAWTools.concat_poly(read_xyn(ds.files[2], fix=true))

        return Outline(gid, outline, splits, outcrops, splits_oc),
          iscalving::Bool, isicecap::Bool, proj::String
    else
        return Outline(gid, outline, splits),
          iscalving::Bool, isicecap::Bool, proj::String
    end
end

function load!(ds::DataSet{DEMData, PeninsulaLoader}, outline::Outline, pp::Phys
               )::Tuple{Gridded, Matrix{Bool}, Matrix{Bool}, Matrix{F}, F, Tuple, String}
    load!(convert(DataSet{DEMData, :ASCIIGridLoader}, ds), outline, pp)
end
function generate_DEM_masks(gid::PeninsulaGlacier, dem, outline)
    if hasoutcrops(outline)
        glaciermask, landmask = maskit(dem, outline, outline.outcrops)
    else
        glaciermask, _ = maskit(dem, outline)
        landmask = convert(Matrix{Bool}, falses(size(glaciermask)))
    end
    return glaciermask, landmask
end


function load!(ds::DataSet{IVData, PeninsulaLoader}, dem, glaciermask, target_proj,
               outline, cropbox, pp::Phys)::Tuple{SomeData, Matrix{Bool}}
    load!(convert(DataSet{IVData, :ASCIIGridLoader}, ds),
          dem, glaciermask, target_proj, outline, cropbox, pp)
end

function load!(ds::DataSet{ThicknessData, PeninsulaLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::SomeData
    @unpack files, opts, daterange = ds
    sigma = opts[:sigma]
    if isfile(files[1])
        out = get_cache_ram!(ds, :thick, ()->readdlm(files[1]))
    else
        return SomeData(nothing)
    end

    if size(out,2)==5
        out = out[out[:,5].==0,:] # drop wall reflection points
    end
    tr = Traj(out[:,1], out[:,2], out[:,4], out[:,3])
    split_traj!(tr, 1e3)
    return SomeData(vals=tr, sigma=sigma)
end

function load!(ds::DataSet{BdotData, PeninsulaLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    bdot = load_asciigrid!(ds, :bdot, cropbox)

    # fill in any gaps inside the glacier
    fillgaps!(bdot, glaciermask)

    # scale from mm w.e./a to m ice / a
    scale_m_ice!(bdot, pp, 1/1000)

    return SomeData(vals=bdot, sigma=ds.opts[:sigma])
end

"This is a synthtetic dh/dt as in H&F 2014"
function load!(ds::DataSet{DhdtData, PeninsulaLoader},
               bdot, terminus_flux,
               dem, glaciermask, target_proj, outline, cropbox,
               iscalving, pp::Phys)
    opts = ds.opts
    cutoff = get(opts, :cutoff, 700)
    base = get(opts, :base, -1)
    sigma = opts[:sigma]

    # TODO: larger dh/dt for glaciers feeding into Larson A & B
    out = (base .- dem.v .* base./cutoff) # m ice/a
    out[dem.v.>cutoff] = 0
    return SomeData(vals=Gridded(dem.x, dem.y, out),
                    sigma=1.0)
end

function load!(ds::DataSet{FslData, PeninsulaLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    error("FSL data not provided in Peninsula.  Use different loader, probably syntetic")
end
function load!(ds::DataSet{TempData, PeninsulaLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    error("Ice-temp data not provided in Peninsula.  Use different loader, probably syntetic")
end
function load!(ds::DataSet{TerminusFluxData, PeninsulaLoader}, dem, glaciermask, iscalving,
               pp::Phys)
    error("Calving flux not provided in Peninsula.  Use different loader, probably synthetic")
end

##############
# helpers
