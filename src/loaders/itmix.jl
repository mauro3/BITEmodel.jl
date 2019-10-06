##########
## ITMIX (of all versions loaders)
#########
#

#####
# ITMIX data
#####
"Translate part of the file-name to DataKind subtypes"
const dk_dict = Dict("RES"=>ThicknessData,
                     "margin"=>OutlineData,
                     "surface"=>DEMData,
                     "mb"=>BdotData,
                     "dhdt"=>DhdtData,
                     "speed"=>IVData,
                     "azimuth"=>nothing,
                     "velocity"=>IVData)

todo = 1
"Whether a glacier is calving, is an ice cap, and its RGI 4.0 number (if there is a single one)."
const iscalving_isicecap = Dict(:Academy=>(true,true,nothing),
                                :Austfonna=>(true,true,nothing),
                                :Devon=>(true,true,"RGI40-03.02443"), # one of the major outlet glaciers in the North
                                :Hellstugubreen=>(false,false,"RGI40-08.02182"),
                                :Mocho=>(false,false,nothing), # glacier complex
                                :SouthGlacier=>(false,false,"RGI40-01.17598"),
                                :Tasman=>(false,false,"RGI40-18.02342"),
                                :Washmawapta=>(false,false,"RGI40-02.04894"),
                                :Aqqutikitsoq=>(false,false,"RGI40-05.00524"), # not sure, size does not match well
                                :Brewster=>(false,false,"RGI40-18.01130"),
                                :Elbrus=>(false,false,nothing), # glacier complex
                                :NorthGlacier=>(false,false,"RGI40-01.18240"),
                                :Starbuck=>(true,false,nothing),
                                :Unteraar=>(false,false,"RGI40-11.01328"),
                                :Columbia=>(true,false,"RGI40-01.11506"),
                                :Freya=>(false,false,"RGI40-05.19983"),
                                :Kesselwandferner=>(false,false,"RGI40-11.00787"),
                                :Urumqi=>(false,false,"RGI40-13.27054"),
                                :Synthetic1=>(false,false,nothing),
                                :Synthetic2=>(false,false,nothing),
                                :Synthetic3=>(false,true,nothing),
                                :AustreGroenfjordbreen=>(false,false,nothing),
                                :ChhotaShigri=>(false,false,nothing),
                                )

# ITMIX glacier with IV data:
const itmix_glaciers_iv = [:Synthetic1,  :Unteraar,
                           :Synthetic2, :Synthetic3,
                           :Tasman, :Austfonna, :Devon,
                           :SouthGlacier, :NorthGlacier,
                           :Brewster, :Hellstugubreen]

# of which these only have point data (at GPS locations)
const itmix_glaciers_iv_point = [:NorthGlacier, :SouthGlacier, :Hellstugubreen, :Brewster]

function LoadPara(gid::ITMIXGlacier; use_glogem=false, pl_kwargs...)
    pl_kwargs = deepcopy(pl_kwargs)
    id = getid(gid)
    data_root = gid.version==1? joinpath(LoadPara(gid=gid).data_root,"ITMIX") : joinpath(LoadPara(gid=gid).data_root,"ITMIX2")
    pl = LoadPara(; gid=gid, data_root=data_root, pl_kwargs...)
    # set the sigmas of uncertainty for the data:
    for (k,v) in Dict(ThicknessData=>10.0, # maybe change to more
                      BdotData=>0.2,  # m/a, 1m/a is what Matthias says about GloGEM
                      DhdtData=>0.2,  # m/a, e.g. 1m/a Nuth et al. 2010, Table 1, for Svalbard
                      IVData=>10.0)
        pl.dataset_opts[k][:sigma]=v
    end
    updated_pm_defaults = Dict(
        :bandsize => 30,
        :window_width_smooth=>100, # can make a noisy tau
        :cap_at_floatation=>false,
        :dist_exp => 0.5,
        #:window => 0
    )
    if iscalving_isicecap[id][2]
        merge!(updated_pm_defaults, Dict(:h_at_top => -1, :shapeF => false))
    end
    if iscalving_isicecap[id][1]
        # relative flux lost through calving
        pl.dataset_opts[TerminusFluxData][:fluxkind] =
            FluxKind.rel_btilde_teminus # compare to prior expectation of btilde integral
        #pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=1.0, sigma=0.5, range_abs=(0.1,3.0))
        pl.dataset_opts[TerminusFluxData][:terminus_flux] =
            SomeData(vals=NaN) # require this to be set individually
    else
        # approximately zero flux at terminus +/- 1% of max bdot flux accumulation
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_max
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=0.0, sigma=0.01)
    end

    if id==:Austfonna
        # This makes a coarser grid as the Austfonna data is on a
        # 50m grid!
        # with 50m DEM:
        #  - 400s without cache
        #  - 35s with cache
        # with 250m DEM:
        #  - 208s without cache
        #  - 35s  with cache
        if gid.version==1
            pl.dataset_opts[DataKind][:grid_downsample_step]=5
        end
        #pl.dataset_LOADERS[DhdtData] = SyntheticBenchLoader # something is off with Austfonna's dh/dt (or bdot)

        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_max # compare to prior expectation of bdot integral
        # Dowdeswell at al. 2008 says this:
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=1/3,
                                                                     sigma=0.4,
                                                                     range_abs=(0.05, 3.0))
    elseif id==:AustreGroenfjordbreen
        # compare to prior expectation of btilde integral because rel_bdot_max is negative
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_btilde_teminus
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=0.0, sigma=0.01)
    elseif id==:Columbia
        # the size of (1040, 920) make it run slow, about 0.5s per forward run
        # reduced by step=2 (4x) it's down to 0.07, in-line with the others.
        pl.dataset_opts[DataKind][:grid_downsample_step] = 2

        pl.dataset_opts[BdotData][:smb_gradient] = 0.02 # maritime
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_max # compare to prior expectation of max-bdot integral
                                                                             # Note, bdot intergral is negative at terminus
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=0.5,
                                                                     sigma=0.3,
                                                                     range_abs=(0.1,3.0))
        merge!(updated_pm_defaults, Dict(:h_at_top => -1)) # more ice cap like at devide
    elseif id==:Aqqutikitsoq
        if gid.version==1
            pl.dataset_opts[DataKind][:grid_downsample_step]=2
        end
    elseif id==:Hellstugubreen
        # pl.dataset_LOADERS[DhdtData] = SyntheticBenchLoader
        # pl.dataset_opts[DhdtData][:ela_elevation] = 2000
    elseif id==:Elbrus
        if gid.version==1
            pl.dataset_opts[DataKind][:grid_downsample_step]=3
        else
            # the size of (601,581) make it run slow, about 0.12s per forward run
            # reduced by step=2 (4x) it's down to 0.025, in-line with the others.
            # (now on 60m resolution vs 30m resolution)
            pl.dataset_opts[DataKind][:grid_downsample_step] = 2
        end


        pl.dataset_LOADERS[BdotData] = SyntheticBenchLoader # original only covers a small fraction of the area
                                                            # TODO: this could be turned into elevation band data:
                                                             # pl.dataset_opts[BdotData][:vals_prior] = 0.0
        # pl.dataset_LOADERS[DhdtData] = SyntheticBenchLoader # doesn't make much sense to then not use this either.
        merge!(updated_pm_defaults, Dict(:h_at_top => -1, :shapeF => false))
    elseif id==:Unteraar
        # something is amiss with dhdt or bdot data as glacier stops at confluence:
        # pl.dataset_LOADERS[DhdtData] = SyntheticBenchLoader
        # pl.dataset_opts[DhdtData][:ela_elevation] = 2900
        pl.dataset_opts[IVData][:sigma] = 4.0 # Vogel et al 2012
        pl.dataset_opts[ThicknessData][:sigma] = 27.0 # Bauder etal 2003
    elseif id==:Freya
        # uses the SyntheticBenchLoader
        pl.dataset_opts[DhdtData][:ela_elevation] = 985 # gets negative fluxes in the middle otherwise
    elseif id==:Urumqi
        pl.dataset_opts[DhdtData][:ela_elevation] = 4200
    elseif id==:Tasman
        merge!(updated_pm_defaults, Dict(:bandsize => 35))
    elseif id==:Devon
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_max # compare to prior expectation of bdot integral
        # https://www.igsoc.org/annals/46/a46a288.pdf says 1/4 calving and 3/4 melt
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=1/4,
                                                                     sigma=0.1,
                                                                     range_abs=(0.1,3.0))
    elseif id==:Synthetic3
        merge!(updated_pm_defaults, Dict(:h_at_top => -1))
    elseif id==:Synthetic1 || id==:Synthetic2 || id==:Synthetic3
        # there are no "measurement" errors here as no noise was put on the data.
        # Note there is still the BITEModel-error in pm.sigma_*
        pl.dataset_opts[IVData][:sigma] = 0.0
        pl.dataset_opts[ThicknessData][:sigma] = 0.0
    elseif id==:ChhotaShigri
        # flux not zero at terminus
    elseif id==:Academy
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_max
        # Dowdeswell at all 2008 says this:
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=1/3,
                                                                     sigma=0.3,
                                                                     range_abs=(0.1,3.0))
    elseif id==:Mocho
        # has negative mass-balance at the top!
        merge!(updated_pm_defaults, Dict(:h_at_top => -1, :shapeF => false))
    elseif id==:Starbuck
        pl.dataset_opts[TerminusFluxData][:fluxkind] = FluxKind.rel_bdot_teminus # compare to prior expectation of bdot integral
        pl.dataset_opts[TerminusFluxData][:terminus_flux] = SomeData(vals=1.0,
                                                                     sigma=0.1,
                                                                     range_abs=(0.1,3.0))
        merge!(updated_pm_defaults, Dict(:h_at_top => -1))  # more ice cap like at devide
    end
    pl.dataset_opts[DEMData][:cropbox] = :from_DEM

    # update pm but keeping any explicitly set (via kwargs) options
    merge!((a,b)->a, pl.dataset_opts[ParaData][:pm], updated_pm_defaults)

    # GloGEM
    if use_glogem && iscalving_isicecap[gid.id][3]!=nothing
        error("glogem loader disabled for ITMIX2 runs")
        pl.dataset_LOADERS[BdotData] = GloGEMLoader
        pl.dataset_LOADERS[DhdtData] = SyntheticBenchLoader
    else
        if use_glogem
            error("Not using GloGEM for $gid")
        end
    end

    # Set loaders and files:
    @unpack dataset_LOADERS, dataset_files, dataset_opts = pl
    # check what files and thus data exist
    name = getname(gid)
    dir = joinpath(pl.data_root, name)
    zip = joinpath(pl.data_root, name*".zip")
    fls = Dict()
    if isfile(zip) # zip-archive
        zipfl = ZipFile.Reader(zip)
        toread = zipfl.files
    elseif isdir(dir)
        toread = readdir(dir)
        append!(toread, joinpath.("shapefiles", readdir(joinpath(dir, "shapefiles"))))
    else
        error("Cannot read $dir or $zip")
    end
    for fl in toread
        out = parse_ITMIX_filename(fl)
        if out!=nothing
            fls[out[1]] = out[2:end]
        end
    end
    isfile(zip) && close(zipfl)
    proj_ = fls[OutlineData][end]
    args = []

    # Set the LOADERS
    for DK in datakinds()
        if !haskey(dataset_LOADERS, DK)
            if haskey(fls,DK) # ITMIX file present
                fln, nr, nameinfile, date_range, proj = fls[DK]
                dataset_LOADERS[DK] = ITMIXLoader
                dataset_files[DK] = [joinpath(dir,fln)]
                dataset_opts[DK][:proj] = proj_
            else
                if DK==ParaData
                    dataset_LOADERS[DK] = ITMIXLoader
                    if iscalving_isicecap[Symbol(name)][2]
                        dataset_opts[DK][:iv_extrapolation_scaling] = :huss
                    end
                # elseif DK==BdotData # TODO: use GloGEM by default
                #     dataset_LOADERS[DK] = GloGEMLoader
                #     dataset_files[DK] = [joinpath(dir,fln)]
                elseif DK==IVData || DK==ThicknessData
                    dataset_LOADERS[DK] = NothingLoader
                else
                    # if there is no file use the synthetic loader
                    dataset_LOADERS[DK] = SyntheticBenchLoader
                end
            end
        else
            # custom set loader & other bits, leave as is
        end
    end

    # bug in shapefile
    if id==:Tasman
        dataset_files[OutlineData] = ["data/ITMIX/Tasman/01_margin_Tasman_1986_UTM-59.csv"]
    end
    return pl
end

# Outline
function load!(ds::DataSet{OutlineData, ITMIXLoader}, gid)
    opts = ds.opts
    proj = opts[:proj]
    # fl, fid = open_file(ds.files[1])
    # outline = readdlm(fl, ';', Float64)'
    # close(fid)
    @assert length(ds.files)==1
    if endswith(ds.files[1], ".shp")
        outline, splits = load_shapefile!(ds, :outline)
    else
        outline = @open_magic readdlm(fl, ';', F)'

        if size(outline,1)==3
            # TODO
            # # process boundary ID stuff:
            # # TODO: not sure this is right yet.
            # boundary_ind = outline[3,:]
            # xy = outline[1:2,:]
            # # all boundary_ind>1 treat as outcrops
            # outline = xy[:,boundary_ind.==boundary_ind[1]]
            # outcropsx = Float64[]
            # outcropsy = Float64[]
            # splits = Int[]
            # for i=2:1000
            #     inds = boundary_ind.==i
            #     length(inds)==0 && break
            #     append!(outcropsx, xy[1,inds])
            #     append!(outcropsy, xy[2,inds])
            # end
            # outline = xy
            outline = outline[1:2,:]
        end
        if outline[:,1] !=outline[:,end]
            # close polygon
            outline = hcat(outline, outline[:,1])
        end
        splits = VAWTools.find_poly_splits(outline)
    end
    # Fix Unteraar
    if getname(gid) == "Unteraar"
        outline = outline[:, splits[end-1]:end-1]
        # remove self-intersections
        outline = outline[:,[1:1483, 1487:1667, 1671:end;]]
        splits = [1,size(outline)[2]]
        @assert outline[:,1]==outline[:,end]
    end

    iscalving,isicecap,rgi = iscalving_isicecap[getid(gid)]
    @assert isa(iscalving, Bool) "Input whether calving or not in iscalving_dict"
    return Outline(gid, outline, splits),
           iscalving::Bool, isicecap::Bool, proj::String
end

"""
opts:
- grid : if specified use a different grid from the grid of the DEM as
  computational grid.
- cropbox -- used for cropping the DEM and other 2D data.  Needs to
  contain all of the outline. `F[xmin, xmax, ymin, ymax]`, defaults
  being calculated during DEM-load

Return:
- dem::Gridded
- glaciermask::Matrix{Bool} -- true on glacier
- landmask::Matrix{Bool} -- true on solid ground (i.e. false on glaciers, ice shelves, water)
- dist_land::Matrix{F} -- distance to closest land
- dist_land_maxdist -- maximum distance to which to calculate above
- cropbox::Vector{Float64} -- box around region of interest (contains glacier + some buffer)
                              `F[xmin, xmax, ymin, ymax]`
"""
function load!(ds::DataSet{DEMData, ITMIXLoader}, outline::Outline, pp::Phys
               )::Tuple{Gridded, Matrix{Bool}, Matrix{Bool}, Matrix{F}, F, Tuple, String}
    @assert ds.opts[:cropbox]==:from_DEM "ITMIX cropbox needs to be: `:from_DEM`"
    return load!(convert(DataSet{DEMData, :ASCIIGridLoader}, ds), outline, pp)
end
function load!(ds::DataSet{IVData, ITMIXLoader}, dem, glaciermask, target_proj,
               outline, cropbox, pp::Phys)::Tuple{SomeData, Matrix{Bool}}
    if splitext(ds.files[1])[2]==".txt"
        # IV at points
        tmp = @open_magic readdlm(fl, header=true)
        speed = F.(tmp[1][:,4])
        x = F.(tmp[1][:,end-2])
        y = F.(tmp[1][:,end-1])
        iv = Traj(x,y,speed)
        ivmask = Bool.(trues(0,0))
        return SomeData(vals=iv, sigma=ds.opts[:sigma]), ivmask
    elseif splitext(ds.files[1])[2]==".asc"
        # IV in 2D
        return load!(convert(DataSet{IVData, :ASCIIGridLoader}, ds), dem, glaciermask, target_proj,
                     outline, cropbox, pp)
    else
        error("Extension of file $file not supported.  Need .txt or .asc")
    end
end
function load!(ds::DataSet{ThicknessData, ITMIXLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::SomeData
    @unpack files, opts, daterange = ds
    @unpack gid = outline
    sigma = opts[:sigma]
    if gid.version==1
        zb = get_cache_ram!(ds, :zb, @open_magic_fn readdlm(fl; header=true)[1])
        x,y = zb[:,1], zb[:,2]
        zb = zb[:,3]
        izs  = Interp.interpolate((dem.x, dem.y), dem.v, Interp.Gridded(Interp.Linear()) )
        zs = [izs[xx,yy] for (xx,yy) in zip(x,y)]
        thick = zs-zb
        out = Traj(x, y, thick, zs)
        split_traj!(out, 1e3) # TODO adjust split-distance
    else
        zb = get_cache_ram!(ds, :zb, @open_magic_fn readdlm(fl; header=true, skipstart=16)[1])
        lines = unique(zb[:,1])
        x,y,thick,zs = zb[:,3], zb[:,4], zb[:,end], zb[:,5]
        splits = []
        for l in lines
            push!(splits, findfirst(zb[:,1], l):findlast(zb[:,1], l))
        end
        out = Traj(x, y, thick, zs, splits)
    end
    return SomeData(vals=out, sigma=sigma)
end
function load!(ds::DataSet{BdotData, ITMIXLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    bdot = load_asciigrid!(ds, :bdot, cropbox)
    # scale to m ice / a
    scale_m_ice!(bdot, pp)
    if haskey(ds.opts, :vals_prior)
        return SomeData(vals=bdot, sigma=ds.opts[:sigma], vals_prior=ds.opts[:vals_prior])
    else
        return SomeData(vals=bdot, sigma=ds.opts[:sigma])
    end
end
function load!(ds::DataSet{DhdtData, ITMIXLoader},
               bdot, terminus_flux,
               dem, glaciermask, target_proj, outline, cropbox,
               iscalving, pp::Phys)
    return load!(convert(DataSet{DhdtData, :ASCIIGridLoader}, ds),
                 dem, glaciermask, target_proj, outline, cropbox, pp)
end

function load!(ds::DataSet{FslData, ITMIXLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    error("FSL data not provided in ITMIX.  Use different loader, probably syntetic")
end
function load!(ds::DataSet{TempData, ITMIXLoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)
    error("Ice-temp data not provided in ITMIX.  Use different loader, probably syntetic")
end
function load!(ds::DataSet{TerminusFluxData, ITMIXLoader}, dem, glaciermask, iscalving,
               pp::Phys)
    error("Calving flux not provided in ITMIX.  Use different loader, probably synthetic")
end

#############################
# Helper

# Each of this folders contains a combination of the following files:
#
# - 01_margin_GGGGGG_YYYY_UTMZZZ.csv
# - 02_surface_GGGGGG_YYYY_UTMZZZ.asc  #
# - 03_RES_GGGGGG.txt                  # GPR file
# - 04_mb_GGGGGG_XXXX-YYYY_UTMZZZ.asc
# - 05_dhdt_GGGGGG_XXXX-YYYY_UTMZZZ.asc
# - 06_speed_GGGGGG_XXXX-YYYY_UTMZZZ.asc
# - 07_azimuth_GGGGGG_XXXX-YYYY_UTMZZZ.asc
# - 08_velocity_GGGGGG_XXXX-YYYY_UTMZZZ.txt
#
# where
#
# GGGGGG = name of the glacier (e.g. "Brewster")
# YYYY   = year the data are referring to (e.g. "2003")
# ZZZZ   = UTM zone of the coordinate system in which the data are projected (e.g. "32") - negative values indicate the southern hemisphere.
# Files 04-08 are mean values over the period "XXXX-YYYY" (e.g. "1997-2009")
#
# The file content is as follows:
# 01: Margin of the glacier (m).
# 02: DEM of the glacier (m a.s.l.) - Note that the DEM might be cropped with the glacier margin or not
# 04: Spatially distributed average (surface) mass balance of the glacier (m w.e./a)
# 05: Spatially distributed average ice thickness change rate (m/a)
# 06: Spatially distributed average surface flow speed of the glacier (m/a)
# 07: Spatially distributed average azimuth (="direction") of the surface flow (clockwise degrees from north)
# 08: Point observations of surface flow velocity (comprising speed and azimuth)
#
# .txt-files are plain text files.
# .csv-files are comma-separated-value files.
# .asc-files are ASCII-grid files (http://en.wikipedia.org/wiki/Esri_grid)
function parse_ITMIX_filename(fln::ZipFile.ReadableFile)
    error("zip files not supported currently")
    dir,fl = splitdir(fln.name)
    if contains(dir,"shapefiles")
        return nothing
    else
        return parse_ITMIX_filename(fl)
    end
end
function parse_ITMIX_filename(fln::String)
    fln=="" && return nothing
    fl, ext = splitext(splitdir(fln)[2])
    splits = find(c->c=='_', fl)
    length(splits)<2 && return nothing
    !(ext in [".txt",".asc",".shp"]) && return nothing # drop the ".csv" outlines in favor of shp
    nr = fl[1:splits[1]-1]
    !all(isnumber, nr) && return nothing
    nr = parse(Int,nr)
    dk = fl[splits[1]+1:splits[2]-1]
    # if dk in ["margin","surface"]
    #     name = fl[splits[2]+1:splits[3]-1]
    #     DK = dk=="margin" ? Outline : DEM
    #     y = parse(Int, fl[splits[3]+1:splits[4]-1])
    #     date_range = (Date(y),Date(y))
    # else
    DK = dk_dict[dk]
    if dk=="RES"
        name = fl[splits[2]+1:end]
        date_range = (Date(),Date())
    else
        name = fl[splits[2]+1:splits[3]-1]
        sp = split(fl[splits[3]+1:splits[4]-1],'-')
        if length(sp)==1
            y1 = parse(Int,sp[1])
            y2 = y1
        else
            y1,y2 = (parse(Int,s) for s in sp)
        end
        date_range = (Date(y1),Date(y2))
    end
    if dk!="RES"
        if name=="Unteraar"
            # Unteraar is not in proper UTM projection!
            zone = fl[splits[end]+1:end]
            proj = "close to UTM zone $zone, but no cigar"
        else
            zone = fl[splits[end]+1:end]
            proj = "+proj=utm +zone=$zone +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
        end
    else
        proj = "unknown"
    end
    return DK, fln, nr, name, date_range, proj
end

"Read a agr file"
_read(name, dirs, num) = read_agr(_getfile(dirs, name, num), F, NA=FILL)


function generate_DEM_masks(gid::ITMIXGlacier, dem, outline, sealevel_cutoff=10.0)
    glaciermask, _ = maskit(dem, outline)
    # also DEM holes from glaciermask
    glaciermask[dem.v.==FILL] = false
    # remove points below sea level
    glaciermask[dem.v.<0] = false

    remove_slivers = true
    landmask = generate_landmask(glaciermask, dem, sealevel_cutoff, remove_slivers)
    return glaciermask, landmask
end

######################
# Special processing (aka hacks)
######################

function special_processing!(gid::ITMIXGlacier, outline, dem, glaciermask, landmask,
                             iv, ivmask, thickness, bdot, dhdt, fsl, temp,
                             pl::LoadPara, pp::Phys, pm::MPara, pn::Num)
    if getname(gid)=="Austfonna"
        # no easy way of figuring land vs sea out for Austfonna, set all to sea.
        landmask[:] = false
        ivmask[iv.vals.v.<=0] = false
    elseif getname(gid)=="Tasman" || getname(gid)=="Devon" || getname(gid)=="Unteraar"
        # not masked at all in the input data
        ivmask[iv.vals.v.<=0] = false
    end
    return nothing
end


###########################
## ITMIX2 Utilities
###########################
module ITMIX2
import ..BITEModel; const BM = BITEModel
import VAWTools
"Number of experiments per glacier"
const nexps = 16
using DataFrames, FileIO, CSVFiles
int2expnr(i) = i>9 ? Symbol("exp$i") : Symbol("exp0$i")

"""
    load_experiment_table(fln="data/ITMIX2/ITMIX2_experiments_key.txt")

Load ITMIX2_experiments_key.txt as a DataFrame
"""
function load_experiment_table(fln)
    tab = DataFrame(load(File(format"CSV", fln),
                         spacedelim=true,
                         header_exists=true,
                         skiplines_begin=11,
                         type_detect_rows=512,
                         #colparsers=[String,Int, fill(x->Bool(parse(Int,x)),16)..., Int]
                         ))
    tab[:glacier] = [BM.ITMIXGlacier(Symbol(gl), 2) for gl in tab[:glacier]]
    for i=1:nexps
        tab[int2expnr(i)] = Bool.(tab[int2expnr(i)])
    end
    # add exp00
    tab[:exp00] = tab[:exp01].*false
    return tab
end

"Get the priority of a glacier"
getpriority(gid) = getexps(gid)[:prio][1]
getpriority() = Dict(((gln, getexps(gln)[:prio][1]) for gln in BM.ITMIXnames))

"Get experiments for a glacier `gl`"
getexps(gid::BM.ITMIXGlacier) = DataFrames.filter(row->row[:glacier]==gid, tab)
getexps(name::Union{Symbol,String}) = getexps(BM.ITMIXGlacier(name,2))

"Get fitting radar-lines of experiment `i` for a glacier `gl`."
get_rlines(gid::BM.ITMIXGlacier, i) = i==0 ? Int[] : find(getexps(gid)[int2expnr(i)])
get_rlines(name::Union{Symbol,String}, i) =  i==0 ? Int[] : get_rlines(BM.ITMIXGlacier(name,2), i)

fln = joinpath(@__DIR__, "../../scripts/data/ITMIX2/ITMIX2_experiments_key.txt")
"The ITMIX2 experiment table"
const tab = isfile(fln) ? load_experiment_table(fln) : nothing

"The ITMIX2 glaciers sorted according to their priorities"
const priorities = tab==nothing ? nothing : [[BM.ITMIXGlacier(gln,2) for gln in keys(getpriority()) if getpriority()[gln]==i] for i=1:4]


# # iteration over the ITMIX 2 experiments
# """
#     experiments([gid])

# Return iterator for ITMIX2 experiments.

# ```
# for (gid, i, radarlines) in experiments([gid])
#    ...
# end
# ```
# """
# function experiments end

# experiments(gid) = ITMIX2Iterator(gid, getexps(gid))
# struct ITMIX2Iterator
#     gid::BM.ITMIXGlacier
# end
# Base.start(::ITMIX2Iterator) = 1
# Base.next(ii::ITMIX2Iterator, state) = ( (gid, int2expnr(state), getexp(ii.gid, state)), state+1)
# Base.done(ii::ITMIX2Iterator, state) = state==17

# experiments() = ITMIX2Iterator_All()
# struct ITMIX2Iterator_All end
# Base.start(::ITMIX2Iterator_All) = (1, BM.ITMIXv2_glaciers[1])
# function Base.next(ii::ITMIX2Iterator_All, state)
#     expnr, gid = state
#     item = (gid, int2expnr(expnr), getexp(gid, expnr))
#     state = if expnr==16 && gid!=BM.ITMIXv2_glaciers[end]
#         (1, BM.ITMIXv2_glaciers[findfirst(BM.ITMIXv2_glaciers, gid)+1])
#     else
#         (expnr+1, gid)
#     end

#     return item, state
# end
# Base.done(ii::ITMIX2Iterator_All, state) = state[1]==17

"""
    make_ITMIX2_glacier(gl, expnr::Number, rm_iv=false)
    make_ITMIX2_glacier(gl, tracks2keep, rm_iv=false)

Take the glacier with all radarlines loaded and return one which only has the
lines of that experiment loaded.  The second interface allows to specify which
tracks to keep and ditch the rest.

Also has the extra feature to make more experiments:
- exp0: no radar lines
- exp 17-32: no IV, otherwise as 1.16

Above two only make sense (in the sense that there is something to fit)
for glaciers with IV data.
"""
make_ITMIX2_glacier(gl, expnr::Number, rm_iv=false) =
    make_ITMIX2_glacier(gl, get_rlines(gl.gid, expnr), rm_iv, expnr)
function make_ITMIX2_glacier(gl, tracks, rm_iv=false, expnr=-1)
    new_h = if length(tracks)>0  # ITMIX2 experiments
        BM.SomeData(gl.h, vals=VAWTools.Traj(gl.h.vals, tracks))
    else # remove all radar-lines: name it exp00.  This is the ITMIX1 run.
        new_h = BM.SomeData(nothing)
    end
    misc = merge(gl.misc, Dict{Symbol,Any}(:expnr=>expnr, :rm_iv=>rm_iv))
    return if rm_iv
        # check that there is IV
        gl.iv isa BM.SomeData{Void} && error("Cannot remove IV from glacier $(BM.getname(gl)) as it has no IV.")
        # remove IV data
        BM.Glacier(gl, h=new_h, iv=BM.SomeData(nothing), misc=misc)
    else
        # normal run
        BM.Glacier(gl, h=new_h, misc=misc)
    end
end

##########################
## Saving
##########################
# needs to be saved on the same grid as the input data
#=
ITMIX2 - instructions for results submission

1) Results format and content
- Submitted results shall be distributed estimates of the ice thickness.
- The unit for the ice thickness is "meters".
- Results are to be submitted as Esri grids, i.e. in the same format of the input data.
- All results shall have same resolution, cellsize, extent, and coordinate system as the input grids.

2) Results nomenclature
- All results shall be contained in a single directory, called "results_[Yourname]" (e.g. "results_Farinotti")
- The directory shall be compressed (.zip) before submission. The name of the uploaded file will thus be: "results_[Yourname].zip"
- The individual files in the folder shall be named
   "thickness_[NameOfTheGlacier]_[Yourname]_exp[XX].asc"
  where
   [NameOfTheGlacier] = Name of the glacier, as defined by the input files (e.g. "NorthGlacier")
   [Yourname] = Your family name (e.g. "Farinotti")
   [XX] = Number of the experiment, as given in the "experiments_key"-file, the result is referring to (e.g. "02")
  A valid filename will thus look like, e.g., "thickness_NorthGlacier_Farinotti_exp02.asc".
- For a test case to be labelled as "considered", results for all 16 experiments will need to be submitted.
  Since there are three test cases defined as "compulsory", the submitted directory will thus contain 3*16=48 files at least.

3) Results upload
- To uplod your results, simply drag-and-drop your .zip-repository to the following address:
 https://drive.switch.ch/index.php/s/hovc5mJjSztoLVd
- When uploaded, please send a notification email to: daniel.farinotti@ethz.ch
=#

"""
    filepaths(gl::BM.Glacier, pl, expnr, rm_iv, dir)

Returns saving filenames including their full path.
The path is `joinpath(pl.data_root, dir).

Returns file names for h, h-sigma, iv, iv-sigma, and jld.
"""
function filepaths(gl::BM.Glacier, pl, expnr, rm_iv, dir)
    expnrs = expnr>9 ? Symbol("exp$expnr") : Symbol("exp0$expnr")
    glaciername = BM.getname(gl)
    fl = rm_iv ?
        joinpath(pl.data_root, dir, "thickness_$(glaciername)_Werder_$(expnrs)_rmIV") :
        joinpath(pl.data_root, dir, "thickness_$(glaciername)_Werder_$(expnrs)")
    # make the directory if it does not exist
    !isdir(pl.data_root) && error("Data-root dir does not exist: $(pl.data_root)")
    !isdir(joinpath(pl.data_root, dir)) && mkpath(joinpath(pl.data_root, dir))
    return fl*".asc", fl*"-sigma.asc", fl*"-iv.asc", fl*"-iv-sigma.asc", fl*".jld"
end


"""
    save_itmix2_run(hs2d, hs2d_sigma, gl, pl, iv2d=nothing, iv2d_sigma=nothing)

Saves the thickness map according to ITMIX2 specs.
"""
function save_itmix2_run(hs2d, hs2d_sigma, gl, pl, ivs2d=nothing, ivs2d_sigma=nothing;
                         dir="results_Werder")
    expnr, rm_iv = gl.misc[:expnr], gl.misc[:rm_iv]
    fln, flns, fln_iv, fln_ivs, fln_jld = filepaths(gl, pl, expnr, rm_iv, dir)

    downsample = get(pl.dataset_opts[BM.DataKind], :grid_downsample_step, 1)
    hs2d, hs2d_sigma = copy(hs2d), copy(hs2d_sigma)
    hs2d[.!gl.glaciermask] = BM.FILL
    hs2d_sigma[.!gl.glaciermask] = BM.FILL
    h = VAWTools.Gridded(gl.dem.x, gl.dem.y, hs2d)
    hs = VAWTools.Gridded(gl.dem.x, gl.dem.y, hs2d_sigma)
    if ivs2d!=nothing
        ivs2d, ivs2d_sigma = copy(ivs2d), copy(ivs2d_sigma)
        ivs2d[.!gl.glaciermask] = BM.FILL
        ivs2d_sigma[.!gl.glaciermask] = BM.FILL
        iv = VAWTools.Gridded(gl.dem.x, gl.dem.y, ivs2d)
        ivs = VAWTools.Gridded(gl.dem.x, gl.dem.y, ivs2d_sigma)
    end
    if downsample!=1
        # up-sample output
        dem_orig = VAWTools.read_agr(pl.dataset_files[BM.DEMData][1])
        h = VAWTools.upsample(h, dem_orig.x, dem_orig.y)
        hs = VAWTools.upsample(hs, dem_orig.x, dem_orig.y)
        if ivs2d!=nothing
            iv = VAWTools.upsample(iv, dem_orig.x, dem_orig.y)
            ivs = VAWTools.upsample(ivs, dem_orig.x, dem_orig.y)
        end
    end
    VAWTools.write_agr(h, fln, NA_g=BM.FILL, NA_agr=-9999.0)
    VAWTools.write_agr(hs, flns, NA_g=BM.FILL, NA_agr=-9999.0)
    if ivs2d!=nothing
        VAWTools.write_agr(iv, fln_iv, NA_g=BM.FILL, NA_agr=-9999.0)
        VAWTools.write_agr(ivs, fln_ivs, NA_g=BM.FILL, NA_agr=-9999.0)
    end
    # TODO make sub-folders with computer-name-date
    return fln, flns, fln_iv, fln_ivs, fln_jld
end
end #module ITMIX2
