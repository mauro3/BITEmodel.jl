"Gets parts of the RGI number as zero-padded stings"
function get_rgi_parts(gid::GlacierID)
    id5 = VAWTools.int2str5(getrgi_nr(gid))
    region = string(getregion(gid))
    region_nr = VAWTools.int2str2(dir2reg[Symbol(region)])
    id5, region, region_nr
end

"""
Reads the attribute file distributed with the RGI dataset.

Available attributes are
"RGIId" "GLIMSId" "BgnDate" "EndDate" "CenLon" "CenLat" "O1Region" "O2Region" 8
"Area" "Zmin" "Zmax" "Zmed" "Slope" "Aspect" "Lmax" "Status" "Connect" "Form" 18
"TermType" "Surging" "Linkages" "Name"
"""
function read_attrib_file(region, glnrs=:all;
                          dir = joinpath(normpath(dirname(@__DIR__)), "../scripts/data/rgi-6.0-huss/attribs/"))
    fls = readdir(dir)
    f = ""
    for f in fls
        if startswith(f, VAWTools.int2str2(region))
            break
        end
    end
    a,h = readcsv(dir*f, header=true)

    glnrs_all = getrgi_nr.(a[:, 1])
    # figure out indices of glnrs in a
    inds = Int[]
    if glnrs==:all
        inds = 1:size(a,1)
        glnrs = getrgi_nr.(a[inds, 1])
    elseif glnrs isa Range
        for (i,g) in enumerate(glnrs_all)
            if g in glnrs
                push!(inds, i)
            end
        end
        glnrs = collect(glnrs)
    else # assume g
        for g in glnrs
            i = searchsorted(glnrs_all, g)
            @assert length(i)==1 "Input gl-nr $glnrs not found in the glacier numbers of RGI region $region"
            push!(inds, i[1])
        end
    end
    ## Ratio of glaciers surveyed in a span of years:
    # @show count(1999.<=a[:,3]./10000.<=2010)./length(a[:,3])
    lon = F.(a[inds, 5])
    lat = F.(a[inds, 6])
    subregion = Int.(a[inds, 8])
    area = F.(a[inds, 9])*1e6 # m^2
    zmin = F.(a[inds, 10])
    zmax = F.(a[inds, 11])
    aspect = F.(a[inds, 14])
    status = Int.(a[inds, 16])
    form = Int.(a[inds, 18])
    termtype = Int.(a[inds, 19])
    surging = Int.(a[inds, 20])

    name = String.(a[inds, end])
    for (i,n) in enumerate(name)
        name[i] = isvalid(n) ? n : "" # the attribute file for Svalbard has some bad encodings (although I fixed it).  Drop those to avoid errors later.
    end
    return glnrs, lon, lat, subregion, area, zmin, zmax, aspect, status, form, termtype, surging, name
end
read_attrib_file(region, glnrs::Number;
                 dir = joinpath(normpath(dirname(@__DIR__)), "../scripts/data/rgi-6.0-huss/attribs/")) =
                     ([a[1] for a in read_attrib_file(region, [glnrs], dir=dir)]...,)


"""
    has_radar(gid::Union{String,RGIGlacier})

Figures out whether a RGI glacier has radar measurements.
"""
function has_radar(gid::RGIGlacier,
                   glathida_dir = joinpath(normpath(dirname(@__DIR__)), "../scripts/data/GlaThiDa_3.0/mhuss"))
    fls = readdir(glathida_dir)
    targetname = VAWTools.int2str2(getrgi_region(gid)) * "_" *
        VAWTools.int2str5(getrgi_nr(gid))
    for f in fls
        if startswith(f, targetname)
            return true
        end
    end
    return false
end
has_radar(rgi::String) = has_radar(RGIGlacier(rgi))

"Get all glacier of a region"
function get_glaciers(region;
                      dir = joinpath(normpath(dirname(@__DIR__)), "../scripts/data/rgi-6.0-huss/attribs/"))
    read_attrib_file(region, dir=dir)[1]
end

"""
    get_glaciers_with_radar(rgi_region::Integer, glathida_dir=usual_dir; use_only_good_ones=true)

Get all glaciers of a RGI region with good radar data in directory glathida_dir.  No effort is made
to check that the returned rgi-numbers are correct.

If use_only_good_ones==false then return all nrs in glathida_dir.
"""
function get_glaciers_with_radar(rgi_region::Integer,
                                 glathida_dir = joinpath(normpath(dirname(@__DIR__)), "../scripts/data/GlaThiDa_3.0/mhuss");
                                 use_only_good_ones=true)
    fls = readdir(glathida_dir)
    rstr = VAWTools.int2str2(rgi_region)

    # if there is a list of glacier with good data use that too
    good_ones = 1:10^6
    if use_only_good_ones
        open(joinpath(glathida_dir, "glaciers_with_useful_data")) do io
            l = readline(io)
            while !isempty(l)
                if !startswith(l, "#")
                    ll = parse.(Int, split(l, ", "))
                    if ll[1]==rgi_region
                        good_ones = ll[2:end]
                        break
                    end
                end
                l = readline(io)
            end
        end
    end

    glsnrs = Int[]
    for f in fls
        if startswith(f, rstr)
            nr = parse(Int, split(splitext(f)[1], '_')[2])
            if nr in good_ones
                push!(glsnrs, nr)
            end
        end
    end

    return glsnrs
end

"""
    LoadPara(gid::RGIGlacier; pl_kwargs...)

For sure sets up to loads the RGI outline.  The rest can be tweaked but this
tries to provide good defaults:
- DEM from Matthias RGI-compatible DEM collection
- bdot from GloGEM
- fsl, temp, TerminusFluxData from Synthetic-loader
- thickeness from GlaThiDa
- dhdt from various sources, e.g. CCI
- iv from various sources, e.g. CCI
"""
function LoadPara(gid::RGIGlacier; pl_kwargs...)
    glacier = getid(gid)
    pl_kwargs = deepcopy(Dict(pl_kwargs))
    use_synthetic_dhdt = get!(Dict(pl_kwargs), :use_synthetic_dhdt, false)
    delete!(pl_kwargs, :use_synthetic_dhdt) # otherwise below errors:
    pl = LoadPara(;gid=gid, data_root=joinpath(LoadPara(gid=gid).data_root,"rgi-6.0-huss"), pl_kwargs...)
    # TODO: add this as opts[IVData][:sigma] = ... below
    for (k,v) in (ThicknessData=>10.0, # m
                  BdotData=>1.0,  # m/a, what Matthias says about GloGEM
                  DhdtData=>1.0,  # m/a, e.g. Nuth et al. 2010, Table 1
                  IVData=>10.0)  # m/a
        pl.dataset_opts[k][:sigma]=v
    end
    updated_pm_defaults = Dict(
        :bandsize => 30, # lowered automatically for smaller glaciers
        :cap_at_floatation=>true,
        :dist_exp => 0.34) # this value is suggested by fitting the two glacier on Svalbard with thickness measurements
    # hacks for Baltoro and Siachen glacier, as bands are not continuous otherwise
    if gid in [RGIGlacier("RGI60-14.06794"), RGIGlacier("RGI60-14.07524")]
        updated_pm_defaults[:bandsize] = 50
    end

    # update pm but keeping any explicitly set (via kwargs) options
    merge!((a,b)->a, pl.dataset_opts[ParaData][:pm], updated_pm_defaults)

    # Set loaders & files
    id5, region_str, region_nr2 = get_rgi_parts(gid)
    region_nr = parse(Int, region_nr2)
    region_strl = lowercase(region_str)
    pl = deepcopy(pl)

    dl = pl.dataset_LOADERS
    fl = pl.dataset_files
    opts = pl.dataset_opts
    # Set the LOADERS
    dl[OutlineData] = RGILoader
    fl[OutlineData] = [joinpath(pl.data_root, "outlines/xyzn/$region_strl/$id5.xyn"), # UTM versions have problems -> use long-lat
                       joinpath(pl.data_root, "grids/$region_strl/gl/gl_$id5.grid")] # the mask contains the projection
    dl[ParaData] = RGILoader

    dl[DEMData] = RGILoader
    fl[DEMData] = [joinpath(pl.data_root, "grids/$region_strl/dem/dem_$id5.grid"),
                   joinpath(pl.data_root, "grids/$region_strl/gl/gl_$id5.grid"), # the mask
                   joinpath(pl.data_root, "outlines/COMMON_BOUNDARY/$region_strl/$id5.dat"), # the COMMON_BOUNDARY
                   joinpath(pl.data_root, "outlines/xyzn/$region_strl/"), # path to outlines of other glaciers
                   ]
    opts[DEMData][:cutoff] = 0.0
    opts[DEMData][:sealevel_cutoff] = 120.0
    opts[DEMData][:cropbox] = :from_DEM # as Matthias has already cropped the DEMs
    opts[DEMData][:masksource] = :from_file # load mask from file
    if region_str=="Svalbard"
        # there are a few slivers where the dem is bad on three glaciers
        opts[DEMData][:mask_unfilled_holes] = true
    end

    #### End of pure RGI stuff.  Now add sensible loader for the rest.

    if parse(Int, id5) in get_glaciers_with_radar(region_nr)
        dl[ThicknessData] = GlaThiDaLoader
        fl[ThicknessData] = [joinpath(pl.data_root, "../GlaThiDa_3.0/mhuss/$(region_nr2)_$id5.dat"),
                             joinpath(pl.data_root, "../GlaThiDa_3.0/mhuss/AVG_$(region_nr2).dat")]
    else
        dl[ThicknessData] = VoidLoader
    end
    # If only an average thickness is available then set the relative mean thickness error:
    opts[ThicknessData][:sigma_mean_h_rel] = 0.3

    dl[FslData] = SyntheticBenchLoader
    dl[TempData] = SyntheticBenchLoader
    dl[TerminusFluxData] = SyntheticBenchLoader

    dl[BdotData] = GloGEMLoader
    fl[BdotData] = [joinpath(pl.data_root, "../GloGEM/$(ucfirst(region_str))/belev_$id5.dat")]

    ### This is where it gets messier as the data is heterogeneous

    # the proper files are figured out later when the projection is known:
    dl[IVData] = CCILoader
    fl[IVData] = [joinpath(pl.data_root, "../CCI/iv_rgi$(region_nr2)")] # add utm zone to that, eg "-utm36.tif"
    if region_nr in [3, 4, 14]
        opts[IVData][:mask_threshold] = (30,500)
    elseif region_nr == 7
        opts[IVData][:mask_threshold] = (50,500)
    elseif region_nr == 11
        opts[IVData][:mask_threshold] = (5,500)
    else
        error("Figure out what a good lower bound on IV is for this region & data")
    end

    if !use_synthetic_dhdt && region_nr in [7,11,14]
        get!(dl, DhdtData, CCILoader)
        fl[DhdtData] = [joinpath(pl.data_root, "../CCI/dhdt_rgi$(region_nr2)")] # add utm zone to that, eg "-utm36.tif"
    else
        dl[DhdtData] = SyntheticBenchLoader
    end
    # TODO
    # - add ice sat dhdt for svalbard

    return pl
end

# Outline
function load!(ds::DataSet{OutlineData, RGILoader}, gid)
    iscalving = gid.attrs.termtype in [1,2,3,5] # this might get modified later as in-accurate
    isicecap =  gid.attrs.form==1

    # projection
    proj = VAWTools.get_utm_asciigrid(ds.files[2])

    # outlines & outcrops
    outline, splits = VAWTools.concat_poly(read_xyn(ds.files[1]))
    outline = VAWTools.transform_proj(outline', "+proj=latlong +datum=WGS84", proj)'

    # # TODO see whether there is a COMMON_BOUNDARY file (but in UTM...)
    # cbf = replace(replace(ds.files[1], "xyzn", "COMMON_BOUNDARY"), ".xyn", ".dat")
    # if isfile(cbf)
    #     cb = readdlm(cbf)
    # end

    return Outline(gid, outline, splits), iscalving::Bool, isicecap::Bool, proj::String
end

# DEM: Arctic DEM, SRTM etc, by Matthias
function load!(ds::DataSet{DEMData, RGILoader}, outline::Outline, pp::Phys
               )::Tuple{Gridded, Matrix{Bool}, Matrix{Bool}, Matrix{F}, F, Tuple, String}
    (dem, glaciermask, landmask,
     dist_land, dist_land_maxdist, cropbox,
     dem_proj) = load!(convert(DataSet{DEMData, :ASCIIGridLoader}, ds), outline, pp)
    dem.v[dem.v.<0] = NaN

    ## Update landmask with other glaciers which are attached
    if isfile(ds.files[3])
        mask = copy(glaciermask)
        # Other glaciers are from the COMMON_BOUNDARY files
        otherglaciers = round.(Int, unique(readdlm(ds.files[3])[:,3]))
        for og in otherglaciers
            id5 = VAWTools.int2str5(og)
            ogf = joinpath(ds.files[4], "$id5.xyn")
            xy, splits = VAWTools.concat_poly(read_xyn(ogf))
            xy = VAWTools.transform_proj(xy', "+proj=latlong +datum=WGS84", dem_proj)'
            m = maskit(dem, xy)[1]
            mask = mask .| m
        end
        # TODO Figure out something clever to do with calving glaciers, the sea
        #      and the landmask...

        # Remove slivers of land:
        sealevel_cutoff = get(ds.opts, :sealevel_cutoff, 10)
        landmask = generate_landmask(mask, dem, sealevel_cutoff, true)

        # recalculate Distance to land
        dist_land_maxdist = ds.opts[:dist_land_maxdist]
        dist_land = dist_to_mask(dem, landmask, glaciermask, dist_land_maxdist)
    end

    return dem, glaciermask, landmask,
           dist_land, dist_land_maxdist, cropbox,
           dem_proj
end

#################
# load all, forward and inverse stuff for RGI glacier

init(rgi::String, fit_vars, fit_target, sigma_of_model, runtyp; plkws...) =
    init(RGIGlacier(rgi), fit_vars, fit_target, sigma_of_model, runtyp; plkws...)
init(nr, region, fit_vars, fit_target, sigma_of_model, runtyp; plkws...) =
    init(RGIGlacier(nr,region), fit_vars, fit_target, sigma_of_model, runtyp; plkws...)

function init(gid::RGIGlacier, fit_vars, fit_target, sigma_of_model, runtyp;
              fit_sigma=false,
              plkws...)
    # load glacier
    print("Time to load glacier:  ")
    @time begin
        pl = LoadPara(gid; plkws...)
        dt, pl = make_datatable(gid, pl);
        dt.para.opts[:pm][:sigma_h_model] = sigma_of_model.sigma_h_model
        dt.para.opts[:pm][:sigma_iv_model] = sigma_of_model.sigma_iv_model

        gl,pp,pm,pn,pl = load_glacier(dt,pl)
        gb,pm = make_bands(gl, pp, pm, pn)
    end;

    # forward model
    Base.Test.@inferred fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

    # Inverse
    th0d = theta0_dict_defaults(pm, fit_vars...;
                                with_sigma=fit_sigma)

    (theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) =
         init_inverse(gb, pp, pm, pn,
                      runtyp=Symbol(runtyp),
                      theta0_dict=th0d,
                      fit_target=fit_target
                      )
    return (gid, gl, gb, pp, pm, pn, pl, th0d,
            theta0, logposterior, logposterior1d, logprior, logprior_debug,
            loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
            pmcmc_defaults, fit_target)
end
