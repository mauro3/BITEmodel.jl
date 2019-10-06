# misc CCI-file loaders

# function load!(ds::DataSet{IVData, CCILoader}, dem, glaciermask, target_proj, outline, cropbox,
#                pp::Phys)::Tuple{SomeData,Matrix{Bool}}
#     # TODO: allow other IV sets
#     @unpack files, opts, daterange = ds

#     zone = split(split(target_proj)[2],'=')[2]
#     flx = files[1]*"utm$zone-x.asc.gz"
#     fly = files[1]*"utm$zone-y.asc.gz"

#     dsx = reconstruct(ds, files=[flx])
#     ivx = load_asciigrid!(dsx, :ivx, cropbox)
#     dsy = reconstruct(ds, files=[fly])
#     ivy = load_asciigrid!(dsy, :ivy, cropbox)
#     fill = ivx.v.==FILL
#     tmp = sqrt.(ivx.v.^2 + ivy.v.^2)
#     # m/d -> m/y
#     tmp .*= 365
#     tmp[fill] = FILL
#     iv = reconstruct(ivx, v=tmp)

#     # # check projection -> now using re-projected files
#     # proj = ""
#     # pfl = splitext(splitext(files[1])[1])[1]*".prj"
#     # if isfile(pfl)
#     #     for st in split(readstring(`gdalsrsinfo iv-2017-y.prj`), '\n')
#     #         if startswith(st, "PROJ.4 :")
#     #             proj = strip(split(st, '\'')[2])
#     #             @assert startswith(proj, "+proj=")
#     #         end
#     #     end
#     #     if proj!=target_proj

#     # end

#     # according to *.xml sheet
#     # <meanE units="m/d">0.00</meanE>
#     # <rmseE units="m/d">0.01</rmseE>
#     # <meanN units="m/d">0.00</meanN>
#     # <rmseN units="m/d">0.02</rmseN>
#     @get sigma = opts | sqrt((0.01*365)^2 + (0.02*365)^2)

#     # get/make masks
#     ivmask = get_cache!(ds, outline.gid, :ivmask,
#                         () -> maskit(iv, outline)[1],
#                         Bool, size(iv))

#     # fill in any gaps inside the glacier
#     fillgaps!(iv, ivmask, mask_unfilled_holes=true)
#     return SomeData(vals=iv, sigma=sigma), ivmask
# end

function load!(ds::DataSet{IVData, CCILoader}, dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::Tuple{SomeData,Matrix{Bool}}
    @unpack files, opts, daterange = ds
    # figure out utm-zone
    zone = split(split(target_proj)[2],'=')[2]
    fln = files[1]*"-utm$zone.tif"
    ds = reconstruct(ds, files=[fln])
    iv = load_geotiff!(ds, Symbol("iv-$zone"), cropbox)

    @assert !any(isnan.(iv.v)) "Loaded IV has NaNs"

    # get/make masks
    ivmask = get_cache!(ds, outline.gid, :ivmask,
                        () -> maskit(iv, outline)[1],
                        Bool, size(iv))

    # I don't think this is a good idea for the rough data we got # fill in any gaps inside the glacier
    # fillgaps!(iv, ivmask, mask_unfilled_holes=true)

    # mask it where it is within threshold
    th = get(opts, :mask_threshold, (0.0, Inf))
    ivmask[.!(th[1].< iv.v .<th[2])] = false
    # mask where FILL value
    ivmask[iv.v.==FILL] = false

    if any(size(iv.v).<2)
        # if there is not enough IV data, don't use it
        return SomeData(nothing), Matrix{Bool}(0,0)
    else
        return SomeData(vals=iv, sigma=opts[:sigma]), ivmask
    end
end

function load!(ds::DataSet{DhdtData, CCILoader},
               bdot, terminus_flux,
               dem,
               glaciermask, target_proj, outline, cropbox,
               iscalving,
               pp::Phys)::SomeData
    @unpack files, opts, daterange = ds
    zone = split(split(target_proj)[2],'=')[2]
    fln = files[1]*"-utm$zone.tif"
    ds = reconstruct(ds, files=[fln])
    dhdt = SomeData(vals=load_geotiff!(ds, Symbol("dhdt-$zone"), cropbox),
                    sigma=ds.opts[:sigma], vals_prior=0.0)
    if length(dhdt.vals.v)==0
        # use SyntheticBenchLoader
        ds = DataSet{DhdtData, SyntheticBenchLoader}(ds)
        dhdt = load!(ds,bdot, terminus_flux,
                     dem,
                     glaciermask, target_proj, outline, cropbox,
                     iscalving,
                     pp)
    end
    return dhdt
end
