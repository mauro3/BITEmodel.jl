# Generic loaders for
# - ascii grid
# - shapefile

import Shapefile, GeoInterface

"""
$(SIGNATURES)

The `symb` is used in caching.

Function which:
- loads an ASCII grid
- crop it (if cropbox!=nothing)
- downsamples it (if requested in ds.opts).
"""
function load_asciigrid!(ds::DataSet, symb::Symbol, cropbox=nothing)

    @unpack files, opts, daterange = ds
    @assert length(files)==1 "Not supported to have several AsciiGrid files with the generic loader"

    out = get_cache_ram!(ds, symb, @open_magic_fn read_agr(fl, F, NA=FILL))

    # Crop
    if cropbox!=nothing
        _, out = crop(cropbox, out)
    end

    # Downsample
    grid_downsample_step = get(opts, :grid_downsample_step, 1)
    if grid_downsample_step>1
        averagemask = out.v.!=FILL
        out = downsample(out, grid_downsample_step, 1, true, averagemask)
    end
    return out
end

"""
$(SIGNATURES)

The `symb` is used in caching.

Function which:
- loads an GeoTiff grid
- crop it (if cropbox!=nothing)
- downsamples it (if requested in ds.opts).
"""
function load_geotiff!(ds::DataSet, symb::Symbol, cropbox=nothing)

    @unpack files, opts, daterange = ds
    @assert length(files)==1 "Not supported to have several GeoTiff files with the generic loader"

    out = get_cache_ram!(ds, symb, () -> VAWTools.read_geotiff(files[1], F, NA=FILL))

    # Crop
    if cropbox!=nothing
        _, out = crop(cropbox, out)
    end

    # Downsample
    grid_downsample_step = get(opts, :grid_downsample_step, 1)
    if grid_downsample_step>1
        averagemask = out.v.!=FILL
        out = downsample(out, grid_downsample_step, 1, true, averagemask)
    end
    return out
end


"""
$(SIGNATURES)


Function which loads a shapefile outline.
"""
function load_shapefile!(ds::DataSet{OutlineData}, symb)
    @unpack files, opts, daterange = ds

    @assert length(files)==1 "Cannot read several (or zero) shapefiles: $files"
    @assert endswith(ds.files[1], ".shp")

    handle = get_cache_ram!(ds, symb, () -> open(files[1], "r") do io
                                              read(io, Shapefile.Handle)
                            end)
    coords = []
    if eltype(handle.shapes) <: Shapefile.Polygon
        for sh in handle.shapes
            push!(coords, hcat(GeoInterface.coordinates(sh)[1][1]...))
        end
    elseif eltype(handle.shapes) <: Shapefile.Polyline
        for sh in handle.shapes
            for cc in GeoInterface.coordinates(sh)
                push!(coords, hcat(cc...))
            end
        end
    elseif eltype(handle.shapes) <: Shapefile.PolygonZ
        for sh in handle.shapes
            push!(coords, hcat(GeoInterface.coordinates(sh.points)...))
        end
    else
        error("Don't know how to handle $(eltype(handle.shapes))")
    end

    outline, splits = VAWTools.concat_poly(coords, close_poly=true)
    return outline::Matrix{Float64}, splits
end


"Generic ASCIIGrid loader."
function load!(ds::DataSet{DK, :ASCIIGridLoader},
               dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::SomeData where DK
    out = load_asciigrid!(ds, datakind2varname[DK], cropbox)
    # What about scaling? Say m/s -> m/a?
    return SomeData(vals=out, sigma=ds.opts[:sigma])
end
function load!(ds::DataSet{DEMData, :ASCIIGridLoader},
               outline, pp::Phys)
    @unpack opts = ds

    dem = if length(ds.files)==1
        load_asciigrid!(ds, :dem, nothing)
    else
        load_asciigrid!(DataSet(ds, files=ds.files[1:1]), :dem, nothing)
    end

    # DEM specific processing:

    # Projection:
    # Note that the DEM cannot be re-projected at the moment as that would also mean
    # defining another grid.  It is also assumed that the outline is in the same
    # projection.
    dem_proj = get(opts, :proj, VAWTools.get_utm_asciigrid(ds.files[1]))

    # Use a grid other than the one of the DEM on file.
    # Note that the grid of the returned DEM is the grid used in
    # the model calculation.
    grid = get(opts, :grid, nothing)
    grid!=nothing && error("Changing DEM/simulation grid is not implemented.")

    # Cropbox
    # - to be made automatically from outline (default)
    # - specified in the options
    # - specified by the extent of the DEM
    cropbox = get(opts, :cropbox, :from_outline)
    if cropbox==:from_outline
        bufferfac = get(opts, :bufferfac, 10)
        cropbox, dem = crop(outline, dem, bufferfac)
    elseif cropbox==:from_DEM
        cropbox = (dem.x[1], dem.x[end], dem.y[1], dem.y[end])
    elseif cropbox isa Tuple
        cropbox, dem = crop(cropbox, dem)
    else
        error("cropbox option not recognized: $cropbox")
    end

    # Get/make masks:
    # - make it from outline (and cache it)
    # - read it from file
    masksource = get(opts, :masksource, :from_outline)
    sealevel_cutoff = get(opts, :sealevel_cutoff, 10)
    remove_slivers = get(opts, :remove_slivers, true)
    if masksource==:from_outline
        glaciermask, landmask = get_cache!(ds, outline.gid, [:glaciermask, :landmask],
                                           ()->generate_DEM_masks(outline.gid, dem, outline, sealevel_cutoff, remove_slivers),
                                           Bool, size(dem))
    elseif masksource==:from_file
        glaciermask = convert(Matrix{Bool}, load_asciigrid!(DataSet(ds, files=ds.files[2:2]), :dem_mask, nothing).v)
        landmask = generate_landmask(glaciermask, dem, sealevel_cutoff, remove_slivers)
    else
        error()
    end

    # cutoff: if dem below then mask it there
    cutoff = get(opts, :cutoff, nothing)
    if cutoff!=nothing
        inds = dem.v.<cutoff
        dem.v[inds] = FILL
    end
    # fill the gaps
    mask_unfilled_holes = get(opts, :mask_unfilled_holes, false)
    fillgaps!(dem, glaciermask, mask_unfilled_holes=mask_unfilled_holes)

    # Distance to land:
    # - make it (and cache it)
    # - read it from file
    dist_land_maxdist = opts[:dist_land_maxdist]
    dist_land = get_cache!(ds, outline.gid, :dist_land,
                           () -> dist_to_mask(dem, landmask, glaciermask, dist_land_maxdist),
                           F, size(dem))

    return dem, glaciermask, landmask,
           dist_land, dist_land_maxdist, cropbox,
           dem_proj
end
function load!(ds::DataSet{IVData, :ASCIIGridLoader},
               dem, glaciermask, target_proj, outline, cropbox,
               pp::Phys)::Tuple{SomeData,Matrix{Bool}}
    iv = load_asciigrid!(ds, :iv, cropbox)

    # get/make masks
    ivmask = get_cache!(ds, outline.gid, :ivmask,
                        () -> maskit(iv, outline)[1],
                        Bool, size(iv))

    # fill in any gaps inside the glacier
    fillgaps!(iv, ivmask, mask_unfilled_holes=true)

    return SomeData(vals=iv, sigma=ds.opts[:sigma]), ivmask
end

"""
$(SIGNATURES)

Generate the masks for a glacier.

This can be overloaded for a particular glacier-type.  This assumes:
- If it has outcrops (`hasoutcrops()==true`) then assume that there is
  no land otherwise, just ice and sea.
- If there are no outcrops assume that there is land outside the glaciermask
  except where the ice is below the `sealevel_cutoff`. (ITMIX is more clever,
  maybe use that everywhere)

Returns: glaciermask, landmask
"""
function generate_DEM_masks(gid, dem, outline,
                            sealevel_cutoff=10.0,
                            remove_slivers=true)
    if hasoutcrops(outline)
        #
        glaciermask, landmask = maskit(dem, outline, outline.outcrops)
    else
        glaciermask, _ = maskit(dem, outline)
        landmask = generate_landmask(glaciermask, dem, sealevel_cutoff, remove_slivers)
    end
    return glaciermask, landmask
end

"""
$(SIGNATURES)

Generate the landmask from the glaciermask and a "seamask".  The seamask is generated
assuming anything lower than sealevel_cutoff is sea.

- sealevel_cutoff=10.0
- remove_slivers=true: remove small slivers of land caught between glacier and sea
"""
function generate_landmask(glaciermask, dem, sealevel_cutoff=10.0, remove_slivers=true)
    landmask = (!).(glaciermask)
    # mark points in the sea
    landmask[FILL.<dem.v.<=sealevel_cutoff] = false
    if remove_slivers
        # remove small slivers of land outside of the glacier
        lm = Float64.(landmask)
        lm = boxcar(lm, 1, (!).(glaciermask)) # if more than half the points around a point
        # are different, then change it.  This will
        # lead to some false positives but that should
        # not matter.
        landmask = Bool.(map(x->round(Int,x), lm)) #.6 Bool.(round.(Int,lm))
    end
    return landmask
end
