# Various processing functions needed in data-loading.jl

"Fill value to avoid using NaNs for performance (TODO: is this necessary?)"
const FILL = F(-9999999)

"""
    proc_dem(dem::Gridded, mask::Matrix{Bool})

Process elevation DEM.

out:     alpha, elemean, elemax, elemin, elemedian
"""
function proc_dem(dem::Gridded, mask::Matrix{Bool})
    m = dem.v[mask]
    return VAWTools.absslope(dem,mask), mean(m), maximum(m), minimum(m), median(m)
end


"""
    maskit(dem::Gridded, outline, outcrops) # assumes only outcrops are solid ground
    maskit(dem::Gridded, outline) # assumes anything outside outline is solid ground
    maskit(dem::Gridded, xy) # assumes anything outside outline is solid ground

Makes the glaciermask and landmask.

The landmask can be refined with generate_landmask.

size(outline) ==(2,n)

Output:
- glaciermask==true on the glacier of interest
- landmask==true on all solid land

TODO: allow taking the common boundary of the RGI setup into account.
"""
function maskit(dem::Gridded, outline::Outline, outcrops)
    glaciermask = Array{Bool}(size(dem)) # true on glacier
    mask_outline = Array{Bool}(size(dem)) # true if inside outline
    mask_outcrops = Array{Bool}(size(dem)) # true if NOT on solid ground
    for (j,y)=enumerate(dem.y), (i,x)=enumerate(dem.x)
        mask_outline[i,j] = inpoly((x,y), outline.outline)
        mask_outcrops[i,j] = !inpoly((x,y), outcrops)
        glaciermask[i,j] = mask_outline[i,j] && mask_outcrops[i,j]
    end
    groom_glaciermask!(glaciermask)
    glaciermask, .!mask_outcrops
end
maskit(dem::Gridded, outline::Outline) = maskit(dem, outline.outline)
function maskit(dem::Gridded, xy)
    glaciermask = Array{Bool}(size(dem)) # true on glacier
    for (j,y)=enumerate(dem.y), (i,x)=enumerate(dem.x)
        glaciermask[i,j] = inpoly((x,y), xy)
    end
    groom_glaciermask!(glaciermask)
    glaciermask, (!).(glaciermask)
end

"""
    groom_glaciermask!(glaciermask)

Removes any points from the glaciermask which are not connected to any
others at at least one of its edges.

Note, at the same time the landmask should be groomed as well but this
would only have a very minor impact on dist_land.
"""
function groom_glaciermask!(glaciermask)
    window = 1
    ni,nj=size(glaciermask)
    for j=1:nj, i=1:ni
        !glaciermask[i,j] && continue
        # check points around:
        glaciermask[max(1,i-1),j] && continue
        glaciermask[min(ni,i+1),j] && continue
        glaciermask[i,max(1,j-1)] && continue
        glaciermask[i,min(nj,j+1)] && continue
        # if we get here then there are no neighbors:
        glaciermask[i,j] = false
    end
    nothing
end

"""
    fillgaps!(g::Gridded, mask, di=1; mask_unfilled_holes=false)

Fills in FILL holes into 1D and 2D gridded data (dem and bdot).

The FILL holes are filled by interpolating around it.  The grid-distance for
interpolation is given by `di`

If mask_unfilled_holes is set, then the mask will be set to false where
there are unfilled holes, otherwise it errors
"""
function fillgaps!(g::Gridded, mask, di=1; mask_unfilled_holes=false)
    sz1,sz2 = size(g.v)
    @inbounds for (j,y)=enumerate(g.y), (i,x)=enumerate(g.x)
        # fill in FILL-gaps inside the glacier:
        if mask[i,j] && g.v[i,j]==FILL
            val = 0.0
            n = 0
            for jj=max(1,j-di):min(sz2,j+di)
                for ii=max(1,i-di):min(sz1,i+di)
                    tmp = g.v[ii,jj]
                    if tmp!=FILL
                        val += tmp
                        n+=1
                    end
                    # ii,jj = VAWTools.around(g, i, j, 1.5*step(g.x))
                    # g.v[i,j] = VAWTools.meanfill(g.v[ii,jj], FILL)
                end
            end
            if n>0
                g.v[i,j] = val/n
            else
                if mask_unfilled_holes
                    mask[i,j] = false
                else
                    error("Unfilled gap present in mask!")
                end
            end
        end
    end
    nothing
end

function fillgaps!(g::Gridded1d, di = 1)
    len = length(g)
    @inbounds for (i,x) in enumerate(g.x)
        val = 0.0
        n = 0
        for ii=max(1,i-di):min(len,i+di)
            tmp = g.v[ii]
            if tmp!=FILL
                val += tmp
                n+=1
            end
        end
        if n>0
            g.v[i] = val/n
        else
            error("Unfilled gap present!")
        end
    end
    nothing
end


"""
    make_cropbox(outline, buffer)

Return a crop-box [xmin, xmax, ymin, ymax]
"""
function make_cropbox(outline, buffer)
    outline = outline.outline
    xymax = maximum(outline,2)+buffer
    xymin = minimum(outline,2)-buffer
    return xymin[1], xymax[1], xymin[2], xymax[2]
end


"""
    crop{T}(cropbox::Tuple, tocrop)
    crop{T}(outline::Outline, tocrop, bufferfac)

Crop 2D data to cropbox or outline

In: outline::Outline/cropbox::Tuple, Vector{Gridded} (, bufferfac)

Out:
cropbox, Gridded, Gridded...

TODO: use Arrayviews
"""
function crop(cropbox::Tuple, tocrop::Vector{<:Gridded})
    tc = tocrop[1]
    for tcc in tocrop
        if tc.x!=tcc.x && tc.y!=tcc.y
            error("All inputs need to be on same grid")
        end
    end
    xmin, xmax, ymin, ymax = cropbox
    ix = searchsortedfirst(tc.x, xmin):searchsortedlast(tc.x, xmax)
    iy = searchsortedfirst(tc.y, ymin):searchsortedlast(tc.y, ymax)
    xs = tc.x[ix]
    ys = tc.y[iy]
    return cropbox, [Gridded(xs,ys,tc.v[ix,iy]) for tc in tocrop]...
end
crop(cropbox::Tuple, tocrop::Gridded) = crop(cropbox, [tocrop])
function crop(outline::Outline, tocrop::Vector{<:Gridded}, bufferfac)
    cropbox = make_cropbox(outline, step(tocrop[1].x)*bufferfac)
    return crop(cropbox, tocrop)
end
crop(outline::Outline, tocrop::Gridded, bufferfac) = crop(outline, [tocrop], bufferfac)

function crop(cropbox::Tuple, tocrop::Traj)
    xmin, xmax, ymin, ymax = cropbox
    inds = Int[]
    splits = Int[1]
    for s in tocrop.splits
        n = 0
        for i in s
            x,y = tocrop.x[i], tocrop.y[i]
            if xmin<=x<=xmax && ymin<=y<=ymax
                push!(inds,i)
                n+=1
            end
        end
        if n>0
            push!(splits, splits[end]+n)
        end
    end
    if length(inds)==1
        splits = [1:0]
    else
        splits = [splits[i]:splits[i+1]-1 for i=1:length(splits)-1]
    end

    err = VAWTools.haserror(tocrop) ? tocrop.err[inds] : tocrop.err
    v = VAWTools.hasvalues(tocrop) ? tocrop.v[inds] : tocrop.v
    Traj(tocrop.x[inds], tocrop.y[inds], v, err, splits, tocrop.proj)
end
crop(cropbox::Tuple, tocrop::Vector{<:Traj}) = map(t->crop(cropbox,t), tocrop)


"""
    dist_to_mask(dem::Gridded, mask4dist::AbstractMatrix, mask4points::AbstractMatrix, dist_max)

    # usual


Calculate approximate distance to the edge of the closest masked
point, i.e. where mask4dist==true.  Only consider points where
mask4points==true.  Don't search for more than dist_max.

TODO: this could be done more cleverly using a spiraling search pattern.
"""
function dist_to_mask(dem::Gridded, mask4dist::AbstractMatrix, mask4points::AbstractMatrix, dist_max)
    @assert size(mask4points)==size(mask4dist)
    sz1,sz2 = size(mask4dist)
    dist = similar(dem.v)
    dx = step(dem.x)
    di = ceil(Int,dist_max/dx)
    # iterate over all points
    @inbounds @fastmath for j=1:sz2, i=1:sz1
        if mask4points[i,j] # only check points inside glacier
            mindist = dist_max
            # Iterate over box of size 2dist_max x 2dist_max
            for jj=max(1,j-di):min(sz2,j+di)
                for ii=max(1,i-di):min(sz1,i+di)
                    if mask4dist[ii,jj]
                        mindist = min(mindist, distance(i,j,ii,jj,dx))
                    end
                end
            end
            dist[i,j] = max(0,mindist) # in case a land-point sneaks into the glaciermask
        else
            dist[i,j] = FILL
        end
    end
    return dist
end
# The `- dx/2` is a rough way to account for the distance to the edge
# and not the cell center.
@fastmath distance(i,j,ii,jj,dx) = (dx * sqrt((i-ii)^2+(j-jj)^2) - dx/2)

####
# Cache in RAM or Disk
####
"""
     get_cache!{DK,LS}(ds::DataSet{DK,LS}, gid, datanames::Vector{Symbol},
                           generate_fn, T::Array, sz;
                           cache2RAM=false, # as dataset is glacier-specific
                           caching2file_treshold=0.4,  # rates then cache.
                           )
     # Single array
     get_cache!{DK,LS}(ds::DataSet{DK,LS}, gid, dataname::Symbol,
                           generate_fn, T::Array, sz;
                           cache2RAM=false,
                           caching2file_treshold=0.4,  # rates then cache.
                           )


Fetch values from RAM-cache or disk-cache; or produce them with
`generate_fn` if there is no cached values.  It can update the
RAM-cache, thus the `!`.

This is *specific* to combo of a glacier and a DataSet.  Typically
this will be used with a function making an expensive product of a
data-set, say a glaciermask.

- `generate_fn` has signature `()->...` and returns a Vector of
  outputs.  All elements need to be of the same `Array` type.

Notes:
- only works for `Array`, all need to have the same size
- no type checking done for RAM cache

TODO: use https://github.com/ExpandingMan/Anamnesis.jl instead or
      https://github.com/simonster/Memoize.jl
"""
function get_cache! end

function get_cache!(ds::DataSet{DK,LS}, gid, datanames::Vector{Symbol},
                           generate_fn, T, sz;
                           cache2RAM=false, # as dataset is glacier-specific
                           caching2file_treshold=0.4,  # rates then cache.
                           verbose=false
                    ) where {DK,LS}
    @assert length(datanames)>0
    force = ds.datacache_dir_update
    cache_dir = ds.datacache_dir
    verbose && print("Checking File-cache for $DK, $LS, $datanames: ")

    # in RAM, use it:
    if !force && mapreduce(n->has_cache_ram(ds,n,sz,gid), &, datanames)
        out = Array{T,length(sz)}[]
        verbose && println("using ram")
        for dn in datanames
            push!(out, _CACHE_RAM[cache_ram_key(ds,dn,sz,gid)])
        end
        return out
    end
    # If in file-cache, use it:

    # To make a unique filenames: hash the size and datatype
    hsh = hash(hash(T), hash(sz))
    filenames = [joinpath(cache_dir,
                          "$LS-$(DK.name.name)-$(getid(gid))-$dataname-$hsh.bin") for dataname in datanames]

    if all(isfile.(filenames)) && !force
        out = Array{T,length(sz)}[read(filename, T, sz) for filename in filenames]
        verbose && println("using files")
        # also add to RAM-cache
        if cache2RAM
            for (i,dn) in enumerate(datanames)
                _CACHE_RAM[cache_ram_key(ds,dn,sz,gid)] = out[i]
            end
        end
        return out
    end
    # Else, generate it and cache as appropriate.
    verbose && println("generating")

    out = Array{T,length(sz)}[]
            time = tic()
    append!(out, generate_fn())
    delta_t = tic()-time
    @assert length(out)==length(datanames)
    delta_mem = mapreduce(sizeof, +, 0, out)
    # cache to RAM
    if cache2RAM
        for (i,dn) in enumerate(datanames)
            _CACHE_RAM[cache_ram_key(ds,dn,sz,gid)] = out[i]
        end
    end
    if delta_mem/delta_t<caching2file_treshold && cache_dir!=""
        if !isdir(cache_dir)
            mkdir(cache_dir)
        end
        for (f,o) in zip(filenames,out)
            write(f, o)
        end
    end
    return out
end
function get_cache!(ds::DataSet{DK,LS}, gid, dataname::Symbol,
                           generate_fn, T, sz;
                           cache2RAM=false, # as dataset is glacier-specific
                           caching2file_treshold=0.4,  # rates then cache.
                           ) where {DK,LS}
    fn = () -> [generate_fn()]
    get_cache!(ds::DataSet, gid, [dataname],
               fn, T, sz;
               cache2RAM=cache2RAM,
               caching2file_treshold=caching2file_treshold
                    )[1]
end

"To cache expensive to load/make data in RAM, so it can be used for several runs."
const _CACHE_RAM = Dict{UInt,Any}()

cache_ram_clear!() = empty!(_CACHE_RAM)

cache_ram_key(ds::DataSet, field::Symbol, size_, gid) = hash(ds) + hash(field) + hash(size_) + hash(gid)
has_cache_ram(ds::DataSet, field::Symbol, size_, gid) = haskey(_CACHE_RAM, cache_ram_key(ds,field,size_,gid))

"""
$(SIGNATURES)

Looks in the dictionary `_CACHE_RAM` to see whether the key `field` exists.
- If yes, return copy of the value of that field, unless ds.datacache_ram_update=true
- If no, generate value with `generate_fn`, store it in `cache_ram` and return copy.

This is not specific to a glacier but only to a DataSet.  Typically
this will be used with a file-loading function.

Note: this supports any datatype, unlike `get_cache!` which can only do arrays.
"""
function get_cache_ram!(ds::DataSet{D,L}, field::Symbol, generate_fn,
                        size_=(), verbose=false,
                        cache2RAM=false) where {D,L}
    force = ds.datacache_ram_update
    hs = cache_ram_key(ds,field,size_,nothing)
    verbose && print("Accessing _CACHE_RAM for $D, $L, $field, $hs")
    if !force && haskey(_CACHE_RAM, hs)
        # get cache
        verbose && println(" ... geting cache.")
        out = deepcopy(_CACHE_RAM[hs])
    else
        # generate it
        out = generate_fn()
        verbose && println(" ... loading/generating.")
        # and cache it
        if cache2RAM
            _CACHE_RAM[hs] = deepcopy(out)
        end
    end
    return out
end

########
## Projections
########

function reproject_outline(outline, outcrops, outline_proj, target_proj)
    @assert outline_proj==target_proj "Re-projecting outlines not implemented yet"
    # check VAWTools.jl
    return outline, outcrops
end

function reproject_DEM(dem, glaciermask, landmask, dist_land, cropbox, dem_proj, target_proj; grid=nothing)
    if target_proj=="from DEM" || target_proj==dem_proj
        target_proj = dem_proj
        return dem, glaciermask, landmask, dist_land, cropbox, target_proj
    else
        @assert grid!=nothing "Need to specify grid, if DEM is reprojected"
        error("Not implemented yet")
        # Interpolate elevation at `grid` points: project each grid
        # point into DEM projection & interpolate elevation.

        return dem, glaciermask, landmask, dist_land, cropbox, target_proj
    end
end

################

"Scale from m w.e./a to m ice/a"
scale_m_ice!(out::Gridded, pp::Phys, extra=1.0) = scale_m_ice!(out.v,pp,extra)
function scale_m_ice!(out::Array, pp::Phys, extra=1.0)
    # m w.e./a -> m ice/a
    @inbounds @fastmath for i in eachindex(out)
        if out[i]!=FILL
            out[i] = out[i] *pp.rhow/pp.rho*extra
        end
    end
    nothing
end

######################3
using ZipFile, CodecZlib

# Opens a file, either normal or within a zip-file
# Read from the first, close the second return
function open_file(fl::String, orig=fl)
    if fl=="/" || fl=="/.zip"
        error("File $orig not readable, nor is it part of a zip-file")
    elseif isfile(fl)
        if splitext(fl)[2]==".zip"
            actual_file = splitdir(orig)[2]
            zfl = ZipFile.Reader(fl)
            for f in zfl.files
                if splitdir(f.name)[2]==actual_file
                    return f, zfl
                end
            end
            close(zfl)
            error("File $actual_file not found in zip file $(fl)")
        elseif endswith(fl, ".gz")
            fid = CodecZlib.GzipDecompressorStream(open(fl))
            return fid, fid
        else
            fid = open(fl, "r")
            return fid, fid
        end
    else
        return open_file(splitdir(fl)[1]*".zip", orig)
    end
    error("something is wrong with this method")
end


"""
    @open_magic some_reader_function(... , fl, ...)

Magically, opens normal files or zip-files & reads them with the
`some_reader_function`.  Returns the output of the
`some_reader_function`.

Note that "fl" needs to be used where the file slots in.
"""
macro open_magic(fn)
    esc(quote
        fl, fid = open_file(ds.files[1])
        out = $fn
        close(fid)
        out
    end)
end
"""
Usage

    @open_magic_fn some_reader_function(... , fl, ...)

As `@open_magic` but creates a function of zero arguments to execute
above.

Note that "fl" needs to be used where the file slots in.
"""
macro open_magic_fn(fn)
    esc(quote
        function ()
          @open_magic $fn
        end
    end)
end
