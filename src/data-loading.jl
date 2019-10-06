using Parameters
# This file sets up the data-loading infrastructure.
export @gid_str, load_glacier

## Note some types are also defined in data-structs.jl

"""
$(TYPEDEF)

Types of different data to be dealt with
"""
abstract type DataKind end
# glacier outline
abstract type OutlineData <: DataKind end
# misc model parameters
abstract type ParaData <: DataKind end
# DEM 2D
abstract type DEMData <: DataKind end
# ice surface velocity
abstract type IVData <: DataKind end
# Surface mass balance
abstract type BdotData <: DataKind end
# dh/dt data
abstract type DhdtData <: DataKind end
# basal sliding fraction
abstract type FslData <: DataKind end
# ice temperature
abstract type TempData <: DataKind end
# ice thickness (e.g. radar or nothing)
abstract type ThicknessData <: DataKind end
# ice flux at terminus
abstract type TerminusFluxData <: DataKind end

"Standard variable names for each DataKind"
datakind2varname = Dict(OutlineData=>:outline,
                        #ParaData=>:para,
                        DEMData=>:dem,
                        IVData=>:iv,
                        BdotData=>:bdot,
                        DhdtData=>:dhdt,
                        FslData=>:fsl,
                        TempData=>:temp,
                        ThicknessData=>:thick)

"""
Returns a dict of default options to individual loaders.  Add opts
here if they are universally used.  For others, modify the dict
returned from this function.
"""
function def_dataset_opts()
    out = Dict{DataType,Any}()
    out[DataKind] = Dict{Symbol,Any}()
    for st in subtypes(DataKind)
        out[st] = Dict{Symbol,Any}()
    end
    out[DEMData][:dist_land_maxdist] = 3e3
    # Set windows to volume-area scaling thickness estimate within below bounds:
    out[ParaData][:pm] = Dict{Symbol,Any}(
        :window_dem_smooth => (0,200),
        :window =>  (0,2000),
        )
    out[ParaData][:pn] = Dict{Symbol,Any}()
    out[ParaData][:pp] = Dict{Symbol,Any}()
    return out
end

function def_dataset_files()
    out = Dict{DataType,Vector{String}}()
    for st in subtypes(DataKind)
        out[st] = [""]
    end
    return out
end
function def_dataset_dateranges()
    out = Dict{DataType,Tuple{Date,Date}}()
    for st in subtypes(DataKind)
        out[st] = (Date(),Date())
    end
    return out
end

# Loading parameters, used to customize the loading process
"""
$(TYPEDEF)

Misc parameters used during loading/synthetic-creation of datasets used for
**customization** of the loading process
Any parameters which are used inside `load_glacier` (with the
exception of some physical parameters) should be in here.

Particular sets of glaciers should overload the constructor of this to
return a suitable pl. Example

    LoadPara(gid::ITMIXGlacier, use_glogem=false; kwargs...) = ...

The flow goes: `LoadPara` instance is used to create the DataSet
instances, which in terms make up the `dt` DataTable instance.  The
`dt` then holds all info needed for data-loading/data-creation:
- Input/output directories and file-paths of data-files.
- The `dataset_opts` present here will be inserted into the `DataSet` instance.
  The options in `dataset_opts[DataKind]` will be added to all via a merge.

    code_root::String=normpath(dirname(dirname(@__FILE__))) # root directory of this package
    data_root::String=normpath(joinpath(code_root, "..", "data")) # root directory of data
    cache_dir::String=joinpath(data_root, "cache") # directory to write cache-files
    update_cache::Bool=false # if true, recalculate all cached files.
    target_proj::String="from DEM" # in what projection to re-project data-sets.  Default "from DEM" uses the projection of the DEM.
    # The following fields mirror the DataSet fields:
    dataset_opts::Dict{DataType,Any}=def_dataset_opts() # holds options for load! functions.
    dataset_LOADERS::Dict{DataType,Symbol}=Dict{DataType,Symbol}() # Dict between the two type parameters for the DataSet instances.
    dataset_files::Dict{DataType,Vector{String}}=def_dataset_files()
    dataset_dateranges::Dict{DataType,Tuple{Date,Date}}=def_dataset_dateranges()
"""
@with_kw struct LoadPara<:APara @deftype String
    gid::GlacierID
    code_root::String=normpath(dirname(@__DIR__)) # root directory of this package
    data_root::String=normpath(joinpath(code_root, "scripts", "data")) # root directory of data
    cache_dir::String=joinpath(data_root, "cache") # directory to write cache-files
    update_cache::Bool=false # if true, recalculate all cached files.
    target_proj::String="from DEM" # in what projection to re-project data-sets.  Default "from DEM" uses the projection of the DEM.
    # The following fields mirror the DataSet fields:
    dataset_opts::Dict{DataType,Any}=def_dataset_opts() # holds options for load! functions.
    dataset_LOADERS::Dict{DataType,Symbol}=Dict{DataType,Symbol}() # Dict between the two type parameters for the DataSet instances.
    dataset_files::Dict{DataType,Vector{String}}=def_dataset_files()
    dataset_dateranges::Dict{DataType,Tuple{Date,Date}}=def_dataset_dateranges()
end

"""
$(TYPEDEF)
$(FIELDS)

Holds information needed to load or create a data-set, which can be a
single data-set or a set of related one (e.g. a DEM and its mask.
Although there the outline comes in too, which is separate...).

Parameterized on the DataKind it holds and the LOADER
it needs to be loaded.

Common opts:
- proj -- expected projection
"""
@with_kw struct DataSet{DK<:DataKind, LOADER}
    "File names, if applicable."
    files::Vector{String}
    "Extra options for loader functions.  Unknown options are (usually) ignored."
    opts::Dict{Symbol,Any}
    "Date range of measurements.  (Date(),Date()) if unknown (yet)."
    daterange::Tuple{Date,Date}
    "File-cache directory"
    datacache_dir::String # File-cache directory.  Note that there may
                          # be several cache-files associated with a
                          # DataSet. Say several masks for a DEM containing
                          # several glaciers.
    "If true, recalculate the cache files"
    datacache_dir_update::Bool
    "If true, do not use the RAM cache"
    datacache_ram_update::Bool
    function DataSet{DK,L}(files, opts, daterange=(Date(),Date()),
                           datacache_dir="", datacache_dir_update=false,
                           datacache_ram_update=false) where {DK,L}
        new{DK,L}(files, opts, daterange, datacache_dir, datacache_dir_update,
                  datacache_ram_update)
    end
end
Base.convert(::Type{DataSet{DK,L}}, ds::DataSet) where {DK,L} =
    DataSet{DK,L}(ds.files, ds.opts, ds.daterange, ds.datacache_dir, ds.datacache_dir_update, ds.datacache_ram_update)
Base.show{DK,LS}(io::IO, ds::DataSet{DK,LS}) =
    println(io, "DataSet{$DK,$LS}: $(ds.files), $(ds.daterange)")

# for caching to RAM use hash:
function Base.hash(ds::DataSet{DK,L}, h::UInt) where {DK,L}
    h += hash(DK,h)
    h += hash(L,h)
    h += hash(ds.files,h)
    h += hash(ds.opts,h)
    h += hash(ds.daterange,h)
    h += hash(ds.datacache_dir,h)
    h += hash(ds.datacache_dir_update,h)
    return h
end


function Base.:(==){DK1,LS1,DK2,LS2}(ds1::DataSet{DK1,LS1}, ds2::DataSet{DK2,LS2})
    DK1==DK2 || return false
    LS1==LS2 || return false
    for fl in fieldnames(ds1)
        getfield(ds1,fl)==getfield(ds2,fl) || return false
    end
    return true
end

"""
A table which holds where/how to load/produce the data for a
particular glacier and a particular run.  (i.e. a glacier may have
several datasets)

Fields:

    gid::GlacierID
    outline
    para
    dem
    iv
    bdot
    dhdt
    fsl
    temp
    thickness

"""
@with_kw_noshow struct DataTable
    gid::GlacierID
    outline::DataSet{OutlineData}
    para::DataSet{ParaData}
    dem::DataSet{DEMData}
    iv::DataSet{IVData}
    bdot::DataSet{BdotData}
    dhdt::DataSet{DhdtData}
    fsl::DataSet{FslData}
    temp::DataSet{TempData}
    thickness::DataSet{ThicknessData}
    terminus_flux::DataSet{TerminusFluxData}
end
function Base.show(io::IO, dt::DataTable)
    println(io, "Data Table for Glacier: $(dt.gid)")
    for fl in fieldnames(dt)
        show(io, getfield(dt, fl))
    end
    println("")
end
"Return all known datakinds in the order used in DataTable"
datakinds() = [_get_first_para(ST)for ST in DataTable.types[2:end]]
_get_first_para(::Type{DataSet{T}}) where {T} = T

"""
$(SIGNATURES)
$(METHODLIST)

Constructs a `DataTable` for a glacier with `GlacierID` and `pl`
loading options.  It will choose reasonable defaults for a glacier for
all the data inputs.  The output, a `DataTable` can later be updated to
choose different data sources for some fields.
"""
function make_datatable(gid::GlacierID, pl::LoadPara, dt=nothing)
    @assert gid==pl.gid "Not matching gid: $gid, $(pl.gid)"
    @unpack dataset_LOADERS, dataset_files,
            dataset_opts, dataset_dateranges = pl
    args = []
    for DK in datakinds()
        push!(args, DataSet{DK,dataset_LOADERS[DK]}(
            dataset_files[DK],
            merge(dataset_opts[DataKind],dataset_opts[DK]),
            dataset_dateranges[DK],
            pl.cache_dir,
            pl.update_cache,
            pl.update_cache
        ))
    end
    return DataTable(gid, args...), pl
end

# """
# Hold data for a model run.
# """
# struct ModelRunData
#     datatable::DataTable
#     outfile::String # path to a JLD file
#     comments::Dict
# end
# function ModelRunData(gid, outfile, comments=Dict())
#     datatable = getdatatable(gid)
#     ModelRunData(datatable, outfile, comments)
# end

#########
# Loader functions
#########
# of form

"""
    load!(ds::DataSet{OutlineData, LS}, gid::GlacierID, target_proj)
    load!(ds::DataSet{DEMData, LS}, gid::GlacierID, target_proj, outline, pp::Phys)
    load!(ds::DataSet{OtherData, LS}, gid::GlacierID, target_proj, outline, cropbox,
          pp::Phys)

Data loader/synthesizing functions.  Overload `load!` functions for
needed combinations of `DataKind` and `LOADER` (some symbol, often
only one combo makes sense).  The `ds.opts` hold options for specific
loaders.

The `load!` function modifies `ds` in that it stores the raw(ish)
loaded data into `ds.datacache`.  This allows caching data like DEMs
which maybe applicable to several glaciers.  If the `ds.datacache` is
already filled then the reading from disc step is skipped.

The primary variable it returns should be the loaded data tailored to
the glacier at hand.
"""
function load! end

"""
    load!(::DataSet{OutlineData, :dummy}, gid::GlacierID)

Loading the outline determines which glacier to run.

Return
- outline: the outline of the glacier including land-islands.
- outcrops: nunataks, or similar outcrops which are not included in
  outline.  Their insides are all treated as land for masking
  purposes.  (This is a bit redundant with the outline, but some data
  comes in this format.)  Return (
- iscalving::Bool -- true for calving glacier (tide-water, fresh-water, hanging glacier)
- proj -- projection as Proj4 string (or something else)
"""
function load!(::DataSet{OutlineData, :dummy}, gid::GlacierID)
    if target_proj=="from DEM"
        proj = use_Outline_proj() # if needed, it get re-projected later
    end

    return outline::Outline, iscalving::Bool, proj::String
end

"""
    load!(::DataSet{ParaData}, gid::GlacierID, outline::Outline, isicecap::Bool)

Load parameters for forward model.  Probably this generic
implementation is fine.  This default loader returns the default
para-structs which are modifed by `ds.opts[:pp]`, etc.  Can
additionally be overloaded for `gid`.

The `outline` and `isicecap` input are used to automatically set the
smoothing windows to approximate mean ice thickness (estimated with
volume area scaling).

Return
- pp::Phys
- pm::MPara
- pn::Num
TODO: add inverse parameters too.
"""
function load!(ds::DataSet{ParaData}, gid::GlacierID, outline::Outline, isicecap::Bool)
    opts = ds.opts
    if isa(opts[:pm][:window],Tuple) || isa(opts[:pm][:window_dem_smooth],Tuple)
        # set bed-DEM smoothing window with volume-area
        hmean = volume_area_mean_h(area(outline), isicecap=isicecap)
        if isa(opts[:pm][:window],Tuple)
            mi,ma = opts[:pm][:window]
            opts[:pm][:window] = min(ma, max(mi,1.5*hmean) )
        end
        if isa(opts[:pm][:window_dem_smooth],Tuple)
            mi,ma = opts[:pm][:window_dem_smooth]
            opts[:pm][:window_dem_smooth] = min(ma, max(mi,hmean) )
        end
    end
    pp = Phys(; opts[:pp]...)
    pm = MPara(; opts[:pm]...)
    pn = Num(; opts[:pn]...)
    return pp::Phys, pm::MPara, pn::Num
end

"""
    load!(ds::DataSet{DEMData, :dummy}, outline,
               pp::Phys)

ds.opts:
- grid : if specified use a different grid from the grid of the DEM as
  computational grid.  Given in the `target_proj`.
- cropbox -- used for cropping the DEM and other 2D data.  Needs to
  contain all of the outline. `F[xmin, xmax, ymin, ymax]`, defaults
  being calculated during DEM-load with some space around the outline
  (needed for interpolations).

Return:
- dem::Gridded
- glaciermask::Matrix{Bool} -- true on glacier
- landmask::Matrix{Bool} -- true on solid ground (i.e. false on glaciers, ice shelves, water)
- dist_land::Matrix{F} -- distance to closest land
- dist_land_maxdist -- maximum distance to which to calculate above
- cropbox::Vector{Float64} -- box around region of interest (contains glacier + some buffer)
                              `F[xmin, xmax, ymin, ymax]`
"""
function load!(::DataSet{DEMData, :dummy}, outline,
               pp::Phys)
    return dem::Gridded, glaciermask::Matrix{Bool}, landmask::Matrix{Bool},
           dist_land::Matrix{F}, dist_land_maxdist::F, cropbox::Vector{Float64},
           dem_proj::String
end
function load!(::DataSet{IVData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    return iv::SomeData, ivmask::Matrix{Bool}
end
function load!(::DataSet{ThicknessData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    return thickness::SomeData
end
function load!(::DataSet{BdotData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    return bdot::SomeData
end
function load!(::DataSet{DhdtData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              iscalving, pp::Phys)
    return dhdt::SomeData
end
function load!(::DataSet{FslData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    return fsl::SomeData
end
function load!(::DataSet{TempData, :dummy}, dem, glaciermask, target_proj, outline, cropbox,
              pp::Phys)
    return temp::SomeData
end
# scalar loaders
"""
    load!(::DataSet{TerminusFluxData}, iscalving, pp::Phys)

Ice flux at terminus.
"""
function load!(::DataSet{TerminusFluxData}, iscalving, pp::Phys)
    error()
end


"""
     special_processing!(gid::GlacierID, outline, dem, glaciermask, landmask,
                             iv, ivmask, thickness, bdot, dhdt, fsl, temp,
                             pl::LoadPara, pp::Phys, pm::MPara, pn::Num)

Allows to do some special processing on a glacier by glacier basis by
overloading this function.
"""
function special_processing!(gid::GlacierID, outline, dem, glaciermask, landmask,
                             iv, ivmask, thickness, bdot, dhdt, fsl, temp,
                             pl::LoadPara, pp::Phys, pm::MPara, pn::Num)
    nothing
end


"""
$(SIGNATURES)

Load/synthesize all data for a model run of a glacier using either:
- a `DataTable`
- a `DataTable` and an already "loaded" `DataTable`
- a glacier ID (TODO: implement)
as input.  For the latter the default dataset is used.  (It
potentially modifies the `_datacache_ram` field of the `DataSet`s, but
this is an internal fields, thus no !."

Input
- glacier-id or datatable
- load-para

Return
- Glacier
- pp, pm, pn
"""
function load_glacier end

function load_glacier(gid::GlacierID, pl::LoadPara)
    dt, pl = make_datatable(gid, pl)
    load_glacier(dt,pl)
end

function load_glacier(datatable::DataTable, pl::LoadPara)
    @assert datatable.gid==pl.gid  "Not matching gid: $gid, $(pl.gid)"
    @unpack gid, outline, para, dem, iv, bdot, dhdt,
            fsl, temp, thickness, terminus_flux = datatable

    @unpack target_proj = pl # TODO do cache dir

    # TODO: maybe reverse loading order: first load DEM with grid, then
    # outline.  With outline come the cropbox & masks, which would be more
    # consistent.  Make DEM & IV masks their separate loader?  dist_land
    # is also generated with the masks.

    outline, iscalving, isicecap, outline_proj = load!(outline, gid)

    # parameters, updated using the optional input parameters
    pp, pm, pn = load!(para, gid, outline, isicecap)

    # note some variables get re-bound:
    (dem, glaciermask, landmask, dist_land,
     dist_land_maxdist, cropbox, dem_proj) = load!(dem, outline, pp)

    @assert !any(dem.v[glaciermask].==FILL) "The gl.dem has FILL values inside the glacier mask, probably mask those points"
    # projections
    target_proj = target_proj=="from DEM" ? dem_proj : target_proj
    @assert dem_proj==outline_proj==target_proj "Need to have same projection"
    # dem, glaciermask, landmask, dist_land, cropbox, target_proj =
    #     reproject_DEM(dem, glaciermask, landmask, dist_land, cropbox, dem_proj, target_proj)
    # outline, outcrops = reproject_outline(outline, outcrops, outline_proj, target_proj)

    ## TODO update about projections
    iv, ivmask = load!(iv, dem, glaciermask, target_proj, outline, cropbox, pp)
    thickness = load!(thickness, dem, glaciermask, target_proj, outline, cropbox, pp)
    terminus_flux = load!(terminus_flux, dem, glaciermask, iscalving, pp)
    bdot = load!(bdot, dem, glaciermask, target_proj, outline, cropbox, pp)
    dhdt = load!(dhdt, bdot, terminus_flux, dem, glaciermask, target_proj,
                 outline, cropbox, iscalving, pp)
    fsl = load!(fsl, dem, glaciermask, target_proj, outline, cropbox, pp)
    temp = load!(temp, dem, glaciermask, target_proj, outline, cropbox, pp)

    special_processing!(gid, outline, dem, glaciermask, landmask,
                        iv, ivmask, thickness, bdot, dhdt, fsl, temp,
                        pl, pp, pm, pn)
    alpha, elemean, elemax, elemin, elemedian = proc_dem(dem, glaciermask)

    return Glacier(gid,
                   iscalving, isicecap,
                   target_proj,
                   dem,
                   outline, cropbox,
                   glaciermask, landmask,
                   dist_land, dist_land_maxdist,
                   alpha, elemean, elemax, elemin, elemedian,
                   bdot, dhdt,
                   temp, fsl,
                   terminus_flux,
                   thickness, iv, ivmask,
                   Dict{Symbol,Any}()),
           pp,pm,pn,pl
end

# function load_glacier(datatable::DataTable, datatable_other::DataTable, pl::LoadPara)
#     @unpack gid, outline, para, dem, iv, bdot, dhdt, fsl, temp, thickness = datatable
#     ot = datatable_other
#     for fl in fieldnames(datatable)
#         fl==:gid && continue
#         if compare_nocache(getfield(datatable, fl), getfield(datatable_other, fl))
#             setfield!(datatable, fl) = getfield(datatable_other, fl)
#         end
#     end
#     load_glacier(datatable, pl)
# end

################
# All together
################

"""
    init_forward(gid::GlacierID, verbose=true; kwargs_for_LoadPara...)

Initialize forward model.  Returns

-  gl,gb,pp,pm,pn,pl
"""
function init_forward(gid::GlacierID, verbose=true; kwargs_for_LoadPara...)
    verbose && print("Loading glacier $gid ... ")
    pl = LoadPara(gid; kwargs_for_LoadPara...)
    t = @elapsed gl,pp,pm,pn,pl = load_glacier(gid,pl)
    if verbose
        println("$(round(t,3)) sec")
        print("Making bands ... ")
    end
    pn.test && Base.Test.@inferred make_bands(gl, pp, pm, pn)
    t = @elapsed gb,pm = make_bands(gl, pp, pm, pn)
    verbose && println("$(round(t,3)) sec")
    return gl,gb,pp,pm,pn,pl
end

#################
# Helper functions see data-loading-misc.jl
################
