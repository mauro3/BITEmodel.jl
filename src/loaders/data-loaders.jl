# All implementations of data loaders/creaters are in this file or loaded from here

export SyntheticGlacier, ITMIXGlacier, RGIGlacier

# generic loader functions
include("generic.jl")

"""
Synthetic glaciers:
- Bench
"""
struct SyntheticGlacier <: GlacierID
    id::Symbol
end
getname(gid::SyntheticGlacier) = string(gid.id)
getregion(::SyntheticGlacier) = Region.Synthetic

const SyntheticBenchLoader = :SynthBench
include("synthetic-bench.jl")

"""
If there is nothing to load, return `SomeData(nothing)`.

Probably only makes sense for iv and thickness, thus only those are implemented.
"""
const NothingLoader = :Nothing
load!(::DataSet{IVData, NothingLoader}, args...) = SomeData(nothing), Matrix{Bool}(0,0)
load!(::DataSet{ThicknessData, NothingLoader}, args...) = SomeData(nothing)


"All ITMIX glaciers ever defined.  Check to code for which versions"
const ITMIXnames = [:Academy, :Aqqutikitsoq, :Austfonna, :AustreGroenfjordbreen, # v2
                    :Brewster, :ChhotaShigri, # v2
                    :Columbia, :Devon, :Elbrus, :Freya, :Hellstugubreen, :Kesselwandferner,
                    :Mocho, :NorthGlacier, :SouthGlacier, :Starbuck, :Synthetic1, :Synthetic2,
                    :Synthetic3, :Tasman, :Unteraar, :Urumqi, :Washmawapta]
const itmix2only = [:AustreGroenfjordbreen, :ChhotaShigri]

"""
ITMIX glacier data
"""
struct ITMIXGlacier <: GlacierID
    id::Symbol # the name as symbol
    version::Int # ITMIX version
    function ITMIXGlacier(name, version)
        @assert 0<version<3
        @assert name in ITMIXnames "$name not in BITEModel.ITMIXnames"
        new(name, version)
    end
end

getname(gid::ITMIXGlacier) = string(gid.id)
getregion(gid::ITMIXGlacier) = Region._Region(gid.version+200)
function getrgi(gid::ITMIXGlacier)
    error("Function not up to date")
    r = iscalving_isicecap[gid.id][3]
    r==nothing && error("Glacier $gid has no RGI number.")
    return r
end
const ITMIXLoader = :ITMIX

ITMIXv1_glaciers = [ITMIXGlacier(n,1) for n in ITMIXnames if !(n in itmix2only)]
ITMIXv2_glaciers = [ITMIXGlacier(n,2) for n in ITMIXnames]
include("itmix.jl")

"""
As defined in the RGI 6.0 technical document.
"""
@with_kw struct RGIattrs
    # bgndate::Date
    # enddate::Date
    lon::F
    lat::F
    subregion::Int
    area::F
    zmin::F
    zmax::F
    # zmed::F
    # slope::F
    aspect::F
    # lmax::F
    status::Int # 0,1 ok, 2 not proper outline
    form::Int
    termtype::Int
    surging::Int
    name::String
end
function RGIattrs(id::Union{String,Symbol})
    reg = getrgi_region(string(id))
    nr = getrgi_nr(string(id))
    RGIattrs(read_attrib_file(reg, nr)[2:end]...)
end

"""
Randolph glacier inventory (RGI) glaciers

http://www.glims.org/RGI/index.html

Use RGI 6.0 numbers, e.g. RGI60-07.01455

Note that RGI itself only provides the outlines.  However Matthias
also prepared DEM files.
"""
struct RGIGlacier <: GlacierID
    id::Symbol
    attrs::RGIattrs
    function RGIGlacier(id::Union{String,Symbol}, attrs::RGIattrs)
        @assert getrgi_version(string(id))==v"6.0" "Need RGI version 6.0 IDs"
        new(uppercase(string(id)), attrs)
    end
    RGIGlacier(id::Union{String,Symbol}) = RGIGlacier(id, RGIattrs(id))
end
function RGIGlacier(nr::Int, region::Int)
    r = VAWTools.int2str2(region)
    id = VAWTools.int2str5(nr)
    RGIGlacier("RGI60-$r.$id")
end
# getid uses default
getrgi(gid::RGIGlacier) = string(getid(gid)) # same as ID in this case, but string
getregion(gid::RGIGlacier) = Region._Region(getrgi_region(gid))
"RGILoader"
const RGILoader = :RGI
getname(gid::RGIGlacier)::String = gid.attrs.name
Base.show(stream::IO, gid::RGIGlacier) = print(stream, "RGIGlacier($(gid.id),...)")

# RGI specifics:
function parse_rgi(rgi::AbstractString)
    vers = VersionNumber(parse(Int,split(rgi,'-')[1][4:end])÷10)
    reg = parse(Int,split(split(rgi,'-')[2],'.')[1])
    nr = parse(Int,split(split(rgi,'-')[2],'.')[2])
    @assert 0<reg<20 "Not a valid RGI region!"
    return vers, reg, nr
end
parse_rgi(gid::RGIGlacier) = parse_rgi(string(gid.id))
getrgi_version(rgi) = parse_rgi(rgi)[1]
getrgi_region(rgi) = parse_rgi(rgi)[2]
getrgi_nr(rgi) = parse_rgi(rgi)[3]

include("rgi.jl")


"""
Antarctic Peninsula glaciers with ID as Huss&Farinotti uses.
"""
struct PeninsulaGlacier <: GlacierID
    id::Int
end
getname(gid::PeninsulaGlacier) = string(gid.id)
getregion(::PeninsulaGlacier) = Region.AntarcticPeninsula
"Data loader from H&F Peninsula"
const PeninsulaLoader = :Peninsula
include("peninsula.jl")


"""
GloGEM bdot data-loader

Glacier ID according to RGI 6.0
"""
const GloGEMLoader = :GloGEM
include("glogem.jl")


"""
GlaThiDa thickness data-loader

Glacier ID according to RGI 6.0
"""
const GlaThiDaLoader = :GlaThiDa
include("glathida.jl")

"""
Misc CCI loaders (IV and dh/dt)

"""
const CCILoader = :CCI
include("CCI.jl")


"""
Glamos data-loader

(No glaciers associated with it)
"""
const GlamosLoader = :GLAMOS

"""
No data loader

Always returns SomeData(nothing)
"""
const VoidLoader = :VOID
load!(ds::DataSet{<:Any, VoidLoader}, varargs...) = SomeData(nothing)
