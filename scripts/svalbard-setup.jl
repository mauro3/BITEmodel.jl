# common setup for Svalbard runs
#
# Todo

println("Time to load BITEModel.jl")
@time import BITEModel
const BM=BITEModel

update_cache = false
verbose = false
plotyes = false

# RGI data
all_svalbard = readdlm("data/rgi-6.0-huss/outlines/xyzn/svalbard/svalbard.dat")
cutoff_km2 = (5,10^10) # km^2
westfonna_ll = [18.5, 79.7] # long lat
austfonna_ll = [20, 79.2] # long lat
svalbard_nrs = Int[all_svalbard[i,1] for i=1:size(all_svalbard,1) if
                   cutoff_km2[1]<=all_svalbard[i,end]<=cutoff_km2[2] && # sizes
                   # remove Austfonna and Westfonna:
                   !(all_svalbard[i,2]>westfonna_ll[1] && all_svalbard[i,3]>westfonna_ll[2]) &&
                   !(all_svalbard[i,2]>austfonna_ll[1] && all_svalbard[i,3]>austfonna_ll[2])]
#with_thick = ["00027", "00408", "00409", "00427", "00428", "00429", "00434"][2:end] # 27 is on west or austfonna
with_thick = ["00409", "00428"] # the others are useless

brucebreen = BM.RGIGlacier("RGI60-07.01455", "Brucebreen") # area 6km^2, land, long lat: 17.291 78.485
rabotbreen = BM.RGIGlacier("RGI60-07.01478", "Rabotbreen") # area 70km^2, land
tunabreen = BM.RGIGlacier("RGI60-07.01458", "Tunabreen") # area 160km^2, tide-water
region_svalbard = 7

init_svalbard_gl(nr::Int; kws...) = init_svalbard_gl(BM.RGIGlacier(nr, region_svalbard); kws...)
init_svalbard_gl(rgi::Union{String,Symbol}; kws...) = init_svalbard_gl(BM.RGIGlacier(rgi); kws...)
function init_svalbard_gl(gid::BM.GlacierID; kws...)
    dt, pl = init_svalbard_gl_dt(gid; kws...)
    init_svalbard_gl(dt, pl)
end

function init_svalbard_gl(dt::BM.DataTable, pl)
    gl,pp,pm,pn,pl = BM.load_glacier(dt,pl)
    gb, pm = BM.make_bands(gl, pp, pm, pn);
    return gl,gb,pp,pm,pn,pl,dt
end
function init_svalbard_gl_dt(gid; kws...)
    pl = BM.LoadPara(gid; kws...)
    dt, pl = BM.make_datatable(gid,pl)

    # set temperature profile:
    dt.temp.opts[:temp] = [-10, -1.0] # colder at tongue than at top
    dt.temp.opts[:ele] = linspace(-30.0, 1717, 2) # a bit below sea to height of highest peak on Svalbard
    return dt, pl
end
