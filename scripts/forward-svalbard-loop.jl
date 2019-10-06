# Run forward model

println("Time to load BM.jl")
@time using BITEModel
const BM=BITEModel
using Parameters
import PyPlot

# Parameters to adjust:
#######################

update_cache = true
verbose = true

glaciertype = [:svalbard, :itmix, :synth, :peninsula, :rgi][1]
plotyes = true

# a few test glaciers:
region_svalbard = 7
all_svalbard = readdlm("data/rgi-6.0-huss/outlines/xyzn/svalbard/svalbard.dat")
cutoff = (0,10^10) # km^2
westfonna_ll = [18.5, 79.7]
austfonna_ll = [20, 79.2]
svalbard_nrs = Int[all_svalbard[i,1] for i=1:size(all_svalbard,1) if
                   cutoff[1]<=all_svalbard[i,end]<=cutoff[2] && # sizes
                   # remove Austfonna and Westfonna:
                   !(all_svalbard[i,2]>westfonna_ll[1] && all_svalbard[i,3]>westfonna_ll[2]) &&
                   !(all_svalbard[i,2]>austfonna_ll[1] && all_svalbard[i,3]>austfonna_ll[2])]
with_radar = ["00027", "00408", "00409", "00427", "00428", "00429", "00434"][2:end] # 27 is on west or austfonna
with_good_radar = [409, 428] # ELFENBEINBREEN and SVEIGBREEN

brucebreen = BM.RGIGlacier("RGI60-07.01455", "Brucebreen") # area 6km^2, land, long lat: 17.291 78.485
ll = (17.291, 78.485)
rabotbreen = BM.RGIGlacier("RGI60-07.01478", "Rabotbreen") # area 70km^2, land
tunabreen = BM.RGIGlacier("RGI60-07.01458", "Tunabreen") # area 160km^2, tide-water
whatbreen = BM.RGIGlacier("RGI60-07.00434") # area 8km^2
gid = [whatbreen, rabotbreen, brucebreen, tunabreen][1]


# Run it
########
probs = [118]
errs = []
# nr =222
#         gid = BM.RGIGlacier(nr,region_svalbard)
#         gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
#                                             update_cache=update_cache)
#         pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01)
#         print("Time to run forward model:"); gc()
#         @time hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
#             BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

#         qtot1d, qd1d = postproc1d(gb, pp, pm)

faulty = [319,440,]
strange = [849,857,892,944] # dem wrong 3x, dem boring 944

for nr in svalbard_nrs
    #nr in merge(faulty,strange) && continue
    println("\n")
    try
        gid = BM.RGIGlacier(nr,region_svalbard)
        @show gid
        gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
                                            update_cache=update_cache)
        pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01)
        print("Time to run forward model:"); gc()
        @time fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

        BM.postproc1d!(fwdsol)

        if plotyes
            print("Time to plot:")
            p = BM.plot_run(fwdsol, reuse=true, plot2d=[:h,:iv,:ivlog][1])
            Plots.savefig(p, "figs/$(gid.id).png")
        end
    catch e
        e==InterruptException() && throw(e)
        @show e
        push!(errs, (nr, e))
    end
end
