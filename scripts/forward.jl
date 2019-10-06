# Run forward model

println("Time to load BITEModel.jl")
@time using BITEModel
const BM=BITEModel
const I2=BM.ITMIX2

# Parameters to adjust:
#######################

glaciertype = [:itmix, :synth, :peninsula, :rgi][2]
if !isdefined(:plotyes)
    plotyes = true
end

if glaciertype==:itmix
    # run forward model for ITMIX glaciers
    glaciers = [:Academy, :Aqqutikitsoq, :Austfonna, :Brewster, :Columbia,
                :Devon, :Elbrus, :Freya, :Hellstugubreen, :Kesselwandferner, :Mocho,
                :NorthGlacier, :SouthGlacier, :Starbuck, :Synthetic1, :Synthetic2,
                :Synthetic3, :Tasman, :Unteraar, :Urumqi, :Washmawapta]
    glacier = :Unteraar
    pl_kws = Dict{Symbol,Any}(:use_glogem => false)
    gid = ITMIXGlacier(glacier)
elseif glaciertype==:synth
    pl_kws = Dict{Symbol,Any}()
    if haskey(pl_kws, :use_glogem)
        pop!(pl_kws, :use_glogem)
    end
    gid = SyntheticGlacier(:bench)
    pl_kws = Dict{Symbol,Any}()
elseif glaciertype==:peninsula
    flaskid,starbuckid = BM.PeninsulaGlacier(848), BM.PeninsulaGlacier(900)
    gid = [flaskid,starbuckid][2]
    pl_kws = Dict{Symbol,Any}()
else
    error("not implemented")
end
update_cache = false
verbose = true

# Run it
########
gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
                                    update_cache=update_cache, pl_kws...)
pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01)
print("Time to run forward model:"); gc()
@time fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

qtot1d, qd1d = postproc1d(gb, pp, pm)

if plotyes
    print("Time to plot:")
    @time display(BM.plot_run(fwdsol, reuse=false, plot2d=[:h,:iv,:ivlog][1]))
end
