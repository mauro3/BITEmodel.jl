# Runs through all ITMIX glaciers, updates cached files and plots all.
using BITEModel
const BM=BITEModel
plotyes = false
global glacier, gl,gb,pp,pm,pn,pl,hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio
gl_with_errors = []
gl_without_errors = []
glns = 1:1606
t_tot = @elapsed for id in glns
    update_cache = false
    pl_kws = Dict{Symbol,Any}()
    verbose = false
    println(" ")
    println("Penisula glacier $id")
    gid = BM.PeninsulaGlacier(id)
    print("Loading ")
    try
        pl = BM.LoadPara(gid, update_cache=update_cache)
        pl.dataset_opts[BM.ParaData][:pn][:calc_boxcar_M] = false
        print("Loading ... ")
        @time gl,pp,pm,pn,pl = load_glacier(gid,pl)
        print("Elevation bands ... ")
        @time gb,pm = BM.make_bands(gl, pp, pm, pn)

        print("Running fwd-model ")
        @time hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
            BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
        mh = round(Int,sum(hs1d[1:end-1].*gb.ws.*gb.ls)/BM.area(gl))
        mhva = round(Int, BM.volume_area_mean_h(BM.area(gl), isicecap=true))
        println("Mean h: $mh ; Volume-area: $mhva; Glacier area: $(BM.area(gl)/1e6)km^2")
        push!(gl_without_errors, (id, mh, mhva, BM.area(gl)/1e6))
        # println("Plotting")
        # display(BM.plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, pp, pm,
        #                     reuse=false, plot2d=[:h,:iv,:ivlog][1]))
    catch e
        e==InterruptException() && break
        push!(gl_with_errors, (id, e))
        warn("ERROR $e")
        println("")
    end
end

println("\n\nTotal time $(round(t_tot/60,2)) minutes; average per glacier $(round(t_tot/length(glns),3))s.")
println("$(length(gl_with_errors)) out of $(length(glns)) glaciers had errors")
