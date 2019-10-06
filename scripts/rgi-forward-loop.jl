# Run forward model

# Parameters to adjust:
#######################

update_cache = true
verbose = true
plotyes = true

@show region = [3, 4, 14, 11, 7][end]

include("rgi-setup.jl")
gl_nrs = get_glaciers(region)[1]

# Run them
##########
errs = []
area_vols = zeros(length(gl_nrs), 4) # number, area, volume-area, volume inferred [mk3

start_t = time()
time_fwd = 0.0
for (i,nr) in enumerate(gl_nrs)
    area_vols[i, 1] = nr
    #nr in merge(faulty,strange) && continue
    println("\n")
    try
        gid = BM.RGIGlacier(nr,region)
        @show gid
        pl = BM.LoadPara(gid; update_cache=update_cache)
        print("Time for load:")
        @time gl,pp,pm,pn,pl = BM.load_glacier(gid,pl)
        area_vols[i, 2] = BM.area(gl) / 1e3^2
        area_vols[i, 3] = BM.volume_area_mean_h(gl) * BM.area(gl) /1e3^3

        print("Time for bands:")
        @time gb,pm = BM.make_bands(gl, pp, pm, pn)

        # gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
        #                                     update_cache=update_cache)
        pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01)
        print("Time to run forward model:"); gc()
        @time time_fwd += @elapsed fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

        BM.postproc1d!(fwdsol)
        area_vols[i, 4] = BM.volume(fwdsol)/1e3^3

        if plotyes
            print("Time to plot:")
            @time Plots.savefig(BM.plot_run(fwdsol, reuse=true, plot2d=[:h,:iv,:ivlog][1]),
                                "figs/$(gid.id).png")
        end
    catch e
        e==InterruptException() && throw(e)
        @show e
        push!(errs, (nr, e))
    end
end
dt = round(Int, (time() - start_t)/60) # time in min
time_fwd = round(Int, time_fwd)
dirname = infos[region][1]
writecsv("figs/err-$dirname-t$dt-f$(time_fwd).csv", errs)
writecsv("figs/vol-$dirname-t$dt-f$(time_fwd).csv", area_vols)

if false
    gid = BM.RGIGlacier(nr,region);
    pl = BM.LoadPara(gid; update_cache=update_cache);
    @time gl,pp,pm,pn,pl = BM.load_glacier(gid,pl);
    @time gb,pm = BM.make_bands(gl, pp, pm, pn);
    pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01);
    @time time_fwd += @elapsed fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp);
    BM.postproc1d!(fwdsol);
    BM.plot_run(fwdsol, reuse=true, plot2d=[:h,:iv,:ivlog][1])
end
