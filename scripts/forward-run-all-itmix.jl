# Runs the forward model of all ITMIX glaciers, updates cached files and can plot all.
@time using BITEModel
const BM=BITEModel
verbose = false

itmixversion = 2
use_glogem = false
update_cache = true
plotting = false

# make them available in case the loop craps out
global glacier, gid, gl,gb,pp,pm,pn,pl,hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio

# these don't load yet
toskip = [] # [:Brewster, :Devon, :Synthetic2, :Synthetic3, :Unteraar, :Urumqi, :Washmawapta #=only v1=#]
for (ii,glacier) in enumerate(BM.ITMIXnames)
    #ii<9 && continue

    glacier in toskip && continue
    itmixversion==1 && glacier in BM.itmix2only && continue
    println(" ")
    println(glacier)
    gid = ITMIXGlacier(glacier, itmixversion)
    print("Loading ")
    # gid = ITMIXGlacier(glacier, 2); pl = BM.LoadPara(gid); dt, pl = BM.make_datatable(gid, pl);  ds = dt.outline;
    # outline, iscalving, isicecap, outline_proj =  BM.load!(ds, gid)
    @time gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
                                              update_cache=update_cache, use_glogem=use_glogem);

    print("Running fwd-model ")
    @time hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
        BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp);
    @show prod(size(gl.dem))
    if plotting
        println("Plotting")
        display(BM.plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, pp, pm,
                            reuse=false, plot2d=[:h,:iv,:ivlog][1]))
    end
end
