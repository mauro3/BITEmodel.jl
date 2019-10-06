# Run forward model demonstrating hand-adjusting some stuff
println("Time to load BITEModel.jl")
@time using BITEModel
const BM=BITEModel
if !isdefined(:plotyes)
    plotyes = true
end

# run forward model for ITMIX glaciers
glaciers = [:Academy, :Aqqutikitsoq, :Austfonna, :Brewster, :Columbia,
            :Devon, :Elbrus, :Freya, :Hellstugubreen, :Kesselwandferner, :Mocho,
            :NorthGlacier, :SouthGlacier, :Starbuck, :Synthetic1, :Synthetic2,
            :Synthetic3, :Tasman, :Unteraar, :Urumqi, :Washmawapta]
# glacier = glaciers[gln]
glacier = :Unteraar
use_glogem = false # use GloGEM bdot, if available, and use the SyntheticBenchLoader for dhdt
update_cache = false

##################
gid = BM.ITMIXGlacier(glacier, 2)

println("Time to to make load-paras")
@time pl = BM.LoadPara(gid, update_cache=update_cache, use_glogem=use_glogem)

# Now customize pl.  This essentially follows what is encoded in the
# functions init_forward and the custom LoadPara method.

# E.g. make sigma larger
pl.dataset_opts[BM.ThicknessData][:sigma]=40.0

if glacier==:Elbrus
    # why not some more down-sampling?
    pl.dataset_opts[BM.DataKind][:grid_downsample_step]=3
end

pl.dataset_opts[BM.ParaData][:pm][:bandsize] = 40 # bigger band

dt, pl = BM.make_datatable(gid,pl)
# swap out the bdot against one from the synthetic bench run
bdot_synth = BM.DataSet{BM.BdotData, BM.SyntheticBenchLoader}([""],
                                                              pl.dataset_opts[BM.BdotData],
                                                              (Date(),Date())
                                                              )
dt = BM.DataTable(dt, bdot=bdot_synth)

gl,pp,pm,pn,pl = load_glacier(dt,pl)
println("Time for loading glacier:")
@time gl,pp,pm,pn,pl = load_glacier(dt,pl)

Base.Test.@inferred BM.make_bands(gl, pp, pm, pn);
println("Time for making bands:")
@time gb,pm = BM.make_bands(gl, pp, pm, pn);

Base.Test.@inferred BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

println("Time to run forward model:")
gc()
@time fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

qtot1d, qd1d = postproc1d(gb, pp, pm)

if plotyes
    println("Time to plot:")
    @time display(BM.plot_run(fwdsol, reuse=false, plot2d=[:h,:iv,:ivlog][1]))
end
