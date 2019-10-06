# Run forward model

# Load it
#########
include("svalbard-setup.jl")
plotyes = true
gl,gb,pp,pm,pn,pl,dt = init_svalbard_gl(1081, update_cache=true) #tunabreen
@show gl.gid
@show tunabreen

# Run it
########
pm = BM.MPara(pm, iv_nye=true, iv_dist_exp2=0.01)
print("Time to run forward model:"); gc()
@time fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

qtot1d, qd1d = BM.postproc1d(gb, pp, pm)

if plotyes
    print("Time to plot:")
    @time display(BM.plot_run(fwdsol, reuse=false, plot2d=[:h,:iv,:ivlog][1]))
end
