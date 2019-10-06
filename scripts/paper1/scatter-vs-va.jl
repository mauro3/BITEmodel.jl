using Plots
outs = (out1,out2);

nrs = ([o isa Dict ? BM.parse_rgi(o[:rgi])[3] : -1 for o in out1],
       [o isa Dict ? BM.parse_rgi(o[:rgi])[3] : -1 for o in out2])

v = ([o isa Dict ? o[:vol_km3][:median] : NaN for o in out1],
     [o isa Dict ? o[:vol_km3][:median] : NaN for o in out2])

a = ([o isa Dict ? o[:area_km2] : NaN for o in out1],
     [o isa Dict ? o[:area_km2] : NaN for o in out2])

h = (1000*v[1]./a[1], 1000*v[2]./a[2])

va = ([o isa Dict ? o[:volarea_km3] : NaN for o in out1],
      [o isa Dict ? o[:volarea_km3] : NaN for o in out2])

hva = (1000*va[1]./a[1], 1000*va[2]./a[2])

ratioh = (h[1]./hva[1], h[2]./hva[2])

# radar/fitting
fls = filter(x->startswith(x, region_str), readdir("data/GlaThiDa_3.0/mhuss/"))
h_nrs = map(x->parse(Int, split(split(x,'.')[1],'_')[2]), fls)


p1 = Plots.scatter(hva[1], h[1], label="with radar", xlabel="mean h of vol-area (m)",
                   ylabel= "mean h BITE (m)", reuse=false, title="Region $region_str") #, ticks=:native
top = max(VAWTools.maximumnan(hva[1]), VAWTools.maximumnan(h[1]))
Plots.plot!(linspace(0, top), linspace(0, top), label="1:1")
Plots.plot!(linspace(0, top), linspace(0, 2*top), label="1:2")
Plots.plot!(linspace(0, top), linspace(0, 0.5*top), label="1:0.5")

p2 = Plots.scatter(hva[2], h[2], label="without radar", xlabel="mean h of vol-area (m)",
                   ylabel= "mean h BITE (m)", reuse=false) # , ticks=:native
top = VAWTools.maximumnan(hva[2])
Plots.plot!(linspace(0, top), linspace(0, top), label="1:1")
Plots.plot!(linspace(0, top), linspace(0, 2*top), label="1:2")
Plots.plot!(linspace(0, top), linspace(0, 0.5*top), label="1:0.5")
Plots.plot!(linspace(0, top), linspace(0, 5*top), label="1:5")
Plots.plot!(linspace(0, top), linspace(0, top/5), label="1:0.2")

#display(plot(p1, p2, layout=(2,1)))
display(p1); display(p2)

#savefig("figs/scatter-$(region_str)_$(split(dir_,'/')[end-1]).png")

### Find outliers above

    fac = 5.0
    outl = []
    outln = []
    for i=1:2
        inds = ratioh[i].>fac
        push!(outl, [(n,ii,o) for (n,ii,o) in zip(nrs[i][inds], (1:length(outs[i]))[inds], outs[i][inds])])
        push!(outln, nrs[i][inds])
    end
    @show length.(outl), length.(nrs)

    # find outliers below
    fac = 0.5
    outlb = []
    for i=1:2
        inds = ratioh[i].<fac
        push!(outlb, [(n,o) for (n,o) in zip(nrs[i][inds], outs[i][inds])])
    end
    @show length.(outlb), length.(nrs)

if false
    mode = 1
    ii = length(outl)
    update_cache=false
    nr, ind, radar, o = outl[mode][ii];
    @show nr, ii
    @show radar, h[mode][ind], hva[mode][ind]

    (gid, gl,gb,pp,pm,pn,pl,th0d,
     theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) = BM.init(o[:rgi], fit_vars, fit_target, sigma_of_model, runtyp, update_cache=update_cache);
    fwdsol = fwdm_fn(o[:thetas_mode]);
    display(BM.plot_run(fwdsol))
    display(Plots.plot(gl, field=gl.dhdt.vals, reuse=false))
    display(Plots.plot(gb.fields1d_original[:dhdt], reuse=false))

end
nothing

# Region 7, mode 1 TOO BIG, mode 1
#  261 odd dhdt
#  437 ?
#  532 synthetic dhdt is bad: positive in ablation area?! -> fixed dhdt.  But still bad...; not fitting to IV does not improve it.
#  657 synthetic dhdt is bad: positive in ablation area?!
#  664 synth, smb>0
#  858 synth but ok dhdt ???
#  887 synth, smb>0
#  890 has no ablation area, not calving?!
# > 893 ??? IV problem?
#  895 synth, smb>0
#  896 no ablation
# > 902 no ablation, odd DEM, --> use other DEM
#  903 high SMB
#  932 SMB>0, iscalving==false
#  951 synthetic dhdt is bad: positive in ablation area?!
#  984 synthetic dhdt is bad: positive in ablation area?!
# > 1396 The fitting is bad here with the inferred dhdt opposite to the measured one?! IV?
# 1464 (Kronebreen) has no radar?! (only one point)
# 1471 fixed(?)
# 1476
# 1492 Julibreen
# 1541
# 1542
# 1558
