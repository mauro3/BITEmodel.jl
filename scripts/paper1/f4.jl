import BITEModel
const BM=BITEModel
using Plots, StatPlots, Measures, NamedTuples

#include("itmix-comp-fit-sigma.jl")
# include("../rgi-setup.jl")
# region = [3, 4, 14, 11, 7][end]
# runtyp = "test"
# update_cache = true
# glnrs1, glnrs2, glnrs_all = get_glaciers(region)
# sigma_of_model = @NT(sigma_h_model = 10.0,
#                      sigma_iv_model = 10.0)
# fit_sigma = false
# nr = glnrs2[end]

# # number of divisions in elevation band variables
# fit_vars = @NT(n_1d_btilde=3,
#                n_1d_fsl=3,
#                n_1d_temp=1)
# const progress_bar = true
# fit_target = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.length][1]
# out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(1), dir_output="", retsol=true, fit_sigma=fit_sigma)


fl = [# not working anymore because of type updates:
    # "figs/f2-5-34-42-h.jls",
      # "figs/f2-5-34-42-h_iv-abserr.jls",
      # "figs/f2-5-34-42-h-abserr.jls",
      # "figs/f2-5-34-42-h_iv-abserr-sigma20.jls",
      # "figs/f2-5-34-42-h_iv-abserr-thin3.jls",
      # "figs/f2-5-34-42-h-abserr-sigmah20-thin3.jls", # looks sweet
      # "figs/f2-5-34-42-h-abserr-thin3.jls",
      # "figs/f2-5-34-42-h-abserr-sigma30-thin1-lengths1e-3.jls",
      # "figs/f2-5-34-42-h-sqrerr-sigma20-thin1-lengths1e-3.jls",
      # "figs/f2-mar20-5-34-42-h-abserr-sigma30-thin1.jls",
      # "figs/f2-mar20-5-34-42-h-sqrerr-sigma40-thin1.jls",
      # "figs/f2-mar21-5-34-42-h-sqrerr-sigma50-thin1.jls",
      # "figs/f2-mar21-5-34-42-h-sqrerr-sigma60-thin1.jls",
      # "figs/f2-mar21-5-34-42-h-sqrerr-sigma70-thin1.jls",
      # "figs/f2-mar21-5-34-42-h-abserr-sigma20-thin1.jls",
      # "figs/f2-mar22-5-34-42-h-abserr-sigma20-thin1.jls",
      "figs/f2-apr01-5-34-42-h-sqrerr-sigma20-thin1.jls",
      # fitting sigma too
      "figs/f2-with_sigma-apr03-5-34-42-h-sqrerr.jls", # with Gaussian sigma prior
      "figs/f2-with_sigma-apr04-5-34-42-h_iv-sqrerr.jls", # Uniform
      "figs/f2-with_sigma-apr04-5-34-42-h-sqrerr.jls", # Uniform: very similar to Gaussian
      "figs/f2-with_sigma-apr04-5-34-42-iv-sqrerr.jls", # Uniform: doesn't actually look that bad...
      # noveau
      "figs/f2-with_sigma-better-apr16-5-34-42-h_iv-sqrerr.jls", # with sigma_btilde = 0.3
      "figs/f2-large-btilde-sigma-apr17-5-34-42-h_iv-sqrerr.jls", # with sigma_btilde = 2.8
      "figs/f2-sigma-btilde-1.4-sigma-apr17-5-34-42-h_iv-sqrerr.jls",  # with sigma_btilde = 1.4 (production value)
      #"figs/f2-sigma-btilde-1.4_sigma-dhdt-26-sigma-apr17-5-34-42-h_iv-sqrerr.jls",  # WRONG with sigma_btilde = 1.4 (production value), =26 (production value)
      "figs/f2-production-sigma-apr17-5-34-42-h_iv-sqrerr.jls",
      "figs/f2-production-length97-sigma-apr17-5-34-42-h_iv-sqrerr.jls", # <--- in publication
      "figs/f2-production-length98-sigma-apr17-5-34-42-h_iv-sqrerr.jls", # BAD, something is off here, and for sure not converged.
      "figs/f2-production-length98-btildes-0.3-apr17-5-34-42-h_iv-sqrerr.jls", # back to small btilde error.  Sigma h-error is large
      "figs/f2-production-length97-btildes-0.3-ball-apr17-5-34-42-h_iv-sqrerr.jls", #  <--- in publication, h-error is small
      "figs/f2-production-length98-btildes-0.3-ball-apr17-5-34-42-h_iv-sqrerr.jls", # Sigma h-error is large
      "figs/f2-production-length995-btildes-0.3-ball-apr17-5-34-42-h_iv-sqrerr.jls", # trying super small length sigmab: looks very good apart from super large sigmas, because the lowest profile is not well fitted
      ][end-2]
@show fl

if true # !isdefined(:sol)
    Revise.track(Base) # base/serialize.jl L971 need to modify it
    sol, sol_prior = open(deserialize, fl, "r")

    gl = sol.theta0.gb.gl
end

if !isdefined(:gl_all_radar)
    glacier = :Unteraar
    gid = BM.ITMIXGlacier(glacier, 2)
    gl_all_radar,gb,pp,pm,pn,pl = BM.init_forward(gid)
end

err, bed, surf, rad = BM.error_h(sol.mh2d, gl_all_radar);
err_max, bed_max = BM.error_h(sol.max_h2d, gl_all_radar);
err_min, bed_min = BM.error_h(sol.min_h2d, gl_all_radar);
_, _,_, _, bed_std = BM.error_h(sol.sh2d, gl_all_radar);

for (ii,sps) in enumerate([[7, 16, 25, 38, 45], [5, 34, 42, 16, 25]]) # fit to 5, 34, 42
    # map view
    p = plot(aspect_ratio=:equal);
    for (ss,s) in enumerate(err.splits)
        plot!(err.x[s], err.y[s], label="", lw=:0.5, lc=:black)
        #annotate!(err.x[s][1], err.y[s][1], Text("$ss"))
    end

    # fitting tracks
    for (ss,s) in enumerate(gl.h.vals.splits)
        plot!(gl.h.vals.x[s], gl.h.vals.y[s], lw=3, label="", lc=:black)
    end


    str = ["A" "B" "C" "D" "E"]
    for (i,ss) in enumerate(sps)
        s = err.splits[ss]
        plot!(gl_all_radar.h.vals.x[s], gl_all_radar.h.vals.y[s], lw=2, label="", lc=:red)
        annotate!(err.x[s][1], err.y[s][1], text("$(str[i])", 15), lc=:green)
    end
    BM.plot_outlines!(gl_all_radar, xyscale=1, lw=0.5, xlabel="", ylabel="", grid=false, framestyle=:none)
    plot!([440e3], [5155e3], xerror=[-1e3,1e3], markercolor=:black, markerstrokecolor=:black, markerstrokewidth=2,
          markersize=10, label="", left_margin=10mm; ((:xlims, (432.2e3,442.7e3)), (:ylims, (5153.8e3, 5162.2e3)))...)
    annotate!(440e3, 5155.4e3, text("2 km", 10))

    #display(p)

    ps = []
    for (i,ss) in enumerate(sps)
        onleft = i==1 || i==4
        onbottom = i>3
        s = err.splits[ss]
        dist = cumsum([0, sqrt.(diff(err.x[s]).^2 + diff(err.y[s]).^2)...])
        push!(ps,
              plot(dist, surf.v[s], lw=2, color=:green,  label= i==length(str) ? "surface" : ""))
        plot!(dist, rad.v[s], lw=2, label= i==length(str) ? "radar bed" : "", markershape=:circle, markersize=1.5,
              color=:red, markercolor=:black)
        plot!(dist, rad.v[s]+27, linecolor=:red, linewidth=0,
              fillrange=rad.v[s]-27, fillcolor=:red, fillalpha=0.3,
              label= i==length(str) ? "   ± 27m" : "")
        plot!(dist, bed.v[s], lw=1.5, linecolor=:blue, label= i==length(str) ? "mean bed" :"",
              xlabel= onbottom ? "distance (m)" : "", ylabel= onleft ? "elevation (m a.s.l.)" : "")
        plot!(dist, bed.v[s]+2bed_std.v[s], linecolor=:blue, linewidth=0,
              fillrange=bed.v[s]-2bed_std.v[s], fillcolor=:blue, fillalpha=0.3,
              label= i==length(str) ? "   2x std" :"")
        plot!(dist, bed_max.v[s], lw=0.5, color=:black, label= i==length(str) ? "    extrema" : "")
        plot!(dist, bed_min.v[s], lw=0.5, color=:black, label="")
        # if i==length(str)
        #     # make labels
        #     x, y = [err.x[1:1],dist[1:1]][2] , [NaN]
        #     # plot!(x,y, lw=2, label="surface")
        #     # plot!(x,y, lw=2, label="radar bed ± 10m")
        #     plot!(x,y, lw=1.5, label= "mean")
        #     plot!(x,y, fillrange=y, linecolor=:blue, linewidth=0,
        #           fillcolor=:blue, fillalpha=0.3,
        #           label="2x std")
        #     plot!(x,y, lw=0.5, label="max")
        #     plot!(x,y, lw=0.5, label="min")
        # end

    end

    l = @layout [grid(1,3)
                 grid(1,3)] # a{0.66w}]
    kwargsj = [(:tickfontsize, 10), (:legendfontsize, 10),
               (:xguidefontsize, 12), (:yguidefontsize, 12), (:titlefontsize, 20)]

    pend = plot(ps..., p, layout=l, title=["A" "B" "C" "D" "E" "F"], title_location=:right, size=(900,600), reuse=false; kwargsj...)
    display(pend)
    if ii==1
        savefig("figs/f4.png")
        savefig("figs/f4.pdf")
    else
        savefig("figs/f4-fitting-profiles.png")
        savefig("figs/f4-fitting-profiles.pdf")
    end

end

#BM.plottheta(sol.thetas, sol.theta0, reuse=false)
