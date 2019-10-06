import BITEModel
const BM=BITEModel
using Plots, StatPlots, Measures

# fl = ["figs/f2-5-34-42-h.jls", "figs/f2-5-34-42-h_iv-abserr.jls"][2]
# sol, sol_prior = open(fl, "r") do io
#     deserialize(io)
# end
# gl = sol.theta0.gb.gl

# glacier = :Unteraar
# gid = BM.ITMIXGlacier(glacier, 2)
# gl_all_radar,gb,pp,pm,pn,pl = BM.init_forward(gid)


erriv = BM.error2d(sol.miv2d, gl, :iv, gl.ivmask)
error()

err, bed, surf, rad = BM.error_h(sol.mh2d, gl_all_radar);
err_max, bed_max = BM.error_h(sol.max_h2d, gl_all_radar);
err_min, bed_min = BM.error_h(sol.min_h2d, gl_all_radar);
_, bed_std = BM.error_h(sol.sh2d, gl_all_radar);

# map view
p = plot(aspect_ratio=:equal);
# make labels
x, y = err.x[1:1], [NaN]
plot!(x,y, lw=1.5, label= "mean")
plot!(x,y, lw=0.5, label="max")
plot!(x,y, lw=0.5, label="min")
plot!(x,y, lw=2, label="radar")
plot!(x,y, lw=2, label="surface")
for (ss,s) in enumerate(err.splits)
    plot!(err.x[s], err.y[s], label="", lw=:0.5, lc=:black)
    #annotate!(err.x[s][1], err.y[s][1], Text("$ss"))
end

# fitting tracks
for (ss,s) in enumerate(gl.h.vals.splits)
    plot!(gl.h.vals.x[s], gl.h.vals.y[s], lw=2, label="", lc=:black)
end

sps = [5, 16, 25, 35, 42]
sps = [5, 34, 42, 16, 25]
str = ["A" "B" "C" "D" "E"]
for (i,ss) in enumerate(sps)
    s = err.splits[ss]
    annotate!(err.x[s][1], err.y[s][1], Text("$(str[i])"), lc=:green)
end
BM.plot_outlines!(gl_all_radar, xyscale=1, lw=0.5, xlabel="", ylabel="", grid=false, framestyle=:none)
#display(p)

ps = []
for ss in sps
    s = err.splits[ss]
    dist = cumsum([0, sqrt.(diff(err.x[s]).^2 + diff(err.y[s]).^2)...])
    push!(ps, plot(dist, bed.v[s], lw=1.5, label= "", xlabel="distance (m)", ylabel="z (m asl)"))
    plot!(dist, bed_max.v[s], lw=0.5, label="")
    plot!(dist, bed_min.v[s], lw=0.5, label="")
    plot!(dist, rad.v[s], lw=2, label="")
    plot!(dist, surf.v[s], lw=2, label="")
end

l = @layout [grid(1,3)
             grid(1,3)] # a{0.66w}]
plot(ps..., p, layout=l, title=["A" "B" "C" "D" "E" "F"], title_location=:left, size=(900,600), reuse=false)

#savefig("figs/f4.png")
