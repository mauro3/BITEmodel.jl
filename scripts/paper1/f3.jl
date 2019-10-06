import BITEModel
const BM=BITEModel
using Plots, StatPlots, Measures, Parameters

# fl = ["figs/f2-5-34-42-h.jls"][1]
# sol, sol_prior = open(fl, "r") do io
#     deserialize(io)
# end
# gl = sol.theta0.gb.gl


function plotinv1d_v(sol; kwargs...)
    @unpack mc_hs1d, mc_ivs1d = sol
    gb = sol.theta0.gb
    ivs1d_measured = gb.fields1d_original[:iv]
    xl = (gb.x[1]/1e3,gb.x[end]/1e3)
    hs1d = mean(mc_hs1d,2)
    hs1d_max = maximum(mc_hs1d,2)
    hs1d_min = minimum(mc_hs1d,2)
    hs1d_std = std(mc_hs1d,2)

    ivs1d = mean(mc_ivs1d,2)
    ivs1d_max = maximum(mc_ivs1d,2)
    ivs1d_min = minimum(mc_ivs1d,2)
    ivs1d_std = std(mc_ivs1d,2)

    # p1 = plot(gb.x/1e3, gb.bands_; label="surface")#, title="Flow line")
    # plot!(gb.x/1e3, gb.bands_-hs1d_max, linecolor=:red, linewidth=0,
    #           fillrange=gb.bands_-hs1d_min, fillcolor=:red, fillalpha=0.5, label="min/max predicted bed",
    #           ylabel="elevation (m)")
    # plot!(gb.x/1e3, gb.bands_-hs1d, label="mean predicted bed", xlims=xl, ylims=ylims1, linecolor=:black)
    # if ele_only
    #     return p1
    # end

    # p2 = plot(gb.x/1e3, hs1d_max,  linecolor=:red, linewidth=0,
    #           fillrange=hs1d_min, fillcolor=:red, fillalpha=0.3,
    #           ylabel="h (m)", label="min/max")
    # plot!(gb.x/1e3, hs1d+hs1d_std,  linecolor=:blue, linewidth=0,
    #       fillrange=hs1d-hs1d_std, fillcolor=:blue, fillalpha=0.5,
    #       label="std")

    # plot!(gb.x/1e3, hs1d, label="mean", xlims=xl, ylims=ylims2, linecolor=:black)

    p3 = plot(gb.x/1e3, ivs1d_max,  linecolor=:red, fillrange=ivs1d_min,
              fillcolor=:red, fillalpha=0.5, label="min/max", linewidth=0,
              ylabel="v (m/a)", xlabel="χ (km)", xlims=xl)
    plot!(gb.x/1e3, ivs1d+ivs1d_std,  linecolor=:blue, linewidth=0,
          fillrange=ivs1d-ivs1d_std, fillcolor=:blue, fillalpha=0.5,
          label="std")
    plot!(gb.x/1e3, ivs1d, label="mean", linecolor=:black)
    # im = copy(ivs1d_measured)
    # im[im.<0] = NaN
    # im[end]=NaN
    # # scatter!(gb.xmid/1e3, im, label="v 1D measured", xlims=xl)


    # plot(p1,p2,p3;layout=(3,1), size=(750,900), kwargs...)
end


p1 = BM.plotinv1d(sol, reuse=false, ele_only=true, xlabel="χ (km)", ylabel="elevation (m a.s.l.)");
annotate!(0.5,2000, text("A", 25), xlabel="χ (km)", ylabel="elevation (m a.s.l.)");

p2 = plotinv1d_v(sol, reuse=false, xlabel="χ (km)");
annotate!(0.5,30, text("B", 25), xlabel="χ (km)", ylabel="speed (ma⁻¹)");

mh = sol.vols./BM.area(gl)
p3 = histogram(mh, xlabel="mean thickness (m)", label="", bins=30,
               xlims=(minimum(mh)-1, maximum(mh+1)), ylabel="frequency",
               grid=false);
yticks!(Int[]);
annotate!(minimum(mh)+1,300, text("C", 25));


p4 = BM.plottheta(sol.thetas, sol.theta0, reuse=false, plt=corrplot2, left_margin=0mm, nsize=5,
                  #labels=["b̃₂", "fₛₗ,₂", "T", "d", "dᵥ", "σₕₘ"],
                  labels=["b̃₂", "fₛₗ,₂", "T", "d", "dᵥ", "σₕₘ", "σᵥₘ"],
#                  xticks=([-0.4:0.2:0.4;], [-0.4:0.2:0.4;], [-0.4:0.2:0.4;], [-0.4:0.2:0.4;], [-0.4:0.2:0.4;], [-0.4:0.2:0.4;], ),
                  );

# pC = plot(right_margin=0mm, left_margin=0mm, grid=:none, framestyle=:none);
# annotate!(0,0.9, text("D", 25));

# l = @layout [grid(3,1){0.4w}  a{0.001w} b{0.5w}]
# p = plot(p1,p2, p3,pC, p4, layout=l, size=(1350,900), reuse=false)


l = @layout [grid(3,1){0.48w} a{0.02w} b{0.5w}]
kwargsj = [(:tickfontsize, 13), (:legendfontsize, 13),
           (:xguidefontsize, 15), (:yguidefontsize, 15),
           (:xtickfontrotation, pi/2)]

p = plot(p1,p2, p3, plot(framestyle=:none), p4, layout=l, size=(1350,900), reuse=true; kwargsj...)

#savefig("figs/f3.png")
#p
