# sol, sol_prior = open("figs/f2-10.jls", "r") do io
#     deserialize(io)
# end

# # ### plots

# if plotyes
#     display(BM.plotinv1d(sol, reuse=false))
#     display(BM.plotinv1d_err(sol, reuse=false))

#     display(BM.plotinv2d(sol, reuse=false))
#     display(BM.plotinv2d_h(sol, reuse=false))
#     display(BM.plotinv2d_iverr(sol, reuse=false))
#     @unpack thetas, theta0 = sol
#     display(BM.plottheta(thetas, theta0))
#     display(BM.plottheta(thetas, theta0, toplot=:btilde, reuse=false))
#     display(BM.plottheta(thetas, theta0, toplot=:fsl, reuse=false))
#     display(BM.plottheta(thetas, theta0, toplot=:temp, reuse=false))

#     display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :btilde, reuse=false, width=1))
#     display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :fsl, reuse=false, width=1))
#     display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :temp, reuse=false, width=1))
# end

using Plots, Measures, Parameters
function plot_run_custom(sol;
                         outlines=true, reuse=false, plot2d=[:h,:iv,:ivlog][1], xyscale=1e3, possibly_flip=true,
                         clims=[:model,:radar][1], radarline_gl=gb.gl, plot1d=[:smb, :flux][1], kwargs...)
    kwargs = [(:lw, 2)]
    # make forward solution
    @unpack theta0 = sol
    theta = sol.theta_expect
    @unpack gb = theta0
    @unpack pp,pn,pm = theta0
    btilde, temp, fsl = BM.make_fields1d(theta, theta0)

    # pre-calculate an update structure for pm
    args_to_MPara = []
    inds = Tuple{Int,Int}[]
    for (ii,name) in enumerate(fieldnames(pm))
        if haskey(theta0.names,name)
            its = theta0.names[name]
            push!(inds, (ii, its[1]))
            push!(args_to_MPara, nothing)
        else
            push!(args_to_MPara, getfield(pm,name))
        end
    end
    for (ias, its) in inds
        args_to_MPara[ias] = theta[its]
    end

    pm_ = BM.MPara(args_to_MPara...)

    fwdsol = BM.fwdm(gb, pp, pm_, pn, btilde, fsl, temp)

    # plot it
    @unpack (gb, hs1d, taus1d, ivs1d, hs2d, taus2d,
             ivs2d_band, taus1d_local, pp, pm, btilde, fsl) = fwdsol

    gl = gb.gl
    qtot1d, qd1d = BM.postproc1d(gb, pp, pm, btilde, fsl)

    p1 = plot(gb.x/xyscale, gb.bands_, label="surface", title="", right_margin=20.0mm; kwargs...)
    plot!(gb.x/xyscale, gb.bands_-hs1d, label="bed",
          xlabel="", ylabel="elevation (m a.s.l.)")
    # make a non-showing line to get the legend for the width
    plot!(gb.x/xyscale, gb.bands_.*NaN, label="width", linecolor=:green, linestyle=:dash; kwargs...)
    annotate!(10, 2800, text("A", 30))

    plot!(twinx(p1), gb.x/xyscale, gb.ws, label="width", linecolor=:green, linestyle=:dash,
          xlabel="", ylabel="width (m)", legend=false, grid=false, framestyle=:box; kwargs...)


    p2 = plot(gb.x/xyscale, taus1d/1e3,
              ylabel="driving stress (kPa)", xlabel="χ (km)",
              label="driving stress averaged"; kwargs...)
    plot!(gb.xmid/xyscale, taus1d_local/xyscale,
          ylabel="driving stress (kPa)", xlabel="",
          label="driving stress local", framestyle=:box; kwargs...)
    plot!(gb.xmid/xyscale, taus1d_local.*NaN,
          linecolor=:green, linestyle=:dash,
          label="sliding fraction"; kwargs...)
    annotate!(10, 150, text("B", 30))

    plot!(Plots.twinx(), gb.xmid/xyscale, fsl,
          xlabel="", ylabel="sliding fraction ()", legend=false,
          linecolor=:green, linestyle=:dash,
          label="", framestyle=:box, grid=false; kwargs...)

    # IV:
    # p3 = plot(gb.xmid/xyscale, ivs1d,
    #           ylabel="iv (m/a)", xlabel="χ (km)",
    #           label="")
    if plot1d==:smb
        # p3 = plot(gb.xmid/xyscale, gb.bdot,
        #           xlabel="χ (km)", ylabel="mass (m /a)",
        #           label="ḃ"; kwargs...)
        # p3 = plot!(gb.xmid/xyscale, gb.bdot - btilde,
        #            xlabel="χ (km)",
        #            label="dh/dt (consistent)"; kwargs...)
        # p3 = plot!(gb.xmid/xyscale, gb.dhdt,
        #            xlabel="χ (km)",
        #            label="dh/dt (observation)"; kwargs...)
        # p3 = plot!(gb.xmid/xyscale, gb.dhdt,
        #            xlabel="χ (km)",
        #            label="dh/dt as in gb (m ice /a)")
        p3 = plot(gb.xmid/xyscale, gb.bdot - gb.dhdt,
                  xlabel="χ (km)", ylabel="b̃ (ma⁻¹)",
                  label="b̃ obs."; kwargs...)
        p3 = plot!(gb.xmid/xyscale, btilde,
                   xlabel="χ (km)",
                   label="b̃ fitted"; kwargs...)
        p3 = plot!(gb.x/xyscale, qtot1d.*NaN,
                   label="q",linecolor=:green,
                   linestyle=:dash; kwargs...)
        p3 = plot!(gb.x/xyscale, qd1d.*NaN,
                   label="q_d", linecolor=:cyan,
                   linestyle=:dash; kwargs...)
        annotate!(10, 0.6, text("C", 30))
        p33 = Plots.twinx()
        plot!(p33, gb.x/xyscale, qtot1d,
              ylabel="flux (m³s⁻¹)", xlabel="χ (km)",
                   label="",linecolor=:green,
                   linestyle=:dash; kwargs...)
        plot!(p33, gb.x/xyscale, qd1d,
                   label="", linecolor=:cyan,
              linestyle=:dash, grid=false,
              framestyle=:box; kwargs...)

    elseif plot1d==:flux
        p3 = plot(gb.x/xyscale, qtot1d,
                   ylabel="q (m³s⁻¹)", xlabel="χ (km)",
                   label="q"; kwargs...)
        p3 = plot!(gb.x/xyscale, qd1d,
                   ylabel="", xlabel="χ (km)",
                   label="q_d"; kwargs...)

        p3 = plot!(gb.x/xyscale, qd1d.*NaN,
                   label="ice speed (m/a)", linecolor=:green, linestyle=:dash; kwargs...)


        p3 = plot!(Plots.twinx(), gb.x/xyscale, ivs1d,
                   ylabel="iv (m/a)", xlabel="χ (km)", linecolor=:green,
                   linestyle=:dash,
                   label="", grid=false, framestyle=:box; kwargs...)
        # p3 = plot!(gb.xmid/xyscale, ivs1d/100,
        #            ylabel="", xlabel="χ (km)",
        #            label="iv (100 m/a)")
    end


    #println("hmax=$(maximum(hs1d)), hmean = $(mean(hs1d))")
    # warn("tau and iv super noisy?")
    ivs2d = copy(ivs2d_band) #;
    ivs2d[ivs2d.==0] = NaN
    flip = if possibly_flip
        abs(gl.dem.x[end]-gl.dem.x[1]) > 1.2* abs(gl.dem.y[end]-gl.dem.y[1])
    else
        false
    end

    kwargs = [ (:aspect_ratio,:equal),
        (:framestyle,:none), (:grid,:none), (:xlabel,""), (:ylabel,""),
               (:ylims, (5153.8, 5162.2)) ]
    p4 = BM.plot2d_h(radarline_gl, hs2d; clims=clims, xyscale=xyscale,
                     outlines=false, fsize=(1000,700), plotiv=false, plotline=true, xlims=(432.2,442.7),
                     flip=flip, left_margin=4mm, colorbar_title="thickness (m)", right_margin=4mm, kwargs...)
    annotate!(441.5, 5159.5, text("D", 30))

    p5 = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, (ivs2d'); size=(600,600), xlims=(432.2,442.7),
                         aspect_ratio=:equal, colorbar_title="ice flow speed (ma⁻¹)", left_margin=4mm, kwargs...)
    annotate!(441.5, 5159.5, text("E", 30))

    outlines && BM.plot_outlines!(gl, xyscale=xyscale, flip=flip, lw=1.5)
    # scale bar
    plot!([440], [5155], xerror=[-1,1], markercolor=:black, markerstrokecolor=:black, markerstrokewidth=2, markersize=10, label="")
    annotate!(440, 5155.2, text("2 km", 15))

    plot(p5;kwargs...)


    #p4 = contour(gl.dem.x, gl.dem.y, taus2d', fill=true, size=(600,600))

    l = @layout [grid(3,1){0.5w}  grid(2,1){0.5w}] #a{0.4w}]
    # l = @layout [grid(3,1) [a{0.4w}
    #                         b{0.4w}]]

    kwargs = [(:tickfontsize, 13), (:labelfontsize, 35), (:legendfontsize, 10),
              (:xguidefontsize, 15), (:yguidefontsize, 15)]

    plot(p1,p2,p3, p4,p5, layout=l, reuse=reuse; size=(1350,900), kwargs...)
end
# Figure 2
p = plot_run_custom(sol, possibly_flip=false)
savefig("figs/f2.png")
savefig("figs/f2.pdf")
