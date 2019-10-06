## misc plotters for inverse model results
import Interpolations
const Interp = Interpolations
using Plots,StatPlots, VAWTools; pyplot()

function getunits(xyscale)
    if xyscale==1
        return " (m)"
    elseif xyscale==1e3
        return " (km)"
    else
        return " ($xyscale m)"
    end
end

"""

Plots the thetas.

plottheta(thetas, theta0; toplot=[:all, :some, :btilde, :temp, :fsl, :pm][2], plt=[corrplot, cornerplot][2])

"""
function plottheta(thetas, theta0::Theta0; toplot=[:all, :some, :btilde, :temp, :fsl, :pm][2], nsize=9,
                   plt=[corrplot, cornerplot][2], nthin=0, fsize=(600,600), labels=nothing,
                   plotoffsets=[:btilde],
                   kwargs...)
    varnames = get_varnames(theta0)
    fields = [:btilde, :fsl, :temp]
    if nthin==0
        # figure out that # points < 50000
        nthin = max(1, size(thetas,2)÷25000)
    end
    thetas = @view thetas[:,1:nthin:end]
    if toplot==:all
        inds = 1:size(thetas,1)
    elseif toplot==:some
        if length(theta0.th0)<nsize
            inds = 1:length(theta0.th0)
        else
            inds = Int[]
            for name in keys(theta0.names)
                # on per name
                push!(inds, round(Int,mean(theta0.names[name])))
            end
        end
    elseif toplot in fields
        inds = theta0.names[toplot]
    elseif toplot==:pm
        inds = Int[]
        for name in keys(theta0.names)
            if !(name in fields)
                 append!(inds, theta0.names[name])
            end
        end
    end
    thetas_plot = thetas[inds,:] # makes a needed copy

    # figure out real value as opposed to offset:
    gb = theta0.gb
    for (ii,i) in enumerate(inds)
        fldn = Symbol.(get_varnames(theta0, pad=false, addnumbers=false))[i]
        fldn in fields || continue
        fldn in plotoffsets && continue
        fld = getfield1d(gb, fldn)
        ifld  = Interp.interpolate((reverse(gb.ele),), reverse(fld), Interp.Gridded(Interp.Linear()) )
        thetas_plot[ii,:] .+= ifld[[get_ele(theta0, i)]]
    end

    if labels==nothing
        labels = varnames[inds]
    end

    return plt(thetas_plot'; size=fsize, label=labels, kwargs...)
end

"""
    plottheta_violin(thetas, theta0::Theta0, toplot=[:btilde, :temp, :fsl][1];
                          plt=[corrplot, cornerplot][2], xvar=[:ele,:x][2], width=10, xyscale=1e3, kwargs...)

Plots the multi-parameters as violin plots.  `thetas` can also be a
tuple `(thetas, thetas_prior)` and then both are plotted
"""
function plottheta_violin(thetas, theta0::Theta0, toplot=[:btilde, :temp, :fsl][1];
                          plt=[corrplot, cornerplot][2], xvar=[:ele,:χ][2], width=10, xyscale=1e3, kwargs...)
    if isa(thetas, Tuple)
        thetas_prior = thetas[2]
        thetas = thetas[1]
        hasprior = true
    else
        hasprior = false
    end
    th = KissMCMC.summarize_run(thetas)
    varnames = get_varnames(theta0)
    inds = theta0.names[toplot]

    # figure out offset:
    gb = theta0.gb
    fld = getfield1d(gb, toplot)
    ifld  = Interp.interpolate((reverse(gb.ele),), reverse(fld), Interp.Gridded(Interp.Linear()) )

    thetas_plot = thetas[inds,:] .+ ifld[theta0.th_ele[toplot]]

    xvar1, xvar2, lab = if xvar==:ele
        (gb.ele,
         theta0.th_ele[toplot],
         "elevation (m)")
    else
        fn = piecewiselinear(gb.ele, gb.xmid)
        (gb.xmid/xyscale,
         fn.(theta0.th_ele[toplot])/xyscale,
         "χ $(getunits(xyscale))")
    end

    width = abs(maximum(diff(xvar2)))*width

    plot(xvar1, fld, label="prior expectation", xlabel=lab, ylabel="$toplot (m/a)"; kwargs...)
    if hasprior
        thetas_prior_plot = thetas_prior[inds,:] .+ ifld[theta0.th_ele[toplot]]
        violin!(xvar2, vec(thetas_prior_plot), fillcolor=:blue, fillalpha=0.15, alpha=0.3, bar_width=width,
                label="prior dist. $(string(toplot))"; kwargs...)
    end
    violin!(xvar2, vec(thetas_plot), fillcolor=:green, fillalpha=0.35, bar_width=width,
            label="posterior dist. $(string(toplot))"; kwargs...)
    scatter!(xvar2, th[:mean][inds].+ ifld[theta0.th_ele[toplot]], label="")
end


"""
    plotinv1d(sol::MCMCSol; ylims1=(), ylims2=(), ylims3=(0,100), ele_only=false, kwargs...)

Plots flow line stuff
"""
function plotinv1d(sol::MCMCSol; ylims1=(), ylims2=(), ylims3=(0,100), ele_only=false, kwargs...)
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

    p1 = plot(gb.x/1e3, gb.bands_; label="surface")#, title="Flow line")
    plot!(gb.x/1e3, gb.bands_-hs1d_max, linecolor=:red, linewidth=0,
              fillrange=gb.bands_-hs1d_min, fillcolor=:red, fillalpha=0.5, label="min/max predicted bed",
              ylabel="elevation (m)")
    plot!(gb.x/1e3, gb.bands_-hs1d, label="mean predicted bed", xlims=xl, ylims=ylims1, linecolor=:black)
    if ele_only
        return p1
    end

    p2 = plot(gb.x/1e3, hs1d_max,  linecolor=:red, linewidth=0,
              fillrange=hs1d_min, fillcolor=:red, fillalpha=0.3,
              ylabel="h (m)", label="min/max")
    plot!(gb.x/1e3, hs1d+hs1d_std,  linecolor=:blue, linewidth=0,
          fillrange=hs1d-hs1d_std, fillcolor=:blue, fillalpha=0.5,
          label="std")

    plot!(gb.x/1e3, hs1d, label="mean", xlims=xl, ylims=ylims2, linecolor=:black)

    p3 = plot(gb.x/1e3, ivs1d_max,  linecolor=:red, fillrange=ivs1d_min,
              fillcolor=:red, fillalpha=0.5, label="min/max", linewidth=0,
              ylabel="v (m/a)", xlabel="χ (km)", xlims=xl, ylims=ylims3)
    plot!(gb.x/1e3, ivs1d+ivs1d_std,  linecolor=:blue, linewidth=0,
          fillrange=ivs1d-ivs1d_std, fillcolor=:blue, fillalpha=0.5,
          label="std")
    plot!(gb.x/1e3, ivs1d, label="mean", linecolor=:black)
    im = copy(ivs1d_measured)
    im[im.<0] = NaN
    im[end]=NaN
    # scatter!(gb.xmid/1e3, im, label="v 1D measured", xlims=xl)


    plot(p1,p2,p3;layout=(3,1), size=(750,900), kwargs...)
end

"""
    plotinv1d_err(sol::MCMCSol; ylims1=(), ylims2=(), ylims3=(), xyscale=1e3, outlines=true, fsize=(1000,700), kwargs...)

Collapses the 2D error (of the expectation value of h and iv vs
observations) onto elevation bands.  Plots the minimum, maximum, mean
and std of those errors against elevation bands.

A negative value means that the model value is bigger than the observation.

- negative == over-prediction
- positive == under-prediction
"""
function plotinv1d_err(sol::MCMCSol; ylims1=(), ylims2=(), ylims3=(), xyscale=1e3, outlines=true, fsize=(1000,700), kwargs...)
    @unpack mc_hs1d, mc_ivs1d, mh2d, sh2d, miv2d, siv2d = sol
    gb = sol.theta0.gb
    gl = gb.gl

    # modeled h (expectation value):
    ih  = Interp.interpolate((gl.dem.x, gl.dem.y), mh2d, Interp.Gridded(Interp.Linear()) )

    tr = deepcopy(gl.h.vals)
    _, bandi = VAWTools.bin_traj(tr, gl.dem, gb.bands, gl.glaciermask)
    # error in 2D of h
    for i=1:length(tr)
        tr.v[i] -= ih[tr.x[i],tr.y[i]] # obs-model
    end
    # collapse back onto bands
    trmin = zeros(length(gb))*NaN
    trmax = zeros(length(gb))*NaN
    trmean = zeros(length(gb))*NaN
    trstd = zeros(length(gb))*NaN
    for (i, ibs) in enumerate(bandi)
        if length(tr.v[ibs])==0
            continue
        end
        trmin[i] = minimum(tr.v[ibs])
        trmax[i] = maximum(tr.v[ibs])
        trmean[i] = mean(tr.v[ibs])
        trstd[i] = std(tr.v[ibs])
    end

    # IV
    iv = gl.iv.vals
    iiv  = Interp.interpolate((gl.dem.x, gl.dem.y), miv2d, Interp.Gridded(Interp.Linear()) )
    iv_err = similar(iv.v)*NaN
    for iy=1:length(iv.y)
        y = iv.y[iy]
        for ix=1:length(iv.x)
            x = iv.x[ix]
            iv_ = iv.v[ix,iy] # obs
            m = gl.ivmask[ix,iy]
            if m && iv_>=0
                iv_err[ix,iy] = iv_ - iiv[x,y] # obs-model
            end
        end
    end
    ive = Gridded(iv.x, iv.y, iv_err)
    # Here it's possible that there are no points in a band because the iv resolution
    # is lower:
    bandi = VAWTools.bandi_for_other_grid(gb.bands, gb.bandi, gb.dem, ive, gl.ivmask)
    # onto bands
    ivmin = zeros(length(gb))*NaN
    ivmax = zeros(length(gb))*NaN
    ivmean = zeros(length(gb))*NaN
    ivstd = zeros(length(gb))*NaN
    for (i, ibs) in enumerate(bandi)
        if length(ive.v[ibs])==0
            continue
        end
        ivmin[i] = minimum(ive.v[ibs])
        ivmax[i] = maximum(ive.v[ibs])
        ivmean[i] = VAWTools.meannan(ive.v[ibs])
        ivstd[i] = VAWTools.stdnan(ive.v[ibs])
    end

    # Plots
    xl = (gb.x[1]/1e3,gb.x[end]/1e3)
    p1 = plotinv1d(sol, ele_only=true, xlims=xl, ylims1=ylims1)

    # h
    p2 = plot(gb.xmid/1e3, trmax,  linecolor=:red, linewidth=0,
              fillrange=trmin, fillcolor=:red, fillalpha=0.3,
              ylabel="error h  (m)", label="min/max", xlims=xl)
    plot!(gb.xmid/1e3, trmean+trstd,  linecolor=:blue, linewidth=0,
          fillrange=trmean-trstd, fillcolor=:blue, fillalpha=0.5,
          label="std", xlims=xl)

    plot!(gb.xmid/1e3, trmean, lw=1, color=:black, label="mean", ylims=ylims2,
          legend=:topleft)

    # IV
    p3 = plot(gb.xmid/1e3, ivmax, linecolor=:red, linewidth=0,
              fillrange=ivmin, fillcolor=:red, fillalpha=0.3,
              ylabel="error v (m/a)", label="min/max")
    plot!(gb.xmid/1e3, ivmean+ivstd,  linecolor=:blue, linewidth=0,
          fillrange=ivmean-ivstd, fillcolor=:blue, fillalpha=0.5,
          label="std", xlims=xl)

    plot!(gb.xmid/1e3, ivmean, marker=:x, lw=1, color=:black, xlabel="χ (km)",
          label="mean", ylims=ylims3,
          legend=:topleft)

    plot(p1,p2,p3,layout=(3,1), size=(750,900); kwargs...)
end

"""
    plotinv2d(sol::MCMCSol; xyscale=1e3, outlines=true, logiv=false,
                   clim1h=(), clim2h=(), clim3h=(), clim4h=(),
                   clim1iv=(), clim2iv=(), clim3iv=(), clim4iv=(),
                   kwargs...)

Plots 2D stuff:
- mean h & iv
- min h & iv
- max h & iv
- std h & iv
"""
function plotinv2d(sol::MCMCSol; xyscale=1e3, outlines=true, logiv=false,
                   clim1h=(), clim2h=(), clim3h=(), clim4h=(),
                   clim1iv=(), clim2iv=(), clim3iv=(), clim4iv=(),
                   kwargs...)
    @unpack mc_hs1d, mc_ivs1d, mh2d, sh2d, miv2d, siv2d,
    min_h2d, max_h2d, min_iv2d, max_iv2d = sol
    gb = sol.theta0.gb
    gl = gb.gl

    out = []
    for f2d in [mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d]
        tmp = copy(f2d)
        f2d[.!gl.glaciermask] = NaN
        push!(out, f2d)
    end
    mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = out

    # mean
    h1 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, mh2d', colorbar_title="mean(h) (m)", fill=true,
                 aspect_ratio=:equal, clims=clim1h)
    outlines && plot_outlines!(gl, xyscale=xyscale)

    # min
    h2 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, min_h2d', colorbar_title="minimum(h) (m)", fill=true,
                 aspect_ratio=:equal, clims=clim2h)
    outlines && plot_outlines!(gl, xyscale=xyscale)

    # max
    h3 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, max_h2d', colorbar_title="maximum(h) (m)", fill=true,
                 aspect_ratio=:equal, clims=clim3h)
    outlines && plot_outlines!(gl, xyscale=xyscale)

    # std
    h4 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, sh2d', colorbar_title="std(h) (m)", fill=true,
                 aspect_ratio=:equal, clims=clim4h)
    outlines && plot_outlines!(gl, xyscale=xyscale)


    fn = logiv ? x->log10(x) : x->x
    st = logiv ? "log" : ""
    # mean
    i1 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, fn(miv2d'), colorbar_title="$st mean(iv) (m/a)", fill=true,
                 aspect_ratio=:equal, clims=clim1iv)
    outlines && plot_outlines!(gl, xyscale=xyscale)
    # min
    # xref https://github.com/tbreloff/Plots.jl/issues/533
    min_iv2d = copy(min_iv2d)
    if length(clim2iv)==2
        min_iv2d[min_iv2d.<clim2iv[1]] = clim2iv[1]
        min_iv2d[min_iv2d.>clim2iv[2]] = clim2iv[2]
    end
    i2 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, fn(min_iv2d'), colorbar_title="$st minimum(iv) (m/a)", fill=true,
                 aspect_ratio=:equal, clims=clim2iv)
    outlines && plot_outlines!(gl, xyscale=xyscale)
    # max
    max_iv2d = copy(max_iv2d)
    if length(clim3iv)==2
        max_iv2d[max_iv2d.<clim3iv[1]] = clim3iv[1]
        max_iv2d[max_iv2d.>clim3iv[2]] = clim3iv[2]
    end
    i3 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, fn(max_iv2d'), colorbar_title="$st maximum(iv) (m/a)", fill=true,
                 aspect_ratio=:equal, clims=clim3iv)
    outlines && plot_outlines!(gl, xyscale=xyscale)

    # std
    siv2d = copy(siv2d)
    if length(clim4iv)==2
        siv2d[siv2d.<clim4iv[1]] = clim4iv[1]
        siv2d[siv2d.>clim4iv[2]] = clim4iv[2]
    end
    i4 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, fn(siv2d'), colorbar_title="$st std(iv) (m/a)", fill=true,
                 aspect_ratio=:equal, clims=clim4iv)
    outlines && plot_outlines!(gl, xyscale=xyscale)

    plot(h1, i1, h2, i2, h3, i3, h4, i4, layout=(4,2); kwargs...)
end

"1D and 2D"
function plotinv(sol::MCMCSol; kwargs...)
    p1 = plotinv1d(sol; kwargs...)

    p2 = plotinv2d(sol; kwargs...)
    plot(p1,p2, layout=(2,1), size=(700,1200); kwargs...)
end


"""
Plots modeled thickness and thickness along radar lines.
"""
function plotinv2d_h(sol::MCMCSol; toplot=[:mean,:max,:min][1], plotiv=false, kwargs...)
    @unpack mc_hs1d, mc_ivs1d, mh2d, sh2d, miv2d, siv2d, min_h2d,
    max_h2d, min_iv2d, max_iv2d = sol
    gb = sol.theta0.gb
    gl = gb.gl
    if toplot==:mean
        top = plotiv ? miv2d : mh2d
    elseif toplot==:min
        top = plotiv ? min_iv2d : min_h2d
    elseif toplot==:max
        top = plotiv ? max_iv2d : max_h2d
    end
    plot2d_h(gl, top; plotiv=plotiv, kwargs...)
end

"""
Plots difference between measured iv and modeled
"""
function plotinv2d_iverr(sol::MCMCSol; xyscale=1e3, outlines=true, fsize=(1000,700),
                         toplot=[:mean,:max,:min][1], kwargs...)
    @unpack mc_hs1d, mc_ivs1d, mh2d, sh2d, miv2d, siv2d,
    min_h2d, max_h2d, min_iv2d, max_iv2d = sol
    gb = sol.theta0.gb
    gl = gb.gl

    if toplot==:mean
        top = miv2d
    elseif toplot==:min
        top = min_iv2d
    elseif toplot==:max
        top = max_iv2d
    end

    iiv  = Interp.interpolate((gl.dem.x, gl.dem.y), top, Interp.Gridded(Interp.Linear()) )
    iv_err = similar(gl.iv.vals.v)*NaN
    for iy=1:length(gl.iv.vals.y)
        y = gl.iv.vals.y[iy]
        for ix=1:length(gl.iv.vals.x)
            x = gl.iv.vals.x[ix]
            iv = gl.iv.vals.v[ix,iy]
            m = gl.ivmask[ix,iy]
            if m && iv>=0
                iv_err[ix,iy] = iv - iiv[x,y]
            end
        end
    end
    if !(length(iv_err)>1)
        warn("Cannot plot this when length(gl.iv.vals.v)==1.")
        h1 = plot()
        outlines && plot_outlines!(gl, xyscale=xyscale)
        return h1
    end
    h1 = contour(gl.iv.vals.x/xyscale, gl.iv.vals.y/xyscale, iv_err', title="", fill=true,
                 aspect_ratio=:equal, colorbar_title="v error (m/a) (negative: model overestimates)"; kwargs...)
    outlines && plot_outlines!(gl, xyscale=xyscale)
    h1
end
