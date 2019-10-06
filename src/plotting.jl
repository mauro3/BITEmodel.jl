# `include` this file if you want to do plotting.
using Plots; pyplot()

# function __init__()
#     pyplot() # https://github.com/JuliaPlots/Plots.jl/issues/745
#     #gr()
#     default(reuse=false)
# end
using StatPlots

## tools.jl
############
# Plotting
function Plots.plot(g::Gridded; xyscale=1e3, kwargs...)
    v = copy(g.v)
    v[v.==FILL] = NaN
    Plots.contourf(g.x/xyscale, g.y/xyscale, v';
                   xlabel="Easting",
                   ylabel="Northing",
                   #colorbar=:none,
                   kwargs...)
end
function Plots.plot!(g::Gridded; xyscale=1e3, kwargs...)
    v = copy(g.v)
    v[v.==FILL] = NaN
    Plots.contourf!(g.x/xyscale, g.y/xyscale, v';
                    xlabel="Easting",
                    ylabel="Northing",
                    #colorbar=:none,
                    kwargs...)
end
function Plots.plot(tr::Traj; xyscale=1e3, ploterr=false, kwargs...)
    Plots.plot(tr.x/xyscale, tr.y/xyscale, tr.v;
              xlabel="Easting",
               ylabel="Northing",
               label="", kwargs...)
    if ploterr
        Plots.plot!(tr.x/xyscale, tr.y/xyscale, tr.err, label="")
    end
end

function VAWTools.plot_bands(gb::Bands; outline=true, kwargs...)
    VAWTools.plot_bands(gb.dem,gb.bands,gb.bandi; kwargs...)
    outline && plot_outlines!(gb.gl, xyscale=1)
end



##################
## data-structs.jl
##################
function Plots.plot(gl::Glacier; xyscale=1e3, field=gl.dem, kwargs...)
    plot(field; xyscale=xyscale, kwargs...)
    plot_outlines!(gl; xyscale=xyscale, kwargs...)
end


## Glaciers
###########

# seriestypes:
#[:none,:line,:path,:steppre,:steppost,:sticks,:scatter,:heatmap,:hexbin,:histogram,:histogram2d,:histogram3d,:density,:bar,:hline,:vline,:contour,:pie,:shape,:image,:path3d,:scatter3d,:surface,:wireframe,:contour3d]

"""
$(SIGNATURES)

kwargs:

     outlines=true, reuse=false, plot2d=[:h,:iv,:ivlog][1], xyscale=1e3, possibly_flip=true,
     clims=[:model,:radar][1], radarline_gl=gb.gl, plot1d=[:smb, :flux][1], kwargs...

Makes summary plot of results

Note, to zoom with updating ticks use `ticks=:native`.
"""
plot_run(fwdsol::FWDSol; kwargs...) = plot_run(fwdsol.gb, fwdsol.hs1d, fwdsol.taus1d,
                                               fwdsol.ivs1d, fwdsol.hs2d, fwdsol.taus2d,
                                               fwdsol.ivs2d_band, fwdsol.taus1d_local, fwdsol.pp,
                                               fwdsol.pm, fwdsol.btilde, fwdsol.fsl;
                                               kwargs...)
function plot_run(sol::MCMCSol; theta = sol.theta_mode, kwargs...)
    @unpack theta0 = sol
    @unpack gb = theta0
    @unpack pp,pn,pm = theta0
    btilde, temp, fsl = make_fields1d(theta, theta0)

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

    pm_ = MPara(args_to_MPara...)

    fwdsol = fwdm(gb, pp, pm_, pn, btilde, fsl, temp)

    plot_run(fwdsol; kwargs...)
end
function plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, pp, pm, btilde, fsl;
                  outlines=true, reuse=false, plot2d=[:h,:iv,:ivlog][1], xyscale=1e3, possibly_flip=true,
                  clims=[:model,:radar][1], radarline_gl=gb.gl, plot1d=[:smb, :flux][1], kwargs...)
    gl = gb.gl
    qtot1d, qd1d = postproc1d(gb, pp, pm, btilde, fsl)

    p1 = plot(gb.x/xyscale, gb.bands_, label="surface", title="$(getname(gb.gl)) $(gl.gid.id)")
    plot!(gb.x/xyscale, gb.bands_-hs1d, label="bed",
          xlabel="χ (km)", ylabel="z (m)")
    # make a non-showing line to get the legend for the width
    plot!(gb.x/xyscale, gb.bands_.*NaN, label="width", linecolor=:green, linestyle=:dash)

    plot!(Plots.twinx(p1), gb.xmid/xyscale, gb.ws, label="width", linecolor=:green, linestyle=:dash,
          xlabel="χ (km)", ylabel="width (m)", legend=false, grid=false, framestyle=:box)


    p2 = plot(gb.x/xyscale, taus1d/1e3,
              ylabel="τ (kPa)", xlabel="χ (km)",
              label="τ averaged")
    plot!(gb.xmid/xyscale, taus1d_l/xyscale,
          ylabel="τ (kPa)", xlabel="χ (km)",
          label="τ local", framestyle=:box)
    plot!(gb.xmid/xyscale, taus1d_l.*NaN,
          linecolor=:green, linestyle=:dash,
          label="fsl")

    plot!(Plots.twinx(), gb.xmid/xyscale, fsl,
          ylabel="fsl", xlabel="χ (km)", legend=false,
          linecolor=:green, linestyle=:dash,
          label="", framestyle=:box, grid=false)

    # IV:
    # p3 = plot(gb.xmid/xyscale, ivs1d,
    #           ylabel="iv (m/a)", xlabel="χ (km)",
    #           label="")
    if plot1d==:smb
        p3 = plot(gb.xmid/xyscale, gb.bdot,
                  ylabel="", xlabel="χ (km)",
                  label="bdot GloGEM")
        p3 = plot!(gb.xmid/xyscale, -(gb.bdot - btilde),
                   ylabel="", xlabel="χ (km)",
                   label="-dh/dt consistent")
        p3 = plot!(gb.xmid/xyscale, -gb.dhdt,
                   ylabel="", xlabel="χ (km)", linestyle=:dot,
                   label="-dh/dt data")
        p3 = plot!(gb.xmid/xyscale, btilde,
                   ylabel="MB (m ice /a)", xlabel="χ (km)",
                   label="b̃")
        p3 = plot!(gb.x/xyscale, qtot1d.*NaN,
                   label="q",linecolor=:green,
                   linestyle=:dash)
        p3 = plot!(gb.x/xyscale, qd1d.*NaN,
                   label="q_d", linecolor=:cyan,
                   linestyle=:dash)

        p33 = Plots.twinx()
        plot!(p33, gb.x/xyscale, qtot1d,
                   ylabel="q (m^3/s)", xlabel="χ (km)",
                   label="",linecolor=:green,
                   linestyle=:dash)
        plot!(p33, gb.x/xyscale, qd1d,
                   label="", linecolor=:cyan,
              linestyle=:dash, grid=false,
              framestyle=:box)

    elseif plot1d==:flux
        p3 = plot(gb.x/xyscale, qtot1d,
                   ylabel="q (m^3/s)", xlabel="χ (km)",
                   label="q")
        p3 = plot!(gb.x/xyscale, qd1d,
                   ylabel="", xlabel="χ (km)",
                   label="q_d")

        p3 = plot!(gb.x/xyscale, qd1d.*NaN,
                   label="iv (m/a)", linecolor=:green, linestyle=:dash)


        p3 = plot!(Plots.twinx(), gb.x/xyscale, ivs1d,
                   ylabel="iv (m/a)", xlabel="χ (km)", linecolor=:green,
                   linestyle=:dash,
                   label="", grid=false, framestyle=:box)
        # p3 = plot!(gb.xmid/xyscale, ivs1d/100,
        #            ylabel="", xlabel="χ (km)",
        #            label="iv (100 m/a)")
    end


    #println("hmax=$(maximum(hs1d)), hmean = $(mean(hs1d))")
    # warn("tau and iv super noisy?")
    ivs2d = copy(ivs2d) #;
    ivs2d[ivs2d.==0] = NaN
    flip = if possibly_flip
        abs(gl.dem.x[end]-gl.dem.x[1]) > 1.2* abs(gl.dem.y[end]-gl.dem.y[1])
    else
        false
    end


    if plot2d==:h
        p4 = plot2d_h(radarline_gl, hs2d; clims=clims, xyscale=xyscale,
                      outlines=false, fsize=(1000,700), plotiv=false, plotline=true,
                      flip=flip, kwargs...)

        # p4 = contour(gl.dem.x/xyscale, gl.dem.y/xyscale, hs2d2', fill=true, size=(600,600),
        #              aspect_ratio=:equal)
    elseif plot2d==:iv
        if flip
            p4 = heatmap(gl.dem.y/xyscale, gl.dem.x/xyscale, (ivs2d), size=(600,600), #, zlims=(0,50))
                         aspect_ratio=:equal, colorbar_title="iv (m/a)", xflip=true, kwargs...)
        else
            p4 = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, (ivs2d'), size=(600,600), #, zlims=(0,50))
                         aspect_ratio=:equal, colorbar_title="iv (m/a)", kwargs...)
        end
    elseif plot2d==:ivlog
        if flip
            p4 = heatmap(gl.dem.y/xyscale, gl.dem.x/xyscale, log10.(ivs2d), size=(600,600),
                         aspect_ratio=:equal, colorbar_title="log(iv) (m)", xflip=true, kwargs...)
        else
            p4 = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, log10.(ivs2d'), size=(600,600),
                         aspect_ratio=:equal, colorbar_title="log(iv) (m)", kwargs...)
        end
    end

    outlines && plot_outlines!(gl, xyscale=xyscale, flip=flip)


    #p4 = contour(gl.dem.x, gl.dem.y, taus2d', fill=true, size=(600,600))

    l = @layout [grid(3,1) a{0.4w}]
    plot(p1,p2,p3,p4, layout=l, reuse=reuse; size=(900,600), kwargs...)
end

function plot_outlines!(gl; xyscale=10^3, plot_outcrops=true, flip=false, label="", kwargs...)
    units = getunits(xyscale)
    xind,yind = flip ? (2,1) : (1,2)
    xl,yl = flip ? ("y $units","x $units") : ("x $units", "y $units")
    # outlines
    splits = gl.outline.splits
    p4 = 1
    inds = splits[1]:splits[2]-1
    p4 = plot!(gl.outline.outline[xind,inds]/xyscale,gl.outline.outline[yind,inds]/xyscale; label=label, linecolor=:green,
               xlabel=xl, ylabel=yl, xflip=flip, aspect_ratio=:equal, kwargs...)
    for i=2:length(splits)-1
        inds = splits[i]:splits[i+1]-2
        p4 = plot!(gl.outline.outline[xind,inds]/xyscale,gl.outline.outline[yind,inds]/xyscale; label="", linecolor=:green,
                   kwargs...)
    end
    if plot_outcrops
        ocs = split_poly(gl.outline.outcrops, gl.outline.outcrops_splits)
        for oc in ocs
            p4 = plot!(oc[xind,:]/xyscale, oc[yind,:]/xyscale; label="", linecolor=:red, kwargs...)
        end
    end
    return p4
    #ocs = split_poly(gl.outcrops[1], gl.outcrops[2])

    #sea: messes up the colorbar
    #contour!(gl.dem.x, gl.dem.y, gl.dem.v', levels=0:1, c=:reds, linewidth=2)
end

"""
   plot_on_bands(gb::Bands, var, cname="Blues", ncolors=100;
                 size=(1000,1000), xyscale=1e3, kwargs...)

Plot a band variable extrapolated to 2D using the simplest extrapolation.

Only works with plotly.
"""
function plot_on_bands(gb::Bands, var, cname="Blues", ncolors=100;
                       size=(1000,1000), xyscale=1e3, kwargs...)
    gl = gb.gl

    nb = length(gb)

    # rng = linspace(minimum(var), maximum(var), ncolors)
    # colors = colormap(cname, ncolors)
    # cols = Array{eltype(colors)}(nb)

    # map onto 2D
    out = similar(gl.dem.v)*NaN
    for (i,bi) in enumerate(gb.bandi)
        out[bi] = var[i]
    end
    heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, out'; clim=(minimum(var),maximum(var)), size=size, kwargs...)
end

"""
    plot_radar(gl; xyscale=1e3)

Plots surface and bed along radar lines as trajectories in 3D.
Note that traj.err needs to hold surface elevation.
"""
function plot_radar(gl; xyscale=1e3)
    tr = gl.h.vals
    kargs = Dict(:label=>"") # slow, :linealpha=>0, :marker=>:o)
    p = Plots.plot(xlabel="Easting",
                   ylabel="Northing",
                   lc=:blue; kargs...)
    #    Plots.plot!(tr.x/xyscale, tr.y/xyscale, tr.err-tr.v, lc=:red; kargs...)
    for inds in gl.h.vals.splits
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.err[inds], lc=:blue; kargs...)
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.err[inds]-tr.v[inds], lc=:red , markershape=:x, markercolor=:red ; kargs...)
    end
    return p
end

"""
As plot_radar but plots just thickness as 3D trajectory.
"""
function plot_radar_thick(gl; xyscale=1e3)
    tr = gl.h.vals
    kargs = Dict(:label=>"") # slow, :linealpha=>0, :marker=>:o)
    p = Plots.plot(xlabel="Easting",
                   ylabel="Northing",
                   lc=:blue; kargs...)
#    Plots.plot!(tr.x/xyscale, tr.y/xyscale, tr.err-tr.v, lc=:red; kargs...)
    for inds in gl.h.vals.splits
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.v[inds], lc=:red , markershape=:x, markercolor=:red ; kargs...)
    end
    return p
end

"""
    plot_radar_h(gl, hs2d; xyscale=1e3, plot_floation=false, pp=Phys())

Plot radar and 2D thickness field.
"""
function plot_radar_h(gl, hs2d; xyscale=1e3, plot_floation=false, pp=Phys())
    itp = Interp.interpolate((gl.dem.x, gl.dem.y), hs2d, Interp.Gridded(Interp.Linear()) )
    itps = Interp.interpolate((gl.dem.x, gl.dem.y), gl.dem.v, Interp.Gridded(Interp.Linear()) )
    tr = gl.h.vals
    # setup plot
    kargs = Dict{Symbol,Any}()
    p = Plots.plot(tr.x/xyscale, tr.y/xyscale, tr.err,
                   xlabel="Easting",
                   ylabel="Northing",
                   label="surf",
                   lc=:blue; kargs...)
    for inds in gl.h.vals.splits
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.err[inds], lc=:blue, label=""; kargs...)
        s = [itps[x,y] for (x,y) in zip(tr.x[inds],tr.y[inds])]
        s[s.<0] = NaN
        s[s.>1000] = NaN
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, s, lc=:purple,
                    label="surf"; kargs...)
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.err[inds]-tr.v[inds], lc=:red,
                    label="radar bed"; kargs...)
        Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, tr.err[inds]-[itp[x,y] for (x,y) in zip(tr.x[inds],tr.y[inds])], lc=:green,
                    label="modeled bed"; kargs...)
        if plot_floation
            Plots.plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale,
                        s-s/((pp.rhosw-pp.rho)/pp.rhosw),
                        lc=:cyan, label="flotation"; kargs...)
        end
        kargs = Dict(:label=>"") # slow, :linealpha=>0, :marker=>:o)
    end
    return p
end

# function plot_radar_h_err(gl, hs2d; xyscale=1e3, plot_floation=false,
#                           pp=Phys(), fsize=(1000,700))
#     itp = Interp.interpolate((gl.dem.x, gl.dem.y), hs2d, Interp.Gridded(Interp.Linear()) )
#     itps = Interp.interpolate((gl.dem.x, gl.dem.y), gl.dem.v, Interp.Gridded(Interp.Linear()) )
#     tr = gl.h[1]
#     max_h = 0.0
#     min_h = 0.0
#     for tr in gl.h
#         h = [itp[x,y] for (x,y) in zip(tr.x,tr.y)]
#         max_h = max(max_h,maximum(tr.v-h))
#         min_h = min(min_h,minimum(tr.v-h))
#     end
#     clims = (min_h,max_h)

#     # setup plot
#     kargs = Dict{Symbol,Any}()
#     h1 = Plots.plot(xlabel="Easting",
#                    ylabel="Northing",
#                    label="surf",
#                    lc=:blue; kargs...)
#     for tr in gl.h
#         h = [itp[x,y] for (x,y) in zip(tr.x,tr.y)]
#         scatter!(tr.x/xyscale, tr.y/xyscale, zcolor=tr.v-h, markerstrokewidth=0, label="", colorbar=false,
#                  clims=clims)
#         # Plots.plot!(tr.x/xyscale, tr.y/xyscale, tr.v-h, lc=:purple,
#         #             label="error h"; kargs...)
#         kargs = Dict(:label=>"") # slow, :linealpha=>0, :marker=>:o)
#     end
#     plot_outlines!(gl, xyscale=xyscale)

#     # fix https://github.com/tbreloff/Plots.jl/issues/412
#     h2 = scatter([0,0], [0,1], zcolor=[clims...], clims=clims,# aspect_ratio=:equal,
#                  xlims=(1,1.1), axis=false, label="", colorbar_title="ice thickness (m)")#, markeralpha=0, markerstrokewidth=0)
#     l = @layout [grid(1,1) a{0.01w}]
#     p=plot(h1,h2,layout=l, palette=:viridis)


#     return p
# end

"""
    plot_traj!(tr; xyscale=1e3, toplot=[:v,:err][1], plotline=true, kwargs...)

Plot a trajectory in 2D with its value encoded in the color.
"""
function plot_traj!(tr; xyscale=1e3, toplot=[:v,:err][1], plotline=true, kwargs...)
    # p = Plots.plot(xlabel="Easting",
    #                ylabel="Northing",
    #                lc=:blue; kwargs...)
    p = 0
    top = toplot==:v ? tr.v : tr.err
    for (i,inds) in enumerate(tr.splits)
        p = Plots.scatter!(tr.x[inds]/xyscale, tr.y[inds]/xyscale; label="", zcolor=top[inds],
                           markerstrokecolor=:white, markerstrokewidth=0, kwargs...)
        if plotline==true
            plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, color=:white, label="", lw=0.5)
        elseif plotline isa AbstractVector
            if i in plotline
                # only plot if line in in plotline
                plot!(tr.x[inds]/xyscale, tr.y[inds]/xyscale, color=:white, label="", lw=1.0)
            end
        elseif plotline==false
        else
            error()
        end
    end
    return p
end


"Radar plot in 2d"
function plot_radar2d(gl; xyscale=1e3, outlines=true, fsize=(1000,700))
    max_h = maximum(gl.h.vals.v)
    clims = (0,max_h)

    h1 = plot(aspect_ratio=:equal)
    plot_traj!(gl.h.vals, label="", colorbar=false, clims=clims)

    outlines && plot_outlines!(gl, xyscale=xyscale, plot_outcrops=false)
    h1
    # fix https://github.com/tbreloff/Plots.jl/issues/412
    h2 = scatter([0,0], [0,1], zcolor=[clims...], clims=clims,# aspect_ratio=:equal,
                 xlims=(1,1.1), axis=false, label="", colorbar_title="h (m)")#, markeralpha=0, markerstrokealpha=0)
    l = @layout [grid(1,1) a{0.01w}]
    plot(h1,h2,layout=l, size=fsize, palette=:viridis)
end


"Plot measured IV in 2D"
function plot_iv2d(gl; xyscale=1e3, outlines=true, fsize=(1000,700), logiv=false, kwargs...)
    iv = copy(gl.iv.vals.v)
    iv[.!gl.ivmask] = NaN
    iv[iv.<0] = NaN
    if logiv
        iv = log10(iv)
    end

    h1 = contour(gl.iv.vals.x/xyscale, gl.iv.vals.y/xyscale, iv', title="", fill=true,
                 aspect_ratio=:equal, colorbar_title="v (m/a)", kwargs...)

    return h1
end


"""
     plot2d_h(gl, hs2d_ivs2d::Matrix; clims=[:model,:radar][1], xyscale=10^3,
                  outlines=false, fsize=(1000,700), plotiv=false, plottraj=true,
                  plotline=!plotiv, flip=false,
                  seriescolor=[:inferno, :pu_or][1],
                  colorbar=true,
                  title = "",
                  colorbar_title = plotiv ? "iv (m/a)" : "h (m)",
                  kwargs...)

Plots modeled thickness and thickness along radar lines.
"""
function plot2d_h(gl, hs2d_ivs2d::Matrix; clims=[:model,:radar][1], xyscale=10^3,
                  outlines=false, fsize=(1000,700), plotiv=false, plottraj=true,
                  plotline=!plotiv, flip=false,
                  seriescolor=[:inferno, :pu_or][1],
                  colorbar=true,
                  title = "",
                  colorbar_title = plotiv ? "iv (m/a)" : "h (m)",
                  kwargs...)

    units = getunits(xyscale)
    if plotiv
        tr = gl.iv.vals
    else
        tr = gl.h.vals
    end
    xs,ys = flip ? (gl.dem.y/xyscale, gl.dem.x/xyscale) : (gl.dem.x/xyscale,gl.dem.y/xyscale)
    hs = flip ? copy(hs2d_ivs2d) : copy(hs2d_ivs2d)'
    xl,yl = flip ? ("y $units","x $units") : ("x $units", "y $units")

    if clims==:radar
        max_h = maximum(tr.v)
        clims = (0,max_h)
    elseif clims==:model
        clims = VAWTools.extremanan(hs2d_ivs2d)
    else
        clims = clims
    end
    clims::Tuple{Number,Number}

    hs[hs.==0] = NaN

    h1 = heatmap(xs, ys, hs; title="", #fill=true,
                 aspect_ratio=:equal,
                 clims=clims,
                 xlabel=xl,
                 ylabel=yl,
                 xflip=flip,
                 seriescolor=seriescolor,
                 colorbar=colorbar,
                 title=title,
                 colorbar_title=colorbar_title,
                 size=fsize,
                 kwargs...)

    if isa(tr, Traj{F}) && length(tr)!=0
        # don't plot more than 10000 points:
        st = max(length(tr.x)÷10000,1)
        xs,ys = flip ? (tr.y/xyscale, tr.x/xyscale) : (tr.x/xyscale, tr.y/xyscale)
        for (i, inds_) in enumerate(tr.splits) # go trough each segement
            inds = inds_[1]:st:inds_[end]
            if plottraj
                # add white border:
                scatter!(xs[inds], ys[inds], zcolor=tr.v[inds], markerstrokewidth=2, label="", colorbar=colorbar,
                         clims=clims, markerstrokecolor=:white)
                scatter!(xs[inds], ys[inds], zcolor=tr.v[inds], markerstrokewidth=0, label="", colorbar=colorbar,
                         clims=clims)
            end
            # plot a white line on top
            if false # plotiv || !(gl.gid==ITMIXGlacier(:SouthGlacier,1) || gl.gid==ITMIXGlacier(:NorthGlacier,1))
                # there are issues with the data there...
                if plotline==true
                    plot!(xs[inds], ys[inds], color=:white, label="", lw=1)
                elseif plotline isa AbstractVector
                    if i in plotline
                        # only plot if line in in plotline
                        plot!(xs[inds], ys[inds], color=:white, label="", lw=1.0)
                    end
                elseif plotline==false
                else
                    error()
                end
            end

            # runs into https://github.com/tbreloff/Plots.jl/issues/412
        end
    end

    outlines && plot_outlines!(gl, xyscale=xyscale, flip=flip)

    return h1
end
