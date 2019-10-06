 # ESA PVIR specific stuff

include("itmix-setup.jl")

sd = "esa-14nov"
data, data_summarized = FileIO.load("data/ITMIX2/results_Werder/$sd/data.jld", "data", "data_summarized");

if !isdefined(:gls)
    gls = open("data/ITMIX2/results_Werder/$sd/gls.jls", "r") do io
        deserialize(io)
    end
end
runit = false
include("itmix2-evaluate.jl")

using Mustache
"""
Turn a DataFrame into a latex table.  Copied from Mustache.jl README.
"""
function df_to_table(df; label="label", caption="caption", fmt=repeat("c", length(df)),
                     digits=4, header=replace(join(names(df), " & ")*"\\\\ \\hline",
                                              "_", "\\_")
                     )
    # row = join(["{{$x}}" for x in map(string, names(df))], " & ")

    # tpl="""
    # \\begin{table}
    #   \\centering
    #   \\begin{tabular}{$fmt}
    #   $header
    # {{#:DF}}    $row\\\\
    # {{/:DF}}  \\end{tabular}
    #   \\caption{$caption}
    #   \\label{tab:$label}
    # \\end{table}
    # """
    # render(tpl, DF=df)

    df_to_table_header(df, fmt=fmt, header=header) *
    df_to_table_body(df, digits) *
    df_to_table_footer(df, label=label, caption=caption)
end

function df_to_table_footer(df; label="label", caption="caption")
    tpl="""
      \\end{tabular}
      \\caption{$caption}
      \\label{tab:$label}
    \\end{table}
    """
    render(tpl, DF=df)
end

function df_to_table_header(df; fmt=repeat("c", length(df)), header="")
    tpl="""
    \\begin{table}
      \\centering
      \\begin{tabular}{$fmt}
        $header
    """
end

function df_to_table_body(df, digits)
    df = signif_df(df, digits)

    row = join(["{{$x}}" for x in map(string, names(df))], " & ")

    tpl="""
    {{#:DF}}    $row\\\\
    {{/:DF}}"""
    render(tpl, DF=df)
end

function signif_df(df, digits=3)
    out = DataFrame()
    colwise(x -> x isa AbstractFloat ? signif(x, digits) : x, df)

    for (n,c) in zip(names(df), DataFrames.columns(df))
        if eltype(c) <: AbstractFloat
            out[n] = signif.(c, (digits,))
        else
            out[n] = c
        end
    end
    out
end


# Compare tables, for exp0, 4, 10, 15 of U, S1, A
# exp4: longitudinal
# exp 10 50% retained
# exp 15 20% retained
exps = [0,4,10,15]
gls_PVIR = [:Unteraar, :Synthetic1, :Austfonna, :Tasman]

"""
Make a table like
```
Glacier & RMSE  & no IV & max err. & no IV  & improv. & ΔRMSE & eff_sample_size \\
```
Note that the RMSE is from the control radar lines only.
"""
function make_table0(data, sig=4)
    error("Remove Rhat")
    first = true
    last = false
    tab = 1
    for e in exps
        !first && println("    Exp $e\\\\")
        for gl in gls_PVIR
            tab = by(data[(data[:expnr].==e) .& (data[:gl].==gl),:], :gl,
                     d -> DataFrame(rmse = signif.(d[d[:rm_iv],:rmse], sig),
                                    rms_rmiv = signif.(d[.!d[:rm_iv],:rmse], sig),
                                    max = signif.(d[d[:rm_iv],:max], sig),
                                    max_rmiv = signif.(d[.!d[:rm_iv],:max], sig),
                                    iv_helps = d[d[:rm_iv],:rmse] .< d[.!d[:rm_iv],:rmse],
                                    iv_rmse_change = signif.(d[d[:rm_iv],:rmse] .- d[.!d[:rm_iv],:rmse], sig),
                                    Rhat = signif.(median(d[:Rhat]), 2)
                                    ))
            if first
                print(df_to_table_header(tab,
                                         fmt="l"*repeat("r", length(tab)-1),
                                         header="Glacier & RMSE & w/o IV & max err. & w/o IV & improv. & \$\\Delta\$ RMSE & \$\\hat{R}\$\\\\ \\hline"
                                         ))
                println("    Exp $e\\\\")
                first = false
            end
            print(df_to_table_body(tab))
        end
    end
    print(df_to_table_footer(tab,
                             label="tab:sum-$e",
                             caption="Summary of results for ITMIX2 experiments $exps for selected glaciers with IV-data.",
                             ))
end

function make_tables1(data, sig=4)
    error("remove Rhat")
    for e in exps
        tab = by(data[(data[:expnr].==e),:], :gl,
                 d -> DataFrame(rmse = signif.(d[d[:rm_iv],:rmse], sig),
                                rms_rmiv = signif.(d[.!d[:rm_iv],:rmse], sig),
                                max = signif.(d[d[:rm_iv],:max], sig),
                                max_rmiv = signif.(d[.!d[:rm_iv],:max], sig),
                                iv_helps = d[d[:rm_iv],:rmse] .< d[.!d[:rm_iv],:rmse],
                                iv_rmse_change = signif.(d[d[:rm_iv],:rmse] .- d[.!d[:rm_iv],:rmse], sig),
                                Rhat = signif.(median(d[:Rhat]), 2)
                                ))
        print(df_to_table(tab,
                          label="tab:sum-$e",
                          fmt="l"*repeat("r", length(tab)-1),
                          caption="Summary of results for ITMIX2 experiment $e for glaciers with IV-data.",
                          header="Glacier & RMSE & w/o IV & max err. & w/o IV & improv. & \$\\Delta\$ RMSE & \$\\hat{R}\$\\\\ \\hline")
              )
    end
end

# esa-14nov
#  expnr = 0
# 11×7 DataFrames.DataFrame
# │ Row │ gl             │ rmse  │ rms_rmiv │ max    │ max_rmiv │ iv_helps │ iv_rmse_change │
# ├─────┼────────────────┼───────┼──────────┼────────┼──────────┼──────────┼────────────────┤
# │ 1   │ Synthetic1     │ 39.1  │ 50.3     │ 145.0  │ 172.0    │ true     │ -11.2          │
# │ 2   │ Unteraar       │ 86.0  │ 145.0    │ 255.0  │ 393.0    │ true     │ -59.0          │
# │ 3   │ Austfonna      │ 226.0 │ 136.0    │ 1050.0 │ 686.0    │ false    │ 90.0           │
# │ 4   │ SouthGlacier   │ 61.5  │ 38.4     │ 164.0  │ 93.7     │ false    │ 23.1           │
# │ 5   │ Synthetic2     │ 131.0 │ 135.0    │ 308.0  │ 318.0    │ true     │ -4.0           │
# │ 6   │ Brewster       │ 56.3  │ 53.6     │ 99.0   │ 92.5     │ false    │ 2.7            │
# │ 7   │ Devon          │ 337.0 │ 523.0    │ 1190.0 │ 1930.0   │ true     │ -186.0         │
# │ 8   │ NorthGlacier   │ 148.0 │ 166.0    │ 375.0  │ 400.0    │ true     │ -18.0          │
# │ 9   │ Synthetic3     │ 64.5  │ 47.1     │ 231.0  │ 160.0    │ false    │ 17.4           │
# │ 10  │ Hellstugubreen │ 78.1  │ 79.1     │ 146.0  │ 146.0    │ true     │ -1.0           │
# │ 11  │ Tasman         │ 295.0 │ 223.0    │ 450.0  │ 411.0    │ false    │ 72.0           │

#  expnr = 4
# 11×7 DataFrames.DataFrame
# │ Row │ gl             │ rmse  │ rms_rmiv │ max    │ max_rmiv │ iv_helps │ iv_rmse_change │
# ├─────┼────────────────┼───────┼──────────┼────────┼──────────┼──────────┼────────────────┤
# │ 1   │ Synthetic1     │ 41.0  │ 40.2     │ 126.0  │ 122.0    │ false    │ 0.8            │
# │ 2   │ Unteraar       │ 82.0  │ 83.3     │ 246.0  │ 248.0    │ true     │ -1.3           │
# │ 3   │ Austfonna      │ 151.0 │ 131.0    │ 804.0  │ 611.0    │ false    │ 20.0           │
# │ 4   │ SouthGlacier   │ 21.4  │ 23.3     │ 66.7   │ 75.1     │ true     │ -1.9           │
# │ 5   │ Synthetic2     │ 25.3  │ 23.8     │ 89.3   │ 80.9     │ false    │ 1.5            │
# │ 6   │ Brewster       │ 26.0  │ 26.3     │ 65.7   │ 66.0     │ true     │ -0.3           │
# │ 7   │ Devon          │ 200.0 │ 259.0    │ 1000.0 │ 1670.0   │ true     │ -59.0          │
# │ 8   │ NorthGlacier   │ 35.5  │ 38.5     │ 98.2   │ 124.0    │ true     │ -3.0           │
# │ 9   │ Synthetic3     │ 21.6  │ 21.6     │ 63.0   │ 66.7     │ false    │ 0.0            │
# │ 10  │ Hellstugubreen │ 28.8  │ 29.1     │ 71.8   │ 72.0     │ true     │ -0.3           │
# │ 11  │ Tasman         │ 180.0 │ 251.0    │ 201.0  │ 398.0    │ true     │ -71.0          │

#  expnr = 10
# 11×7 DataFrames.DataFrame
# │ Row │ gl             │ rmse  │ rms_rmiv │ max   │ max_rmiv │ iv_helps │ iv_rmse_change │
# ├─────┼────────────────┼───────┼──────────┼───────┼──────────┼──────────┼────────────────┤
# │ 1   │ Synthetic1     │ 28.2  │ 27.8     │ 54.7  │ 53.8     │ false    │ 0.4            │
# │ 2   │ Unteraar       │ 62.1  │ 62.1     │ 164.0 │ 165.0    │ false    │ 0.0            │
# │ 3   │ Austfonna      │ 72.8  │ 73.2     │ 288.0 │ 291.0    │ true     │ -0.4           │
# │ 4   │ SouthGlacier   │ 20.8  │ 21.2     │ 65.4  │ 66.8     │ true     │ -0.4           │
# │ 5   │ Synthetic2     │ 25.9  │ 26.0     │ 61.2  │ 61.5     │ true     │ -0.1           │
# │ 6   │ Brewster       │ 27.5  │ 27.3     │ 50.8  │ 49.9     │ false    │ 0.2            │
# │ 7   │ Devon          │ 137.0 │ 138.0    │ 499.0 │ 487.0    │ true     │ -1.0           │
# │ 8   │ NorthGlacier   │ 35.1  │ 34.9     │ 88.1  │ 88.7     │ false    │ 0.2            │
# │ 9   │ Synthetic3     │ 22.6  │ 22.7     │ 70.3  │ 70.8     │ true     │ -0.1           │
# │ 10  │ Hellstugubreen │ 28.6  │ 28.7     │ 96.8  │ 96.7     │ true     │ -0.1           │
# │ 11  │ Tasman         │ 206.0 │ 230.0    │ 280.0 │ 255.0    │ true     │ -24.0          │

#  expnr = 15
# 11×7 DataFrames.DataFrame
# │ Row │ gl             │ rmse  │ rms_rmiv │ max   │ max_rmiv │ iv_helps │ iv_rmse_change │
# ├─────┼────────────────┼───────┼──────────┼───────┼──────────┼──────────┼────────────────┤
# │ 1   │ Synthetic1     │ 31.7  │ 31.7     │ 120.0 │ 119.0    │ false    │ 0.0            │
# │ 2   │ Unteraar       │ 52.7  │ 52.5     │ 184.0 │ 184.0    │ false    │ 0.2            │
# │ 3   │ Austfonna      │ 66.6  │ 66.8     │ 288.0 │ 291.0    │ true     │ -0.2           │
# │ 4   │ SouthGlacier   │ 19.9  │ 20.4     │ 67.6  │ 66.7     │ true     │ -0.5           │
# │ 5   │ Synthetic2     │ 23.5  │ 23.3     │ 52.7  │ 52.5     │ false    │ 0.2            │
# │ 6   │ Brewster       │ 26.0  │ 23.7     │ 54.7  │ 50.0     │ false    │ 2.3            │
# │ 7   │ Devon          │ 152.0 │ 160.0    │ 479.0 │ 451.0    │ true     │ -8.0           │
# │ 8   │ NorthGlacier   │ 37.3  │ 37.9     │ 99.7  │ 103.0    │ true     │ -0.6           │
# │ 9   │ Synthetic3     │ 22.7  │ 22.9     │ 79.2  │ 78.1     │ true     │ -0.2           │
# │ 10  │ Hellstugubreen │ 40.2  │ 41.0     │ 73.8  │ 74.6     │ true     │ -0.8           │
# │ 11  │ Tasman         │ 224.0 │ 220.0    │ 350.0 │ 313.0    │ false    │ 4.0            │


using Plots

function plot_PVIR_2x8(name, expnr, subdir=sd;
                       reuse=true,
                       xyscale=1e3,
                       kwargs...)
    gid = BM.ITMIXGlacier(name,2)
    gl, gb, pp, pm, pn, pl = gls[gid];
    ps = Array{Any}(2,4)

    # with IV
    for i=1:2
        xf = i==2 ? x->round(Int, x) : x->""
        xl = i==2 ? "χ (km)" : ""

        # https://discourse.julialang.org/t/avoid-repeating-tick-labels-with-subplots-in-plots-jl/351/4
        data_out, hs2d, hs2ds, iv2d, iv2ds, r_lines = load_one(gid, expnr, [true,false][i]);

        ps[i,1] = BM.plot2d_h(gl, hs2d.v, xyscale=xyscale, colorbar=false,
                              xformatter=xf, xlabel=xl)
        ps[i,2] = BM.plot2d_h(gl, hs2ds.v, xyscale=xyscale, colorbar=false,
                              xformatter=xf, yformatter=x->"", xlabel=xl, ylabel="")
        ps[i,3] = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, iv2d.v', size=(600,600), #, zlims=(0,50))
                          aspect_ratio=:equal, colorbar=false,
                          title= i==1 ? "iv" :"",
                          xformatter=xf, xlabel=xl,
                          kwargs...);
        ps[i,4] = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, iv2ds.v', size=(600,600), #, zlims=(0,50))
                          aspect_ratio=:equal, colorbar=false,
                          title= i==1 ? "std(iv)" :"",
                          xformatter=xf, yformatter=x->"", xlabel=xl, ylabel="",
                          kwargs...);
    end
    c1 = scatter([0,0], [0,1], zcolor=[clims...], clims=clims,# aspect_ratio=:equal,
                 xlims=(1,1.1), axis=false, label="", colorbar_title=tit,
                 seriescolor=seriescolor
                 )#, markeralpha=0, markerstrokewidth=0)

    l = @layout [grid(2,2) a{0.01w} grid(2,2) a{0.01w}]
    plot(permutedims(ps, [2,1])[1:4]..., c1, permutedims(ps, [2,1])[5:8]..., c2, layout=l, reuse=reuse; size=(900,600), kwargs...)
end

# ## hack to place a colorbar willy-nilly
# ## https://github.com/tbreloff/Plots.jl/issues/412
# h1 = contour(rand(3,3), cbar=false, clims=(0,3), fill=true)
# scatter!(1:3,1:3, zcolor=1:3, cbar=false, clims=(0,3))
# h3 = contour(rand(3,3), cbar=false, clims=(0,3), fill=true)
# h2 = scatter([0,0], [0,1], zcolor=[0,3], clims=(0,3),
#              xlims=(1,1.1), axis=false, grid=false, label="", colorbar_title="cbar")
# l = @layout [grid(2,1) a{0.01w}]
# plot(h1,h3,h2,layout=l)



#     h2 = scatter([0,0], [0,1], zcolor=[clims...], clims=clims,# aspect_ratio=:equal,
#                  xlims=(1,1.1), axis=false, label="", colorbar_title=tit,
#                  seriescolor=seriescolor
#                  )#, markeralpha=0, markerstrokewidth=0)
#     l = @layout [grid(1,1) a{0.01w}]
#     plot(h1,h2,layout=l, size=fsize, palette=:viridis; kwargs...)

function plot_PVIR_2x2(name, expnr, subdir=sd;
                       plotiv=false,
                       reuse=true,
                       xyscale=1e3,
                       seriescolor=[:inferno, :pu_or][1],
                       fsize=(1800,900),
                       kwargs...)
    gid = BM.ITMIXGlacier(name,2)
    gl, gb, pp, pm, pn, pl = gls[gid];
    ps = Array{Any}(2,3)

    clims_h = (1, 1.0)
    clims_sh = (0, 1e-5)
    clims_err = (-1e-5, 1e-5)
    # just to the plot limits:
    for rm_iv = [true,false]
        data_out, hs2d, hs2ds, iv2d, iv2ds, r_lines_fitting = load_one(gid, expnr, rm_iv);
        d, ds = plotiv ? (iv2d, iv2ds) : (hs2d, hs2ds)
        rd = plotiv ? gl.iv.vals : gl.h.vals
        err, err_ = if plotiv
            err = BM.error2d(d.v, gl, :iv, gl.ivmask)
            (err, VAWTools.extremanan(err))
        else
            err = BM.error_h(d.v, gl)
            (err, VAWTools.extremanan(err.v))
        end
        h = VAWTools.extremanan(d.v)
        sh = VAWTools.extremanan(ds.v)
        clims_h = (min(clims_h[1], h[1]), max(clims_h[2], h[2]))
        clims_sh = (min(clims_sh[1], sh[1]), max(clims_sh[2], sh[2]))
        clims_err = (min(clims_err[1], err_[1]), max(clims_err[2], err_[2]))
    end

    labels = ["Exp. $expnr with IV", "Exp. $expnr w/o IV"]
    letters = [["A", "B", "C"], ["D", "E", "F"]]
    # plot with IV and without:
    for i=1:2 # rows, corresponding to rm_iv false/true
        xf = i==2 ? x->round(Int, x) : x->""
        xl = i==2 ? "χ (km)" : ""
        # https://discourse.julialang.org/t/avoid-repeating-tick-labels-with-subplots-in-plots-jl/351/4
        data_out, hs2d, hs2ds, iv2d, iv2ds, r_lines_fitting = load_one(gid, expnr, [true,false][i]);
        d, ds = plotiv ? (iv2d, iv2ds) : (hs2d, hs2ds)
        rd = plotiv ? gl.iv.vals : gl.h.vals
        err, rmse, linf = if plotiv
            err = BM.error2d(d.v, gl, :iv, gl.ivmask)
            vv = err[.!isnan.(err)]
            rmse = round(Int, sqrt(sum(vv.^2)/length(vv)) )
            (err, rmse, round(Int, maximum(abs.(vv))) )
        else
            err = BM.error_h(d.v, gl)
            vv = err.v[.!isnan.(err.v)]
            # RMSE, max-error on control lines
            rmse, maxerr = BM.compare_to_radar(gl, d, ds, r_lines_fitting, thin= (name==:Tasman ? 1 : 10)) # TODO hack...
            # maxerr = (maximum(abs.(vv)),)

            (err, round(Int, rmse[1]), round(Int, maxerr[1]))
        end

        ## Plot value
        ps[i,1] = BM.plot2d_h(gl, d.v, xyscale=xyscale, colorbar=false,
                              plotiv=plotiv,
                              xformatter=xf, xlabel=xl, clims=clims_h,
                              seriescolor=seriescolor,
                              plotline=plotiv ? false : I2.get_rlines(gid, expnr),
                              label=labels[i],
                              outlines=true
                              )
        # add a legend
        scatter!([d.x[1]/xyscale], [d.y[1]/xyscale], label="_", legendtitle=labels[i], markeralpha=0, markerstrokewidth=0,
                 xformatter=xf, xlabel=xl)

        ## Plot std
        ps[i,2] = BM.plot2d_h(gl, ds.v, xyscale=xyscale, colorbar=false,
                              plotiv=plotiv,
                              xformatter=xf, yformatter=x->"", xlabel=xl, ylabel="",
                              clims=clims_sh,
                              plotline=plotiv ? I2.get_rlines(gid, expnr) : false,
                              seriescolor=seriescolor,
                              #plotline=true,
                              plottraj=false,
                              outlines=true)

        ## Plot error
        lt = if rd isa Traj
            # plot error on the radar tracks
            ps[i,3] = Plots.plot();
            BM.plot_traj!(err; xyscale=xyscale, plotline=I2.get_rlines(gid, expnr),
                          colorbar=false, aspect_ratio=:equal,
                          yformatter=x->"", ylabel="",
                          xformatter=xf, xlabel=xl,
                          clims=clims_err)
           "RMSE = $(rmse)m\nMax err = $(linf)m"
        else
            # this does not work for point measurements, such as SouthGlacier
            err[.!gl.ivmask] = NaN
            ps[i,3] = heatmap(gl.dem.x/xyscale, gl.dem.y/xyscale, err',
                              colorbar=false, aspect_ratio=:equal,
                              clims=clims_err,
                              kwargs...);
            "RMSE = $(rmse)m/a\nMax err = $(linf)m/a"
        end
        # make a label
        scatter!([d.x[1]/xyscale], [d.y[1]/xyscale], label="_", legendtitle=lt,
                 markeralpha=0, markerstrokewidth=0,
                 xformatter=xf, xlabel=xl, yformatter=x->"")
        BM.plot_outlines!(gl, xyscale=xyscale, ylabel="", xlabel=xl)
    end
    xls,yls = xlims(), ylims()
    for p in ps
        xlims!(p, xls)
        ylims!(p, yls)
    end
    # ## hack to place a colorbar willy-nilly
    # ## https://github.com/tbreloff/Plots.jl/issues/412

    c1 = scatter([0,0], [0,1], zcolor=[clims_h...], clims=clims_h,# aspect_ratio=:equal,
                 xlims=(1,1.1), axis=false, label="", colorbar_title=plotiv ? "iv (m/a)" : "h (m)",
                 seriescolor=seriescolor, grid=false
                 )#, markeralpha=0, markerstrokewidth=0)
    c2 = scatter([0,0], [0,1], zcolor=[clims_sh...], clims=clims_sh,# aspect_ratio=:equal,
                 xlims=(1,1.1), axis=false, label="", colorbar_title=plotiv ? "std(iv) (m/a)" : "std(h) (m)",
                 seriescolor=seriescolor, grid=false
                 )#, markeralpha=0, markerstrokewidth=0)
    c3 = scatter([0,0], [0,1], zcolor=[clims_err...], clims=clims_err,# aspect_ratio=:equal,
                 xlims=(1,1.1), axis=false, label="", colorbar_title=plotiv ? "error(iv) (m/a)" : "error(h) (m)",
                 seriescolor=seriescolor, grid=false
                 )#, markeralpha=0, markerstrokewidth=0)

    l = @layout [grid(2,1) a{0.01w} grid(2,1) b{0.01w} grid(2,1) b{0.01w}]
    # l = @layout [grid(2,1) [a{0.01w}
    #                         b{0.01w}
    #                         c{0.01w}] grid(2,1) d{0.01w}]
    #l = @layout [grid(2,1) grid(3,1,widths=[0.01,0.01,0.01]) grid(2,1) d{0.01w}]
    plot(ps[1:2]..., c1, ps[3:4]..., c2, ps[5:6]..., c3, layout=l, reuse=reuse; size=fsize, kwargs...)
end

# dir_ = "/home/mauro/projects/ESA-ice-volume/BITEModel/scripts/data/ITMIX2/results_Werder/esa-14nov/plots/esa/"
# error()
# plotiv = false
# for name in gls_PVIR
#     @show name
#     for e in exps
#         ee = e<10 ? "0$e" : "$e"
#         @show ee
#         fln = plotiv ? "iv" : "h"
#         plot_PVIR_2x2(name, e, plotiv=plotiv)
#         savefig(dir_*"$(fln)_$name-$ee.png")
# #        error()
#     end
# end


############
# get Matthias' results
function get_matthias()
    data_huss = make_data()

    # ITMIX2
    dir = "/home/mauro/projects/ESA-ice-volume/reporting/n6_Product_Validation_Intercomparison/itmix2-results-matthias"
    for name in map(BM.getid, keys(gls))
        gid = BM.ITMIXGlacier(name,2)
        for e in 1:16
            push!(data_huss, read_others_results(gid, e, dir, "Huss"))
        end
    end

    # ITMIX1
    dir = "/home/mauro/projects/ESA-ice-volume/reporting/n6_Product_Validation_Intercomparison/itmix1-results-matthias"
    for name in map(BM.getid, keys(gls))
        name in [:Austfonna,:Devon] && continue
        gid = BM.ITMIXGlacier(name,2)
        for e in 0
            push!(data_huss, read_others_results(gid, e, dir, "Huss"))
            # if name==:Austfonna
            #     push!(data_huss, read_others_results(gid, e, dir, "Huss", thin=6, start=(190, 523), stop=(2927, ))
            # else
        #     push!(data_huss, read_others_results(gid, e, dir, "Huss"))
            # end
        end
    end
    sort!(data_huss)
end

data_huss = get_matthias()
## Compare
glids0 = [BM.getid(k) for k in keys(gls) if !(BM.getid(k) in [:Austfonna, :Devon])]
glids = [BM.getid(k) for k in keys(gls)]

data_huss0 = data_huss[in.(data_huss[:gl], (glids0,)), :];
huss = summarize_data(data_huss0, expnrs=0:16, person=true);
huss0 = summarize_data(data_huss0, expnrs=0, person=true);
huss1_16 = summarize_data(data_huss, expnrs=1:16, person=true);
huss4 = summarize_data(data_huss, expnrs=4, person=true);

data0 = data[in.(data[:gl], (glids0,)), :];
me = summarize_data(data0[data0[:rm_iv].==false,:], expnrs=0:16, person=true);
me0 = summarize_data(data0[data0[:rm_iv].==false,:], expnrs=0, person=true);
me1_16 = summarize_data(data[data[:rm_iv].==false,:], expnrs=1:16, person=true);
me4 = summarize_data(data[data[:rm_iv].==false,:], expnrs=4, person=true);

me_rmiv = summarize_data(data0[data0[:rm_iv].==true,:], expnrs=0:16, person=true);
me0_rmiv = summarize_data(data0[data0[:rm_iv].==true,:], expnrs=0, person=true);
me1_16_rmiv = summarize_data(data[data[:rm_iv].==true,:], expnrs=1:16, person=true);

function compare_two(summarized_huss, summarized_me,
                     median_mean_h = (sort(summarized_huss)[:median_mean_h] .+ sort(summarized_me)[:median_mean_h])/2)
    huss, me = sort(summarized_huss), sort(summarized_me)
    fls = [(:rel_mean_rmse, :mean_rmse)]
    out = DataFrame()
    out[:gl] = huss[:gl]
    out[:median_mean_h_huss] = huss[:median_mean_h]
    out[:median_mean_h_me] = me[:median_mean_h]
    out[:median_mean_h] = median_mean_h

    out[:rel_mean_rmse_huss] = huss[:mean_rmse]./median_mean_h
    out[:rel_mean_rmse_me] = me[:mean_rmse]./median_mean_h

    for (flo,fl) in fls
        rmse_diff = (huss[fl] .- me[fl])./median_mean_h
        out[flo] = rmse_diff
        out[:better] = significatly_better.(rmse_diff, 0.05)
    end
    out
end

function significatly_better(diff, thresh)
    if diff>thresh
        return "true"
    elseif diff<-thresh
        return "false"
    else
        return "-"
    end
end

# tables of Intercomparison
df_to_table(compare_two(huss1_16, me1_16)[[:gl, :median_mean_h_huss, :median_mean_h_me, :rel_mean_rmse_huss, :rel_mean_rmse_me, :better]]
                                           , label="comp1-16", caption="Comparing ITMIX2 results between Huss and Werder",
            fmt="lrrrrr", digits=3,
            header="Glacier & \$\\bar{h}_{Huss}\$ & \$\\bar{h}_{BITE}\$ & \$\\overline{\\rm RMSE}_{H\\&F}\$ & \$\\overline{\\rm RMSE}_{BITE}\$ & better\\\\ \\hline") |> print

df_to_table(compare_two(huss4, me4)[[:gl, :median_mean_h_huss, :median_mean_h_me, :rel_mean_rmse_huss, :rel_mean_rmse_me, :better]]
                                           , label="comp4", caption="Comparing ITMIX2 results between Huss and Werder for Exp04",
            fmt="lrrrrr", digits=3,
            header="Glacier & \$\\bar{h}_{Huss}\$ & \$\\bar{h}_{BITE}\$ & \$\\overline{\\rm RMSE}_{H\\&F}\$ & \$\\overline{\\rm RMSE}_{BITE}\$ & better\\\\ \\hline") |> print


df_to_table(compare_two(huss0, me0,
                        (sort(huss0)[:median_mean_h] .+ sort(me0)[:median_mean_h])/2
                        )[[:gl, :median_mean_h_huss, :median_mean_h_me, :rel_mean_rmse_huss, :rel_mean_rmse_me, :better]]
                                           , label="comp4", caption="Comparing ITMIX results between Huss and Werder (i.e. Exp00)",
            fmt="lrrrrr", digits=3,
            header="Glacier & \$\\bar{h}_{Huss}\$ & \$\\bar{h}_{BITE}\$ & \$\\overline{\\rm RMSE}_{H\\&F}\$ & \$\\overline{\\rm RMSE}_{BITE}\$ & better\\\\ \\hline") |> print

df_to_table(compare_two(huss0, me0_rmiv,
                        (sort(huss0)[:median_mean_h] .+ sort(me0)[:median_mean_h])/2
                        )[[:gl, :median_mean_h_huss, :median_mean_h_me, :rel_mean_rmse_huss, :rel_mean_rmse_me, :better]]
                                           , label="comp4", caption="Comparing ITMIX results between Huss and Werder without using IV (i.e. Exp00)",
            fmt="lrrrrr", digits=3,
            header="Glacier & \$\\bar{h}_{Huss}\$ & \$\\bar{h}_{BITE}\$ & \$\\overline{\\rm RMSE}_{H\\&F}\$ & \$\\overline{\\rm RMSE}_{BITE}\$ & better\\\\ \\hline") |> print


function get_mean_thicknesses(data)
    d = summarize_data(data)[[:gl, :median_mean_h]]
    Dict(k=>v for (k,v) in zip(d[:gl], d[:median_mean_h]))
end


function plot_whiskers(data_huss, data_me, data_mean_h=vcat(data_huss, data_me); kw...)
    # join
    all = vcat(data_huss, data_me);
    all[:gl] = CategoricalArray(string.(all[:gl]));

    # get median mean thicknesses
    mts = get_mean_thicknesses(data_mean_h)

    # explode (manual stack)
    # need [person, gl, expnr, err_h]
    #    T = CategoricalArrays.CategoricalArray{String,1,UInt32,String,CategoricalArrays.CategoricalString{UInt32},Union{}}
    T = CategoricalArrays.CategoricalString{UInt32}
    newer = DataFrame([T, T, Int64, Float64], [:person, :gl, :rm_iv, :err_h], 0)

    for g in groupby(all, [:person, :gl, :rm_iv])
        person, gl, rm_iv = g[:person][1], g[:gl][1], g[:rm_iv][1]
        person = person==:Huss ? "H&F" : "BITE"
        gl = gl=="Hellstugubreen" ? "Hellstugub." : gl
        mt = mts[Symbol(g[:gl][1])]
        err_h = filter(x->!isnan(x), vcat([gg.v for gg in g[:err_h]]...))
        newer = vcat(newer, DataFrame(person=person, gl=gl, rm_iv=rm_iv, err_h=err_h))
    end

    StatPlots.@df newer boxplot(:gl, :err_h, group=:person; kw...)
end

dir_ = "/home/mauro/projects/ESA-ice-volume/BITEModel/scripts/data/ITMIX2/results_Werder/esa-14nov/plots/esa/"

plot_whiskers(data_huss0[data_huss0[:expnr].==0,:],
              data0[(data0[:expnr].==0) .&  (data0[:rm_iv].==false),:], vcat(data_huss, data);
              ylims=(-250,250), ylabel="rel. error h (m)",
              size=(900,600), yticks=-500:50:500);
savefig(dir_*"intercomp-huss-exp00.png")

plot_whiskers(data_huss0[data_huss0[:expnr].==0,:],
              data0[(data0[:expnr].==0) .&  (data0[:rm_iv].==true),:], vcat(data_huss, data);
              ylims=(-250,250), ylabel="rel. error h (m)",
              size=(900,600), yticks=-500:50:500);
savefig(dir_*"intercomp-huss-exp00-rmiv.png")

plot_whiskers(data_huss[data_huss[:expnr].==4,:],
              data[(data[:expnr].==4) .&  (data[:rm_iv].==false),:], vcat(data_huss, data);
              ylims=(-250,250), ylabel="rel. error h (m)",
              size=(900,600), yticks=-500:50:500);
savefig(dir_*"intercomp-huss-exp04.png")


dh = data_huss[ in.(data_huss[:expnr], (1:16,)), :];
dm = data[ in.(data[:expnr], (1:16,))  .&  (data[:rm_iv].==false), :];
plot_whiskers(dh, dm,
              vcat(data_huss, data);
              ylims=(-250,250), ylabel="rel. error h (m)",
              size=(900,600), yticks=-500:50:500);
savefig(dir_*"intercomp-huss-exp01-16.png")
