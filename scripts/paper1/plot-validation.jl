# Plots and tables of the validation runs for Svalbard
using Missings, DataFrames, Measures
using Plots, StatPlots

saveyes = true
scaleyes = true

include("load-validation.jl")
#######
# Table helpers
include("../mustache-tables.jl")
#######
#
# use median absolute deviation MAD https://en.wikipedia.org/wiki/Median_absolute_deviation

function proc_validation(region,
                         out_fit, out_test, out_test_fitgls)
    # mean ice thicknesses
    scales = Dict(3=>235,
                  4=>172,
                  7=>210,
                  11=>68)
    scale = scaleyes ? scales[region] : 1


    fit_MAD = OrderedDict()
    fit_meanerr = OrderedDict()
    for ft in keys(out_fit)
        fit_MAD[ft] = [o[:abserr_h][:median] for o in out_fit[ft] if o isa Dict]./scale
        fit_meanerr[ft] = [o[:err_h][:mean] for o in out_fit[ft] if o isa Dict]./scale
    end

    test_MAD = OrderedDict()
    test_meanerr = OrderedDict()
    for k in keys(out_test)
        test_MAD[k] = [o[:abserr_h][:median] for o in out_test[k] if o isa Dict]./scale
        test_meanerr[k] = [o[:err_h][:mean] for o in out_test[k] if o isa Dict]./scale
    end

    test_fitgls_MAD = OrderedDict()
    test_fitgls_meanerr = OrderedDict()
    for k in keys(out_test_fitgls)
        test_fitgls_MAD[k] = [o[:abserr_h][:median] for o in out_test_fitgls[k] if o isa Dict]./scale
        test_fitgls_meanerr[k] = [o[:err_h][:mean] for o in out_test_fitgls[k] if o isa Dict]./scale
    end

    return fit_MAD, fit_meanerr, test_MAD, test_meanerr, test_fitgls_MAD, test_fitgls_meanerr
end

function plot_whiskers(labels, data::Vector; kwargs...)
    @assert length(labels)==length(data)

    labels = vcat([fill(l, length(d)) for (l,d) in zip(labels,data)]...)
    dat = vcat(data...)
    boxplot(labels, dat; sorted=false, label="", kwargs...)
end
function plot_whiskers(data::Associative; kwargs...)
    labels = map(string, keys(data))
    plot_whiskers(labels, collect(values(data)))
#     labels = vcat([fill(string(k), length(v)) for (k,v) in data]...)
#     dat = vcat(collect(values(data))...)
# #    dat = hcat([v[1:72] for v in  values(data)]...)
#     boxplot(labels, dat; kwargs...)
end

function plot_supplement(region,
                         fit_MAD, fit_meanerr, test_MAD, test_meanerr, test_fitgls_MAD, test_fitgls_meanerr)

    # labels = ["F:lh-", "F:lhv", "F:l-v", "F:l--", "P:---" , "P:lh-", "P:l-v",
    #           "P:lhv", "PF:---", "PF:lhv"]
    labels = ["C:lh-", "C:lhv", "C:l-v", "C:l--", "V:l--;\n  ---" , "V:l--;\n  lh-", "V:l--;\n  l-v",
              "V:l--;\n  lhv", "V:l-v;\n  ---", "V:l-v;\n  lhv"]


    p1 = plot_whiskers(labels, collect(values(merge(fit_MAD, test_MAD))),
                       xtickfontfamily="monospace", ylabel="norm. MAD (m)");

    fs = 20
    if region==7
        plot!([4,4],[-0.1,2.2], lc=:black, lw=2, label="", ls=:dash)
        annotate!(0.25, 350, text("A", fs))
    elseif region==11
        plot!([4,4],[-0.1,2.4], lc=:black, lw=2, label="", ls=:dash, ylims=(-0.1,2.4))
        annotate!(0.25, 120, text("A", fs))
    elseif region==4
        plot!([4,4],[-0.1,2.75], lc=:black, lw=2, label="", ls=:dash) #, ylims=(-5,150))
        annotate!(0.25, 120, text("A", fs))
    elseif region==3
        plot!([4,4],[-0.1,2.1], lc=:black, lw=2, label="", ls=:dash) #, ylims=(-5,150))
        annotate!(0.25, 120, text("A", fs))
    end


    p2 = plot_whiskers(labels, collect(values(merge(fit_meanerr, test_meanerr))),
                       ylabel="norm. ME (m)", xtickfontfamily="monospace");

    if region==7
        plot!([4,4],[-1.1,2.2], lc=:black, lw=2, label="", ls=:dash)
        annotate!(0.25, 350, text("B", fs))
    elseif region==11
        plot!([4,4],[-2,2.1], lc=:black, lw=2, label="", ls=:dash, ylims=(-2,2.1))
        annotate!(0.25, 100, text("B", fs))
    elseif region==4
        plot!([4,4],[-1.1,2.4], lc=:black, lw=2, label="", ls=:dash) #, ylims=(-105,120))
        annotate!(0.25, 100, text("B", fs))
    elseif region==3
        plot!([4,4],[-2,1.9], lc=:black, lw=2, label="", ls=:dash) #, ylims=(-105,120))
        annotate!(0.25, 100, text("B", fs))
    end

    plot!([0.4,9.6],[0,0], lc=:black, lw=1, label="", ls=:dash, alpha=0.5)


    p = plot(p1,p2, layout=grid(2,1), size=(700,700))

    saveyes && savefig("figs/all_error-$(VAWTools.int2str2(region)).png")
    saveyes && savefig("figs/all_error-$(int2str2(region)).pdf")
    return p
end

############
# For paper make the plot easier by showing less
############

# Following https://github.com/JuliaPlots/StatsPlots.jl/pull/192#issuecomment-430239298 and @macroexpand
# gives
# StatPlots.add_label([], boxplot, xlabels, yvalues, group=grouplabels)

# labels = OrderedDict("F:lh-" =>1,
#                      "F:lhv" => 2,
#                      "F:l-v" => 3,
#                      "F:l--" => 4,
#                      "P:---"  => 5,
#                      "P:lh-" => 6,
#                      "P:l-v" => 7,
#                      "P:lhv" => 8,
#                      "PF:---" => 9,
#                      "PF:lhv" => 10]

function plot_validation(region,
                         fit_MAD, fit_meanerr, test_MAD, test_meanerr, test_fitgls_MAD, test_fitgls_meanerr;
                         saveyes=saveyes)
    # (index into data, xlabel (fit-type), group (prior-type)
    labels = [(2, "lhv", "calibration"),
              (4, "l--", "calibration"),
              (3, "l-v", "calibration"),
              (8, "l--", "validation"),
              (10, "l-v", "validation")]

    # make long form
    data_MAD = collect(values(merge(fit_MAD, test_MAD)))
    data_ME = collect(values(merge(fit_meanerr, test_meanerr)))
    datas_MAD = Float64[]
    datas_ME = Float64[]
    groups = String[]
    xlabs = String[]
    for (i, x, g)  in labels
        len = length(data_MAD[i])
        @assert len==length(data_ME[i])
        append!(datas_MAD, data_MAD[i])
        append!(datas_ME, data_ME[i])
        append!(groups, fill(g, len))
        append!(xlabs, fill(x, len))
    end
    p1  = StatPlots.add_label([], boxplot, xlabs, datas_MAD, group=groups, sorted=false, xgrid=false,
                              xtickfontsize=14, ytickfontsize=12, guidefontsize=14);

    fs = 20
    fs_reg = 12
    left, right = 0.3, 2.8
    ss = ((0.0,2.25), (-1.5, 2.2))
    if region==3
        annotate!(left, scaleyes ? ss[1][2]*0.75 : 260, text("A", fs), ylabel="norm. MAD ()",ylims=ss[1], legend=:topleft, title="  Arctic Canada N. (03)") # three outliers ylims=(-10/scale,350/scale)
    elseif region==4
        annotate!(left, scaleyes ? ss[1][2]*0.9 : 320, text("B", fs), ylims=ss[1], legend=false, title="  Arctic Canada S. (04)")
    elseif region==7
        annotate!(left, scaleyes ? ss[1][2]*0.9 :370, text("C", fs), ylims=ss[1], legend=false, title="Svalbard (07)")
    elseif region==11
        annotate!(left, scaleyes ? ss[1][2]*0.9 :130, text("D", fs), ylims=ss[1], legend=false, title="  Central Europe (11)")  # one outlier , ylims=(-10/scale,155/scale)
    end

    p2 = StatPlots.add_label([], boxplot, xlabs, datas_ME, group=groups, sorted=false, xgrid=false,
                             legend=false, xtickfontsize=14, ytickfontsize=12,
                             guidefontsize=14);

    if region==03
        annotate!(left, scaleyes ? ss[2][2]*0.85 :260, text("E", fs), ylabel="norm. ME ()", ylims=ss[2], legend=false) # two outliers, ylims=(-300/scale,300/scale)
    elseif region==04
        annotate!(left, scaleyes ? ss[2][2]*0.85 :300, text("F", fs), ylims=ss[2], legend=false)
    elseif region==7
        annotate!(left, scaleyes ? ss[2][2]*0.85 :350, text("G", fs), ylims=ss[2], legend=false)
    elseif region==11
        annotate!(left, scaleyes ? ss[2][2]*0.85 :100, text("H", fs), ylims=ss[2], legend=false) # one outlier ylims=(-105/scale,140/scale)
    end
    # zero line
    plot!([0.4,2.8],[0,0], lc=:black, lw=1, label="", ls=:dash, alpha=0.5, right_margin=3.0mm, top_margin=3.0mm);

    label_region = if region==7
        "\nSvalbard (07)"
    elseif region==11
        "\nC. Europe (11)"
    elseif region==3
        "\nA. Canada N. (03)"
    elseif region==4
        "\nA. Canada S. (04)"
    end

    #plot!(twinx(p2), [NaN], [NaN], ylabel=label_region, grid=false, showaxis=false, legend=false)

    # # region
    # p3 = plot()
    # t = text("Svalbard (07)", fs_reg)
    # t.font.rotation = pi/2
    # annotate!(right+0.2, 100, t)


    fac=1.0
    p = plot(p1,p2, layout=grid(2,1)) #, size=(fac*300,fac*600));

    saveyes && savefig("figs/select-errors-$(VAWTools.int2str2(region)).png")
    saveyes && savefig("figs/select-errors-$(VAWTools.int2str2(region)).pdf")

    return p
end

if ~isdefined(:data_regions)
    data_regions = [load_validation(rs) for rs in ["03", "04", "07", "11"]]
    proc_regions = [proc_validation(rs, data...) for (rs,data) in zip([3, 4, 7, 11], data_regions)]
end

ps = []
for (r,pr) in zip([3, 4, 7, 11], proc_regions)
    push!(ps, plot_validation(r, pr...; saveyes=false))
end

plot(ps..., layout=grid(1,4), size=(1000, 600))
saveyes && savefig("figs/select-errors-all-regions.png")
saveyes && savefig("figs/select-errors-all-regions.pdf")


for (r,pr) in zip([3, 4, 7, 11], proc_regions)
    plot_supplement(r, pr...)
end
