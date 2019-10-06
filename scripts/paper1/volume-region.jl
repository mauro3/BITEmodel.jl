# Calculate to total volume of regions

using Plots
include("../mustache-tables.jl")
include("load-one.jl")

# regions = [11]

# # assume independence: gives super small error
# samples = 200
# inds = 1:100
# totvol = zeros(length(inds))
# for out in [out1, out2]
#     for o in out
#         !(o isa Dict) && continue
#         vs = o[:mc_vols_above]
#         for i=inds
#             totvol[i] += vs[i]./1e9 #km3
#         end
#     end
# end
# Plots.histogram(totvol)

# # do not assume independence
# totvol = 0.0;
# totvol_std = 0.0;
# for out in [out1, out2]
#     for o in out
#         !(o isa Dict) && continue
#         # @show o[:vol_above_km3][:mean], std(o[:mc_vols_above])/1e9
#         totvol += o[:vol_above_km3][:mean]
#         totvol_std += std(o[:mc_vols_above])/1e9
#     end
# end
# println("Region $region total vol $(round(Int, totvol)) ± $(round(Int, totvol_std)) km3")

using Measurements, DataStructures

"Values from the G2TI paper"
G2TI = OrderedDict(3 => 28.33±7.35,
                   4 => 8.61±2.23,
                   7 => 7.47±1.94,
                   11 =>0.13±0.03,
#                   14 =>2.87±0.74)  # note, this is for both subregions
                   14 => 2.19±0.56) # this is for subregion 2 only (error is just adjusted from whole)

"""
Total volume in 10^3km^3

The errors (stdev) are just added.  I.e. no independence is assumed
For glaciers which errored it uses the Volume-Area volume instead (if
use_VA=true, default).

Returns:
- total volume with error
- the error-glaciers' V-A volume (should be small)
- total number of glaciers
- number of error glaciers
- total area
- area of error glaciers
"""
function get_vols(o1s, o2s, gnr1, gnr2; use_VA=true)
    out = Dict()
    va_out = Dict()
    nerrs = Dict()
    aerrs = Dict()
    ntot = Dict()
    atot = Dict()
    for reg in keys(o1s)
        glnrs, areas = BM.read_attrib_file(reg)[[1, 5]]
        area = Dict(g=>a/1e6 for (g,a) in zip(glnrs, areas))

        totvol = 0.0;
        totvol_std = 0.0;
        va_vol = 0.0;
        va_vol_std = 0.0;
        nerrs[reg] = 0
        aerrs[reg] = 0.0
        ntot[reg] = length(o1s[reg]) + length(o2s[reg])
        atot[reg] = sum(values(area))
        if reg==14
            # only subregion 2
            atot[reg] = 22871
        end

        for (oo,nr) in zip([o1s[reg], o2s[reg]], [gnr1[reg], gnr2[reg]])
            for (o,n) in zip(oo,nr)
                if !(o isa Dict)
                    a = area[n]
                    if use_VA
                        v = BM.volume_area_mean_h(a*1e6)*a*1e6 /1e9/1000
                        stdv = 0.2*v
                        totvol += v
                        totvol_std += stdv
                        va_vol += v
                        va_vol_std += stdv
                    end
                    nerrs[reg] += 1
                    aerrs[reg] += a
                else
                    # @show o[:vol_above_km3][:mean], std(o[:mc_vols_above])/1e9
                    totvol += o[:vol_above_km3][:mean]/1000
                    totvol_std += std(o[:mc_vols_above])/1e9/1000
                end
            end
        end
        out[reg] = totvol ± totvol_std
        va_out[reg] = va_vol ± va_vol_std
    end
    return out, va_out, ntot, nerrs, atot, aerrs
end

# https://github.com/JuliaPlots/Plots.jl/issues/143
# https://github.com/JuliaPlots/StatsPlots.jl/issues/32

"""
Compares the regional G2TI estimates with ours.
"""
function plot_compare_G2TI(vols; kws...)
    ks = [k for k in keys(G2TI) if haskey(vols, k)]
    vs = hcat([[vols[k].val, G2TI[k].val] for k in ks]...)
    errs = hcat([[vols[k].err, G2TI[k].err] for k in ks]...)
    ind11 = findfirst(ks.==11)
    # p1= StatPlots.groupedbar(map(VAWTools.int2str2, repeat(ks,inner=2)), vs', err=errs', marker=stroke(2, :black),
    #                      xlabel="RGI region", ylabel="V (10^3 km^3)", label="") #group=repeat(["BITE", "G2TI"], inner=length(ks)))

    inds1 = [findfirst(ks.==r) for r in [3,4,7,14] if findfirst(ks.==r)>0]
    inds2 = [findfirst(ks.==r) for r in [11] if findfirst(ks.==r)>0]
    p1= StatPlots.groupedbar(map(VAWTools.int2str2, repeat(ks[inds1],inner=2)), vs[:,inds1]', err=errs[:,inds1]',
                             markercolor=:black, markerstrokecolor=:black, markerstrokewidth=2,
                             xlabel="RGI region", ylabel="Volume (10³ km³)", label="", bar_width=0.7); #group=repeat(["BITE", "G2TI"], inner=length(ks)))

    p0 = StatPlots.groupedbar!([1,2], NaN*vs[:,1:2]', err=NaN*errs[:,1:2]', group=["BITE", "G2TI"]);

    p2 = StatPlots.groupedbar(map(VAWTools.int2str2, repeat(ks[inds2],inner=2)), vs[:,inds2]', err=errs[:,inds2]',
                              markercolor=:black, markerstrokecolor=:black, markerstrokewidth=2,
                              xlabel="", ylabel="", bar_width=0.7, label=""); #, group=repeat(["BITE", "G2TI"], inner=length(ks[inds2])))

    plot(p1, p2, p0, layout=@layout [a{0.8w} c{0.15w} b{0.05w}]; kws...)
#        plot(p1, p2, layout=@layout [a{0.8w} c{0.15w}])
end

function plot_compare_G2TI_2(vols, volsdhdt, atot; kws...)
    ks = [k for k in keys(G2TI) if haskey(vols, k)]
    kks = ["$r ($k)" for (r,k) in zip(["A. Canada N.", "A. Canada S.", "Svalbard", "C. Europe", "Karakoram"], keys(G2TI)) if haskey(vols, k)]
    vs = hcat([[vols[k].val, volsdhdt[k].val, G2TI[k].val] for k in ks]...)
    errs = hcat([[vols[k].err, volsdhdt[k].err, G2TI[k].err] for k in ks]...)
    ind11 = findfirst(ks.==11)

    hmean = 1e6*hcat([[vols[k].val/atot[k] ± vols[k].err/atot[k], volsdhdt[k].val/atot[k] ± volsdhdt[k].err/atot[k], G2TI[k].val/atot[k] ± G2TI[k].err/atot[k]] for k in ks]...)
#    hmean_errs = hcat([[vols[k].err/atot[k], volsdhdt[k].err/atot[k], G2TI[k].err] for k in ks]...)

    # p1= StatPlots.groupedbar(map(VAWTools.int2str2, repeat(ks,inner=2)), vs', err=errs', marker=stroke(2, :black),
    #                      xlabel="RGI region", ylabel="V (10^3 km^3)", label="") #group=repeat(["BITE", "G2TI"], inner=length(ks)))

    inds1 = [findfirst(ks.==r) for r in [3,4,7,14] if findfirst(ks.==r)>0]
    inds2 = [findfirst(ks.==r) for r in [11] if findfirst(ks.==r)>0]
    p1= StatPlots.groupedbar(repeat(kks[inds1],inner=3), vs[:,inds1]', err=errs[:,inds1]',
                             markercolor=:black, markerstrokecolor=:black, markerstrokewidth=1.5,
                             xlabel="RGI region", ylabel="Volume (10³ km³)", label="", bar_width=0.7); #group=repeat(["BITE", "G2TI"], inner=length(ks)))

    off = 0.22
    font = text("").font
    font.rotation = pi/2
    font.pointsize = 10
    for (i,ind) in enumerate(inds1)
        for j=1:3
            isnan(hmean[j,ind]) && continue
            txt = text("$(round(Int, hmean[j,ind]))±$(round(Int, hmean[j,ind].err))m", font)
            annotate!((i-0.5) + off*(j-2), 15, txt)
        end
    end

    #p0 = StatPlots.groupedbar!([1,2,3], NaN*vs[:,1:3]', err=NaN*errs[:,1:3]', group=["BITE", "BITE-dhdt", "G2TI"]);
    p0 = StatPlots.groupedbar([1,2,3], NaN*vs[:,1:3]', err=NaN*errs[:,1:3]', group=["BITE", "BITE-dhdt", "G2TI"], frame=:none);


    p2 = StatPlots.groupedbar(repeat(kks[inds2],inner=3), vs[:,inds2]', err=errs[:,inds2]',
                              markercolor=:black, markerstrokecolor=:black, markerstrokewidth=1.5,
                              xlabel="", ylabel="", bar_width=0.7, label=""); #, group=repeat(["BITE", "G2TI"], inner=length(ks[inds2])))

    for (i,ind) in enumerate(inds2)
        for j=1:3
            isnan(hmean[j,ind]) && continue
            txt = text("$(round(Int, hmean[j,ind]))±$(round(Int, hmean[j,ind].err))m", font)
            annotate!((i-0.5) + off*(j-2), 0.07, txt)
        end
    end

#    plot(p1, p2, layout=@layout [a{0.8w} c{0.15w}]; kws...)
    plot(p1, p0, p2, layout=@layout [a{0.8w} b{0.001w} c{0.15w}]; kws...)
#        plot(p1, p2, layout=@layout [a{0.8w} c{0.15w}])
end


""
function split_camel_case(cc)
    out = []
    for c in cc
        if isupper(c)
            push!(out, ' ')
        end
        push!(out, c)
    end
    out = "$(out...)"
    strip(out)
end

"Total area of regions"
areas = Dict(3=>104920,
             4=>40860,
             7=>33932,
             11=>2091,
             14=>33561)

function table_compare_G2TI(vols)
    ks = [k for k in keys(G2TI) if haskey(vols, k)]

    # names = [:Region, :Volume, :Error, :Volume_G2TI, :Error_G2TI]
    # vals = hcat([Any[k, vols[k].val, vols[k].err, G2TI[k].val, G2TI[k].err] for k in ks]...)'
    # df = DataFrame(vals, names)

    header = ["RGI Region", "N", "A (km²)", "V (10³ km³)", "\$\\bar\{\\mathrm\{h\}\}\$ (m)", "G2TI V (10³ km³)", "G2TI \$\\bar\{\\mathrm\{h\}\}\$ (m)"]
    names = Symbol.(["RGI Region", "N", "A (km²)", "V (10³ km³)", "h", "G2TI V (10³ km³)", "G2TI h (m)"])
    Vals = [Any[replace(replace(split_camel_case(string(BM.reg2dir[k])), "North", "N"), "South", "S") * " ($k)",
                length(get_glaciers(k)[3]),
                areas[k],
                map(x->signif(x.val, 3)±signif(x.err, 3), vols[k]),
#                signif.(vols[k], 3),
                round.(vols[k]*1e6./areas[k]),
                signif.(G2TI[k], 3),
                round.(G2TI[k]*1e6./areas[k])
                ] for k in ks]
    vals = permutedims(hcat(vals...), [2,1])
    return DataFrame(vals, names), header
end

"""
Plot the distribution of parameters
"""
function plot_paras(out1; kws...)
    exp_dist = hcat([o[:thetas_expect] for o in out1 if o isa Dict]...)
    mean_exp = calculate_prior_update(out1, false)[3]
    println(mean_exp)
    pls = []
    for (n,inds) in out1[1][:theta0][:names]
        for i in inds
            yl = i==1 ? BM.getrgi_region(out1[1][:rgi]) : ""
            push!(pls, Plots.histogram(exp_dist[i,:], xlabel=n, label="", ylabel=yl))
            Plots.scatter!([mean_exp[i]], [0], label="")
        end
    end
    Plots.plot(pls...; layout=grid(1,length(mean_exp)), kws...)
end

function plot_paras(o1s::Dict)
    pls = []
    for (k,o) in o1s
        if !isempty(o)
            push!(pls, plot_paras(o))
        end
    end
    Plots.plot(pls..., layout=grid(length(pls),1))
end


## Takes TIME
##       3,4      4,7,11      14
# dirbase = "output/rasters"
# dirs = ["mar-23", "apr-02", "apr-03"]
# o1s, o2s = Dict(), Dict()
# gnr1s, gnr2s = Dict(), Dict()
# for dir in dirs
#     @show dir
#     o1, o2, gnr1, gnr2 = load_all(dir, dirbase)
#     merge!(o1s, o1)
#     merge!(o2s, o2)
#     merge!(gnr1s, gnr1)
#     merge!(gnr2s, gnr2)
# end

# nosigma, withsigma = nothing, nothing
withsigma = load_all("production-without-dhdt", "output/rasters-fit_sigma");#, regions=[11])
withsigma_dhdt_ = load_all("production-with-dhdt", "output/rasters-fit_sigma");#, regions=[11])

# o1s, o2s, gnr1s, gnr2s = [nosigma, withsigma, withsigma_dhdt_][2];


volsdhdt, va_vols, ntot, nerrs, atot, aerrs = get_vols(withsigma_dhdt_...);
vols, va_vols, ntot, nerrs, atot, aerrs = get_vols(withsigma...);

volsdhdt[3] = volsdhdt[7]*NaN
volsdhdt[4] = volsdhdt[7]*NaN
if !haskey(volsdhdt, 14)
    volsdhdt[14] = volsdhdt[7]*NaN
end
# df_vols, header = table_compare_G2TI(vols)
# print(df_to_table(df_vols, header=header, fmt="lrrcrcr"))

display(plot_compare_G2TI_2(vols,volsdhdt,atot))

#display(plot_compare_G2TI(vols))
#savefig("paper1/figs/g2ti_comp.png")
savefig("figs/g2ti_comp.pdf")

error()

## get expectation of region 14
em = []
ss = []
for r in [3,4,7,11]
    if r in [3,4]
        o1s = withsigma[1]
    else
        o1s = withsigma_dhdt_[1]
    end
    @show r
    @show o1s[r][1][:thetas_expect]
    @show o1s[r][1][:theta0][:names]
    if !haskey(o1s, r)
        println("Not using region $r")
        continue
    end
    _,expect_median,_,stdev = calculate_prior_update(o1s[r], false)
    push!(em, expect_median)
    push!(ss, stdev)
end
println(median(hcat(em...),2)[:])
println(mean(hcat(em...),2)[:])
println(std(hcat(ss...),2)[:])
