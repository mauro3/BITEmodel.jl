# plots a map of the radarlines of all glacier with radar

using VAWTools, ProgressMeter
import BITEModel
const BM = BITEModel
# include("../rgi-setup.jl")

region_strs = ["03", "04", "07", "11", "14"]

for region_str in region_strs
    @show region_str

    region = parse(Int, region_str)

    glnrs1 = BM.get_glaciers_with_radar(region)
    @showprogress 0.5 "Plotting..." for g in glnrs1
        try
            gid = BM.RGIGlacier(g, region)
            pl = BM.LoadPara(gid)
            gl,pp,pm,pn,pl = BM.load_glacier(gid,pl)
            # BM.plot_radar2d(gl)
            # Plots.savefig("paper1/figs/radar/$(region_str)_$(VAWTools.int2str5(g)).png")

            Plots.plot(gl)
            Plots.savefig("paper1/figs/dem/$(region_str)_$(VAWTools.int2str5(g)).png")
        catch e
            if e isa InterruptException
                rethrow(e)
            end
            println("$g Error $e")
        end
    end
end


# get sorted
function get_sorted(reg)
    out = Int[]
    for f in readdir("paper1/figs/radar/sorted/$(VAWTools.int2str2(reg))")
        push!(out, parse(Int, split(splitext(f)[1], '_')[2]))
    end
    out
end
