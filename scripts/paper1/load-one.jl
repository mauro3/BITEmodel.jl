include("../rgi-setup.jl")
using BSON

function load_one(region, dir, dirbase="output/rasters-fit_sigma")
    base = dirname(@__FILE__)
    region_str = VAWTools.int2str2(region)
    dir_ = joinpath(base, "..", dirbase, dir)
    loaded1 = BSON.load(dir_*"/$(region_str)_out1.bson");
    loaded2 = BSON.load(dir_*"/$(region_str)_out2.bson");
    out1 = loaded1[:out1];
    out2 = loaded2[:out2];
    glnrs1, glnrs2 = try
        (loaded1[:glnrs1],
         loaded1[:glnrs2])
    catch
        (nothing, nothing)
    end
    # there was a mistake in glnrs1 for region=3
    if region==3
        for (i,o) in enumerate(out1)
            if (o isa Dict) && (o[:rgi]=="RGI60-03.00452" || o[:rgi]=="RGI60-03.00410")
                out1[i] = ErrorException("no h-data")
            end
        end
    end
    return out1, out2, glnrs1, glnrs2
end

function load_all(dir, dirbase="output/rasters-fit_sigma", regions=[3,4,7,11,14])
    out1s, out2s = Dict(), Dict()
    gnr1, gnr2 = Dict(), Dict()
    for r in regions
        print("Loading region $r...")
        try
            o1,o2,g1,g2 = load_one(r, dir, dirbase)
            out1s[r] = o1
            out2s[r] = o2
            gnr1[r] = g1
            gnr2[r] = g2
            println("done.")
        catch e
            println("error: $e.")
        end
    end
    return out1s, out2s, gnr1, gnr2
end


nothing
