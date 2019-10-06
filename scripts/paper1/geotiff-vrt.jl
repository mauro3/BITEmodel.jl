# Makes VRT of all the geotiffs in a folder and in a region.  One vrt per UTM zone.

using VAWTools
import ArchGDAL
const AG = ArchGDAL

regions = [11]
dir = "../output/rasters-fit_sigma/production-with-dhdt/"

for region in [3,4,7,11,14][end:end]
    @show region
    fls = Dict() # UTM zones => list of files
    for f in readdir(dir)
        ff, ext = splitext(f)
        ext!=".tif" && continue
        reg, nr = parse.(Int, split(split(ff, '-')[end], '.') )
        reg!=region && continue
        # get all projections
        #p = readstring(`gdalsrsinfo -o proj4 $(joinpath(dir,f))`)
        # Faster with ArchGDAL?
        p = AG.registerdrivers() do
            AG.read(joinpath(dir,f)) do dataset
                proj4 = strip(AG.toPROJ4(AG.importWKT(AG.getproj(dataset))))
            end
        end
        zone = [[nr,parse(Int, split(s,'=')[2])] for s in split(p) if split(s, '=')[1]=="+zone"][1]
        push!(get!(fls, zone[2], []), joinpath(dir,f))
    end

    for (z,fs) in fls
        @show z
        @show vrt = joinpath(dir, "composite-$(VAWTools.int2str2(region))-zone_$z.vrt")
        run(`gdalbuildvrt $vrt $(fs...)`)
    end
end
