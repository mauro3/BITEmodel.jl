# Loads SMB results from Matthias' GloGEM model
#
# Format is Elevation bands vs year total-balance & winter-balance
# Elev  1980  1981  1982  1983  1984  1985  1986  1987  1988  1989  1990  1991  1992  1993  1994  1995  1996  1997  1998  1999  2000  2001  2002  2003  2004  2005  2006  2007  2008  2009  2010  2011  2012  1980  1981  1982  1983  1984  1985  1986  1987  1988  1989  1990  1991  1992  1993  1994  1995  1996  1997  1998  1999  2000  2001  2002  2003  2004  2005  2006  2007  2008  2009  2010  2011  2012
# 2495  -0.23  -0.55  -0.79  -0.61  -0.21  -0.90  -1.19  -0.82  -1.17  -0.79  -1.20  -0.87  -1.45  -0.91  -1.30  -0.35  -1.18  -0.97  -1.44  -1.30  -1.41  -0.28  -1.72  -2.36  -0.59  -1.93  -1.43  -2.05  -1.33  -1.73  -0.98  -2.06  -1.68   0.72   0.71   0.76   0.72   0.55   0.53   0.68   0.50   0.64   0.44   0.53   0.52   0.70   0.70   0.71   0.65   0.06   0.62   0.47   0.78   0.65   0.93   0.32   0.62   0.71   0.23   0.47   0.11   0.56   0.59   0.53   0.48   0.54
# 2505  -0.21  -0.52  -0.76  -0.59  -0.18  -0.87  -1.16  -0.80  -1.13  -0.76  -1.17  -0.85  -1.42  -0.88  -1.27  -0.32  -1.14  -0.94  -1.41  -0.91  -1.38  -0.26  -1.68  -2.33  -0.57  -1.52  -1.40  -2.01  -1.30  -1.70  -0.95  -1.73  -1.65   0.72   0.71   0.76   0.73   0.56   0.54   0.69   0.51   0.65   0.45   0.53   0.53   0.70   0.71   0.72   0.66   0.07   0.63   0.48   0.78   0.65   0.94   0.32   0.62   0.72   0.24   0.47   0.12   0.56   0.59   0.53   0.48   0.54
# ...


dir2reg  = Dict(
    :Alaska=>1,
    :WesternCanada=>2,
    :ArcticCanadaNorth=>3,
    :ArcticCanadaSouth=>4,
    :Greenland=>5,
    :Iceland=>6,
    :Svalbard=>7, # and Jan Mayen
    :Scandinavia=>8,
    :RussianArctic=>9,
    :NorthAsia=>10,
    :CentralEurope=>11,
    :CaucasusMiddleEast=>12,
    :CentralAsia=>13,
    :SouthAsiaWest=>14,
    :SouthAsiaEast=>15,
    :LowLatitudes=>16,
    :SouthernAndes=>17,
    :NewZealand=>18,
    :SubAntarctic=>19
)

reg2dir = Dict(r=>n for (n,r) in dir2reg)

function load!(ds::DataSet{BdotData, GloGEMLoader}, dem, glaciermask, target_proj,
               outline, cropbox, pp::Phys)
    @unpack files, opts, daterange = ds
    @unpack sigma = opts
    make2D = get(opts, :make2D, true)

    rgi = getrgi(outline.gid)
    vers,reg,nr = parse_rgi(rgi)
    if vers!=v"6"
        error("GloGEM runs on RGI version 6 but here it is version $vers")
    end
    # make filename
    #base = "/home/mauro/netfs/vaw/root/scratch_net/sermeq/mhuss/r6_global_results"
    #ff = joinpath(base, string(reg2dir[reg]), "PAST/mb_elevation", "belev_$(VAWTools.int2str5(nr)).dat")
    ff = ds.files[1]
    #
    out = get_cache_ram!(ds, :bdot, ()->readdlm(ff, F, header=true))
    # only keep first half as second half is winter-balance:
    years = map(x->parse(Int,x),out[2][2:end÷2])
    bdots = out[1][:,2:end÷2] # values at elevation bands and years
    for i in eachindex(bdots)
        if bdots[i]==-99
            bdots[i]=FILL
        end
    end

    ele_ = out[1][:,1] # middle of elevation band
    # convert to range
    st = diff(ele_)[1]
    ele = ele_[1]:st:ele_[end]
    @assert collect(ele)==ele_

    # average over some years
    if daterange==(Date(),Date()) # use all
        bdot = squeeze(mean(bdots,2),2)
    else
        inds = Dates.year(daterange[1]).<=years.<=Dates.year(daterange[2])
        bdot = squeeze(mean(bdots[:,inds],2),2)
    end

    # scale to m ice / a
    scale_m_ice!(bdot, pp)
    vals = Gridded1d(ele,bdot)

    # fill gaps
    fillgaps!(vals, 3)

    if make2D
        # make a 2D bdot field:
        bands, bandi = VAWTools.bin_grid(dem, vals.x, glaciermask) #, min_bands=0)
        bdot_2D = VAWTools.map_back_to_2D(size(dem.v), bandi, vals.v)
        return SomeData(vals=Gridded(dem.x, dem.y, bdot_2D), sigma=sigma)
    else
        return SomeData(vals=vals, sigma=sigma)
    end
end
