# GlaThiDa ice thickness loader
#
# This is a bit of a mess as it is hard to find the right glacier.
# Strategy is to just use Matthias' compiled data.

"""
Returns thickness data if there is any.  Uses Matthias's compiled data instead of directly GlaThiDa.
"""
function load!(ds::DataSet{ThicknessData, GlaThiDaLoader}, dem,
               glaciermask, target_proj, outline::Outline{RGIGlacier}, cropbox,
               pp::Phys)
    if isfile(ds.files[1])
        # get radar lines
        out, head = readdlm(ds.files[1], F, header=true)
        x = out[:,2]
        y = out[:,3]
        z = out[:,4]
        z[z.==-9999] = FILL
        bed = out[:,5]
        bed[bed.==-9999] = FILL
        thick = out[:,6]
        thick[thick.==-9999] = FILL
        if all(z.==FILL) && all(thick.!=FILL) && all(bed.!=FILL)
            z .= thick.+bed
        elseif all(thick.==FILL)
            thick .= z.-bed
        end
        keep = Int[]
        for (i,t) in enumerate(thick)
            t!=FILL && push!(keep,i)
        end

        # coordinate transform
        for (i,(xx,yy,zz)) in enumerate(zip(x,y,z))
            x[i],y[i],z[i] = transform_proj([xx,yy,zz], "+proj=longlat +datum=WGS84" , target_proj)
        end
        tr = Traj(x=x[keep],y=y[keep],v=thick[keep],err=z[keep])

        # some radar data has all zero, if so discard
        if all(abs.(tr.v).<1)
            println("All thickness data less than 1m, not using GlaThiDa measurements.")
            return SomeData(nothing)
        end

        # Sort it.  Do this a few times, makes things better
        if length(tr)>3
            for i = 1:3
                if !VAWTools.issorted_traj(tr, 1:min(length(tr), 15)) && length(tr)<2*10^4
                    tr = VAWTools.sort_traj(tr)
                end

                # purge points closer than some distance apart
                len = sqrt(outline.gid.attrs.area)
                window = len/100
                tr = VAWTools.downsample(tr, window)
            end
        end

        SomeData(vals=tr, sigma=ds.opts[:sigma])
    elseif isfile(ds.files[2])
        return SomeData(nothing)

        # # Get a mean thickness in AVG_*.dat
        # out, head = readdlm(ds.files[2], F, header=true)
        # rginr = out[:,3]
        # @assert length(unique(rginr))==length(rginr)
        # i = findfirst(rginr.==getrgi_nr(outline.gid))
        # if i==0
        #     return SomeData(nothing)
        # else
        #     meanh = out[i,end]
        #     return SomeData(vals=meanh, sigma=ds.opts[:sigma_mean_h_rel]*meanh)
        # end
    else
        return SomeData(nothing)
    end
end


# ##
# ttt = CSV.read(GZip.open("/home/mauro/julia/mypkgs/BITEModel/scripts/data/GlaThiDa_2016/TTT_2016_corr.csv.gz"), datarow=5,  rows_for_type_detect=10^6, header=3)
# pu = by(ttt, :GlaThiDa_ID, d->DataFrame(c=first(d[:POLITICAL_UNIT]),
#                                    lat=mean(d[:POINT_LAT]),
#                                    lon=mean(d[:POINT_LON])
#                                    ))
# pu[pu[:c].=="NO",:] # two in Svalbard


# t = CSV.read("/home/mauro/julia/mypkgs/BITEModel/scripts/data/GlaThiDa_2016/T_2016_corr.csv", datarow=3,  rows_for_type_detect=10^6, header=1)
# no = t[t[:POLITICAL_UNIT].=="NO",:]
# no = no[no[:LAT].>75,:]
