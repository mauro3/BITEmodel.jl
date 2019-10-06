export postproc1d, postproc2d

function postproc1d!(fwdsol::FWDSol)
    @unpack gb, pp, pm, btilde, fsl = fwdsol
    qtot1d, qd1d = postproc1d(gb, pp, pm,
                              btilde, fsl)
    fwdsol.misc[:qtot1d] = qtot1d
    fwdsol.misc[:qd1d] = qd1d
    nothing
end
"1D post-processing"
function postproc1d(gb, pp, pm,
                    btilde1d=gb.bdot-gb.dhdt,
                    fsl=gb.fsl
                    )
    qtot1d = flux_n(gb, btilde1d)
    qd1d = flux_deform(qtot1d, on_nodes(fsl), pp.r)
    return qtot1d, qd1d
end

"2D post-processing"
function postproc2d()
end

"""
2D errors of a field vs a measurement field of gl (:h, :iv)

See also iv2d_err.
"""
function error2d(h_or_iv::Matrix, gl, field, mask=nothing; fillval=NaN)
    obs_ = getfield(gl, field).vals
    samesize = size(obs_.v)==size(h_or_iv)
    obs  = samesize ? obs_ :
        Interp.interpolate((obs_.x, obs_.y), obs_.v, Interp.Gridded(Interp.Linear()) )[gl.dem.x, gl.dem.y]
    out = obs.v .- h_or_iv
    if mask==nothing
        return out
    else
        out[.!mask] = fillval
        return out
    end
end

"""
    error_h(hs2d::Matrix, gl)

Error along the radar trajectories: model - measurement

Returns the trajectories with the errors on them:

err, bed_model, surf, radar_bed, thick_model

See also compare_to_radar, `thick_err`.
"""
function error_h(hs2d::Matrix, gl, rlines_to_drop=[], in_glacier=in_glacier(gl))
    tr = gl.h.vals
    hi = Interp.interpolate((gl.dem.x, gl.dem.y), hs2d, Interp.Gridded(Interp.Linear()) )
    si = Interp.interpolate((gl.dem.x, gl.dem.y), gl.dem.v, Interp.Gridded(Interp.Linear()) )
    err = deepcopy(gl.h.vals)
    bed = deepcopy(gl.h.vals)
    radar_bed = deepcopy(gl.h.vals)
    surf = deepcopy(gl.h.vals)
    thick = deepcopy(gl.h.vals)
    for i in eachindex(tr.v)
        x,y = tr.x[i], tr.y[i]
        if in_glacier[x,y]==1 # only inside glacier
            err_ = hi[x,y] - tr.v[i]
            err.v[i] = err_
            err.err[i] = 0 # hi[tr.x[i],tr.y[i]] # if the interpolation values are also desired
            bed.v[i] = si[x,y] - hi[x,y]
            bed.err[i] = 0
            surf.v[i] = si[x,y]
            surf.err[i] = 0
            radar_bed.v[i] = si[x,y] - tr.v[i]
            radar_bed.err[i] = 0
            thick.v[i] = hi[x,y]
            thick.err[i] = 0
        else
            # otherwise drop the point
            err.v[i] = NaN
            err.err[i] = NaN
            bed.v[i] = NaN
            bed.err[i] = NaN
            surf.v[i] = NaN
            surf.err[i] = NaN
            radar_bed.v[i] = NaN
            radar_bed.err[i] = NaN
            thick.v[i] = NaN
            thick.err[i] = NaN
        end
    end
    return err, bed, surf, radar_bed, thick
end
