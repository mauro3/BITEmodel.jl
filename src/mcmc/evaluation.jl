import StatsBase
square_error2RMSE(sqrerr, n) = sqrt(sqrerr/n)
abs_error2MAE(abserr, n) = sqrerr/n

"""
$(SIGNATURES)

Returns tuples (control-lines, fitting-lines, all-lines):
- statistics of absolute errors (StatsBase.summarystats)
- statistics of square errors
- cumulative h error
- TODO ratio of track-length with errors smaller than:
  - 20%
  - std of h
  - 2std of h
"""
function compare_to_radar(sol::MCMCSol, radarlines_fitting::AbstractVector{<:Integer}=Int[],
                          radarlines_control=nothing; # nothing means all but the radarlines_fitting
                          thin=1, meters_thin=0.0)
    compare_to_radar(sol.theta0.gb.gl, sol.mh2d, sol.sh2d, radarlines_fitting, radarlines_control;
                     thin=thin, meters_thin=meters_thin)
end
function compare_to_radar(gl, hs2d, hs2ds, radarlines_fitting::AbstractVector{<:Integer}=Int[],
                          radarlines_control=nothing; # nothing means all but the radarlines_fitting
                          thin=1, meters_thin=0.0)
    hs2d = hs2d isa Gridded ? hs2d.v : hs2d
    hs2ds = hs2ds isa Gridded ? hs2ds.v : hs2ds
    if gl.h==SomeData(nothing)
        return NaN, NaN, NaN, (), (), ()
    end
    radarlines_all = 1:length(gl.h.vals.splits)
    if radarlines_control==nothing
        radarlines_control = [i for i=1:length(gl.h.vals.splits) if !(i in radarlines_fitting)]
    end
    rls = (radarlines_control, radarlines_fitting, radarlines_all)
#     RMSE = map(x->square_error2RMSE(x...), map(r->thick_err(gl, hs2d, r, thin, sqrerr), rls))
#     max_err = map(r->thick_err(gl, hs2d, r, thin, abserr, +)[1], rls)

# @show    StatsBase.summarystats(thick_err(gl, hs2d, rls[1], thin, abserr, push!, Float64[])[1])
#     @show map(r->isnan(r) StatsBase.summarystats(thick_err(gl, hs2d, r, thin, abserr, push!, Float64[])[1]), rls)
#     MAE = map(x->x[1]/x[2], map(r->thick_err(gl, hs2d, r, thin, abserr, +), rls))
#     quantiles = ()

    SEs = []
    AEs = []
    errs = []
    for r in rls
        err, n = thick_err(gl, hs2d, r, thin, -, push!, Float64[])

        if isequal(err, NaN) || isempty(err)
            push!(AEs, NaN)
            push!(SEs, NaN)
            push!(errs, NaN)
        else
            push!(AEs, StatsBase.summarystats(abs.(err)))
            push!(SEs, StatsBase.summarystats(err.^2))
            push!(errs, StatsBase.summarystats(err))
        end
    end

    # # errors as Traj:
    # err = error_h(hs2d, gl, radarlines_fitting, in_glacier(gl))

    # smaller than
    smaller_20 = () # map(r->/(thick_err(gl, hs2d, r, thin,
                    #              (model,data) -> abs(model-data)/abs(data)<0.2, +)...), rls)
    smaller_std = ()
    smaller_2std = ()

    # remember the other `return`
    return AEs, SEs, errs, smaller_20, smaller_std, smaller_2std
end

"""
$(SIGNATURES)

"""
function compare_to_iv(gl, iv2d, iv2ds, fwdsol)
    error("TODO: Update to same API as compare_to_radar")
    iiv  = Interp.interpolate((dem.x, dem.y), ivs2d, Interp.Gridded(Interp.Linear()) )
    iivs  = Interp.interpolate((dem.x, dem.y), ivs2ds, Interp.Gridded(Interp.Linear()) )

    thin = 1
    if gl.iv.vals isa Gridded{F}
        ivmask = gl.ivmask .& fwdsol.ivs2d_band_mask
        RMSE = square_error2RMSE(iv2d_err(gl, iv2d, ivmask, thin)...)
        max_err = iv2d_err(gl, iv2d, ivmask, thin, abserr, max)
        # produces garbage, would need some abstol:
        # max_rel_err = iv2d_err(gl, iv2d, ivmask, thin, rel_abserr, max)

    elseif gl.iv.vals isa Traj{F}


    end
    quantiles = ()
    return RMSE, max_err, MAE, quantiles, smaller_20, smaller_std, smaller_2std
end


"""
    summarize_mcmc_sols(sols::Vector{<:MCMCSol})

Returns the quartiles of the quartiles of the absolute error
of the solutions ins `sol`.

Return: q1, median, q3

TODO:
- allow using particular fitting lines?
"""
function summarize_mcmc_sols(sols::Vector{<:MCMCSol})
    # do stats on quartiles
    q1 = Float64[]
    q2 = Float64[]
    q3 = Float64[]
    for sol in sols
        ae = compare_to_radar(sol)[1][1]
        push!(q1, ae.q25)
        push!(q2, ae.median)
        push!(q3, ae.q75)
    end
    # now summarize the summaries

    return StatsBase.summarystats.((q1, q2, q3))
end
function summarize_mcmc_sols(gls::Vector, hs2d::Vector, hs2ds::Vector)
    # do stats on quartiles
    q1 = Float64[]
    q2 = Float64[]
    q3 = Float64[]
    for (gl, h, hs) in zip(gls, hs2d, hs2ds)
        ae = compare_to_radar(gl, h, hs)[1][1]
        push!(q1, ae.q25)
        push!(q2, ae.median)
        push!(q3, ae.q75)
    end
    # now summarize the summaries

    return StatsBase.summarystats.((q1, q2, q3))
end
