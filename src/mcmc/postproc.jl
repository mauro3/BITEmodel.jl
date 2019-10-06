# MCMCSol post-proc

"""
$(SIGNATURES)

Return h and iv at the locations where the 1D parameters
are by averaging (weighted by elevation band area).
"""
function h_iv_at_thetas(theta0::Theta0, mc_hs1d, mc_ivs1d; elebins=Float64[])
    gb = theta0.gb
    if elebins==Float64[]
        eles = []
        for k in [:btilde, :temp, :fsl]
            if haskey(theta0.th_ele, k)
                push!(eles, theta0.th_ele[k])
            end
        end
        ele = eles[findmax(length.(eles))[2]]
        elebins = [gb.bands_[1]; ele[1:end-1] + diff(ele)/2; gb.bands_[end]]
    end
    # find indices of bins
    inds = UnitRange{Int}[]
    for (b1,b2) in zip(elebins[1:end-1],elebins[2:end])
        push!(inds, findfirst(gb.ele.<=b1):findlast(gb.ele.>b2))
    end
    @assert [inds...;]==collect(1:length(gb))
    @assert all(length.(inds).>0)

    h_at_th = Array{F}(length(inds),size(mc_ivs1d,2))
    iv_at_th = Array{F}(length(inds),size(mc_ivs1d,2))
    for i=1:size(mc_ivs1d,2)
        for (j,ind) in enumerate(inds)
            h_at_th[j,i] = sum(gb.ws[ind].*gb.ls[ind].*mc_hs1d[ind,i]) /
                sum(gb.ws[ind].*gb.ls[ind])
            iv_at_th[j,i] = sum(gb.ws[ind].*gb.ls[ind].*mc_ivs1d[ind,i]) /
                sum(gb.ws[ind].*gb.ls[ind])
        end
    end
    return h_at_th, iv_at_th
end

"""
$(SIGNATURES)

Return mean h and iv weighted by elevation band area.
"""
function h_iv_mean(mc_hs1d, mc_ivs1d, theta0::Theta0)
    gb = theta0.gb

    h_mean = Vector{F}(size(mc_ivs1d,2))
    iv_mean = Vector{F}(size(mc_ivs1d,2))
    for i=1:size(mc_ivs1d,2)
        h_mean[i] = sum(gb.ws.*gb.ls.*mc_hs1d[:,i]) /
            sum(gb.ws.*gb.ls)
        iv_mean[i] = sum(gb.ws.*gb.ls.*mc_ivs1d[:,i]) /
            sum(gb.ws.*gb.ls)
    end
    return h_mean, iv_mean
end

"""
$(SIGNATURES)

Condenses the MCMCSol into something manageable.

If store2d==true, then also store 2D fields (uses lots of space)

Uses only basic Julia types for storage to enable easy writing to a
file and future-proofing it.
"""
function condense_MCMCSol(sol::MCMCSol{T}; store2d=false, store_thetas=false) where T
#    error("add h error when there is radar")

    @unpack bands_, gl = sol.theta0.gb
    @unpack x, y = gl.dem
    dx = step(x)
    out = Dict{Symbol,Any}()

    # glacier
    out[:rgi] = getrgi(gl.gid)
    out[:area_km2] = area(gl) / 1e6
    out[:volarea_km3] = volume_area_mean_h(gl) * area(gl) /1e9
    out[:proj] = gl.dem.proj
    out[:fit_target] = string(sol.misc[:fit_target])

    # collect the ranges for future-proofing:
    out[:x] = collect(gl.dem.x)
    out[:y] = collect(gl.dem.y)
    out[:bands] = collect(bands_)

    if !(T<:Void)
        # 0D
        out[:vol_km3] = type2dict(StatsBase.summarystats(sol.vols./1e9))
        out[:vol_above_km3] = type2dict(StatsBase.summarystats(sol.vols_above./1e9))

        # reduce number of samples of vols
        nsamples = length(sol.vols)
        thin = if nsamples==sol.solverstats[:nsamples]
            max(sol.solverstats[:τ],1)
        else
            1
        end
        inds = thin==-1 ? (1:nsamples) : (1:thin:nsamples)
        # we want at most between 250 and 500 samples
        if length(inds)>250
            inds = inds[1:(length(inds)÷150):end]
        end
        out[:mc_vols] = sol.vols[inds]
        out[:mc_vols_above] = sol.vols_above[inds]

        # 1D vars
        out[:hs1d] = mean(sol.mc_hs1d,2)
        out[:hs1d_std] = std(sol.mc_hs1d,2)
        out[:hs1d_min] = minimum(sol.mc_hs1d,2)
        out[:hs1d_max] = maximum(sol.mc_hs1d,2)

        out[:ivs1d] = mean(sol.mc_ivs1d,2)
        out[:ivs1d_std] = std(sol.mc_ivs1d,2)
        out[:ivs1d_min] = minimum(sol.mc_ivs1d,2)
        out[:ivs1d_max] = maximum(sol.mc_ivs1d,2)

        # 2D vars
        if store2d
            out[:hs] = sol.mh2d
            out[:hs_std] = sol.sh2d
            out[:hs_min] = sol.min_h2d
            out[:hs_max] = sol.max_h2d

            out[:ivs] = sol.miv2d
            out[:ivs_std] = sol.siv2d
            out[:ivs_min] = sol.min_iv2d
            out[:ivs_max] = sol.max_iv2d
        end
    end

    # thetas
    if store_thetas
        out[:thetas] = sol.thetas
        out[:reshape_revert] = sol.reshape_revert
    end
    out[:thetas_expect] = sol.theta_expect
    out[:thetas_std] = sol.thetas_std
    out[:thetas_minmax] = sol.thetas_minmax
    out[:thetas_mode] = sol.theta_mode
    out[:thetas_quartiles] = [quantile(sol.thetas[i,:], [0.25, 0.5, 0.75]) for i=1:size(sol.thetas,1)]
    out[:thetas_quant_95] = [quantile(sol.thetas[i,:], [0.05, 0.95]) for i=1:size(sol.thetas,1)] # as in ITMIX paper plots

    # others
    out[:theta0] = Dict(:th0 => sol.theta0.th0,
                        :names => [k=>collect(v) for (k,v) in sol.theta0.names],
                        :th_ele => [k=>collect(v) for (k,v) in sol.theta0.th_ele])

    out[:solverstats] = sol.solverstats
    out[:pmcmc] = Parameters.type2dict(sol.pmcmc)
    out[:pm] = Parameters.type2dict(sol.theta0.pm)
    out[:pn] = Parameters.type2dict(sol.theta0.pn)
    out[:pp] = Parameters.type2dict(sol.theta0.pp)
    pop!(out[:pmcmc], :sample_ppdf)

    # errors statistics
    ae, _, err = compare_to_radar(sol)
    out[:abserr_h] = type2dict(ae[1])
    out[:err_h] = type2dict(err[1])
    # TODO iv errors

    return out
end

function retrace_samples_condensed_MCMCSol(o::Dict)
    (gid, gl, gb, pp, pm, pn, pl, th0d,
     theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) = init(o[:rgi])

    _retrace_samples(logposterior, o[:thetas], pmcmc::MCMCNum, demsize,
                          nsamples, verbose, use_progress_meter,
                          logdensity_vals,  # we want the density of the original pdf not the resampled one
                          solverstats)
end


"""
$(SIGNATURES)

Makes
"""
function make_geotiff(sol::MCMCSol, fln)

    # make gridded
    @unpack bands_, gl = sol.theta0.gb
    @unpack x, y, proj = gl.dem
    @unpack glaciermask = gl
    dx = step(x)
    z = Matrix{Float32}(0,0)

    exp = Gridded{Float32}(x, y, sol.mh2d, z, true, proj)
    std = Gridded{Float32}(x, y, sol.sh2d, z, true, proj)
    min = Gridded{Float32}(x, y, sol.min_h2d, z, true, proj)
    max = Gridded{Float32}(x, y, sol.max_h2d, z, true, proj)

    VAWTools.write_geotiff([exp, std, min, max], fln, nodataval=0.0)
end
