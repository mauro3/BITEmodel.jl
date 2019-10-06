#########
# use KissMCMC
#########
"""
Numerical parameters for MCMC

```
struct MCMCNum<:APara @deftype Int
    alg=[:met, :emcee, :emceep][2] # which MCMC stepper to use
    niter = 10^5 # number of MC steps
    nthin = 10     # only store theta every so often
    nburnin = niter÷3 # how much burn-in
    nchains = max(min(niter÷1000,100),2) # number of chains used in the emcee sampler
    sigma_ppdf::Vector{F}
    sample_ppdf::Function = theta_vec -> [randn()*sigma_ppdf[i]+theta_vec[i] for i=1:length(theta_vec)]

    # emcee
    drop_low_accept_ratio::Bool = true
    ball_radius::Any = 0.1 # can be a number or vector

    # output
    plotyes::Bool = false
    verbose::Bool = false
end
```
"""
@with_kw struct MCMCNum<:APara @deftype Int
    alg::Symbol=[:metro, :emcee, :emceep][2] # which MCMC stepper to use
    niter = 10^5 # number of MC steps
    nthin = 10     # only store theta every so often
    nburnin = niter÷2 # how much burn-in
    @assert nburnin<niter
    nchains = max(min(niter÷3000,100),11) # number of chains used in the emcee sampler
    sigma_ppdf::Vector{F} = F[] # needs to be length(theta0.th0)
    sample_ppdf::Function = theta_vec -> [randn()*sigma_ppdf[i]+theta_vec[i] for i=1:length(theta_vec)]

    # emcee
    drop_low_accept_ratio::Bool = true
    ball_radius::Any = 0.1 # can be a number or vector
    a_scale::Float64 = 2.0 # jump size, probably needn't be adjusted

    # output
    plotyes::Bool = false
    verbose::Bool = false
    use_progress_meter::Bool = true
end



"""
Results of an MCMC run.

Parameterised on T to allow nothing, when only the prior was sampled.

```
@with_kw struct MCMCSol{T} @deftype T
    thetas::Matrix{Float64}
    reshape_revert::Any
    theta_expect::Vector{Float64}
    theta_mode::Vector{Float64}
    thetas_std::Vector{Float64}
    thetas_minmax::Matrix{Float64}
    blobs::Any
    logdensity_vals::Vector{Float64}

    vols::Union{Void,Vector{Float64}}
    vols_above::Union{Void,Vector{Float64}}
    mc_hs1d
    mc_ivs1d
    mh2d
    sh2d
    miv2d
    siv2d
    min_h2d
    max_h2d
    min_iv2d
    max_iv2d

    # solver performance
    solverstats::Dict{Symbol,Any} # wall_time, eff_sample_size, τ, τ_convergence, accept_ratio, # of samples in thetas

    # inputs
    theta0::Theta0{Float64}
    pmcmc::MCMCNum
    logdensity::Function # density which was sampled

    # Holds various things, such as:
    # misc=Dict(:logposterior=>logposterior,
    #       :logprior=>logprior,
    #       :logprior_debug=>logprior_debug,
    #       :loglikelihood=>loglikelihood,
    #       :fwdm_fn=>fwdm_fn,
    #       :fit_target=>fit_target))
    misc::Dict{Symbol,Any}
end
```
"""
@with_kw struct MCMCSol{T} @deftype T
    thetas::Matrix{Float64}
    reshape_revert::Any
    theta_expect::Vector{Float64}
    theta_mode::Vector{Float64}
    thetas_std::Vector{Float64}
    thetas_minmax::Matrix{Float64}
    blobs::Any
    logdensity_vals::Vector{Float64}

    vols::Union{Void,Vector{Float64}}
    vols_above::Union{Void,Vector{Float64}}
    mc_hs1d
    mc_ivs1d
    mh2d
    sh2d
    miv2d
    siv2d
    min_h2d
    max_h2d
    min_iv2d
    max_iv2d

    # solver performance
    solverstats::Dict{Symbol,Any}  # wall_time, eff_sample_size, τ, τ_convergence, accept_ratio, # of samples in thetas

    # inputs
    theta0::Theta0{Float64}
    pmcmc::MCMCNum
    logdensity::Function # density which was sampled

    # Holds various things, such as:
    # misc=Dict(:logposterior=>logposterior,
    #       :logprior=>logprior,
    #       :logprior_debug=>logprior_debug,
    #       :loglikelihood=>loglikelihood,
    #       :fwdm_fn=>fwdm_fn,
    #       :fit_target=>fit_target))
    misc::Dict{Symbol,Any}
end

function bite(logdensity, theta0::Theta0, pmcmc::MCMCNum;
              resample_logdensity=nothing,
              n_resamples=300,
              misc=Dict{Symbol,Any}())

    sol = mcmc(logdensity, theta0, pmcmc; misc=misc)
    if resample_logdensity!=nothing
        sol = BM.retrace_samples(resample_logdensity, sol, pmcmc, nsamples=n_resamples, misc=sol.misc)
    end

    return sol
end

"""
    mcmc(logposterior, theta0::Theta0, pmcmc::MCMCNum; verbose=true, use_progress_meter=true)

Runs the MCMC.

Input:
- logdensity: log-density to sample
- theta0: one value for metropolis, tuple for emcee (theta0, ball_radius)
KW:
- misc : if passed in it will be stored in the sol.misc dictionary

Output:
- thetas, theta_expect, theta_mode, accept_ratio, blobs, logdensity_value: results
- mc_hs1d, mc_ivs1d : all 1D results
- mh2d, sh2d, miv2d, siv2d : 2D mean and standard dev
- eff_sample_size, τ, τ_convergence: from KissMCMC.eff_samples
- reshape_revert: to undo the chain squash
"""
function mcmc(logdensity, theta0::Theta0, pmcmc::MCMCNum;
              verbose=pmcmc.verbose, use_progress_meter=pmcmc.use_progress_meter,
              misc=Dict{Symbol,Any}())
    tic = time()
    @unpack alg,a_scale = pmcmc
    gl = theta0.gb.gl
    p = pmcmc
    blob_reduce! = if KissMCMC.hasblob(logdensity, theta0.th0)
        _blob_reduce!
    else
        KissMCMC.default_blob_reduce!
    end

    if alg==:emcee
        res_ec = emcee(logdensity, (theta0.th0, p.ball_radius), niter=p.niter,
                       nburnin=p.nburnin, nchains=p.nchains,
                       blob_reduce! = blob_reduce!,
                       use_progress_meter=use_progress_meter,
                       a_scale=a_scale)

        # get metrics
        eff_sample_size, τ, τ_convergence, eff_sample_sizes = eff_samples(res_ec[1])

        thetas, accept_ratio, blobs, logdensity_vals, reshape_revert = squash_chains(res_ec... ,
                                                                                     drop_low_accept_ratio=p.drop_low_accept_ratio,
                                                                                     blob_reduce! = blob_reduce!,
                                                                                     verbose=verbose)
    elseif alg==:metro
        thetas, accept_ratio, blobs, logdensity_vals = metropolis(logdensity, p.sample_ppdf, theta0.th0,
                                                                  niter=p.niter, nthin=p.nthin, nburnin=p.nburnin,
                                                                  blob_reduce! = blob_reduce!)
        eff_sample_size, τ, τ_convergence, nthin, reshape_revert = nothing, nothing, nothing, nothing, nothing
    else
        error("alg unknown.  Note parallel algs may not work currently due to closures.")
    end
    mc_vol, mc_vol_above, mc_hs1d, mc_ivs1d,
    mh2d, sh2d, miv2d, siv2d,
    min_h2d, max_h2d, min_iv2d, max_iv2d #=,hquant, ivquant=# = unpack_blob(blobs, size(gl.dem.v))

    theta_expect = squeeze(mean(thetas,2),2)
    theta_mode = thetas[:,findmax(logdensity_vals)[2]]

    solverstats = Dict{Symbol,Any}()
    nsamples = size(thetas, 2)
    @pack solverstats = eff_sample_size, τ, τ_convergence, eff_sample_sizes, accept_ratio, nsamples

    solverstats[:wall_time] = time() - tic

    return MCMCSol(thetas, reshape_revert,
                   theta_expect, theta_mode, std(thetas,2)[:], hcat([[e...] for e in extrema(thetas,2)]...)',
                   blobs, logdensity_vals,
                   mc_vol, mc_vol_above,
                   mc_hs1d, mc_ivs1d,
                   mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d,
                   solverstats,
                   theta0, pmcmc, logdensity, misc)
end

"""
    retrace_samples(logdensity, sol_in, pmcmc::MCMCNum;
                         nsamples = 200,
                         verbose=true, use_progress_meter=true)

When just sampling the prior, this function will re-sample a certain number of
samples with running the forward model (thus logdensity will tipically be logposterior).
This thus allows to do the predictions
of ice thickness whilst only sampling the priors (which is much cheaper).

Note:
- that of the logdensity, only the blob is used.
- this means also that the statistics, e.g. theta_expect, is calculated from sol_in.theta
  and not from just the thetas used in re-sampling.  This is maybe slightly inconsistent.
"""
function retrace_samples(logdensity, sol_in, pmcmc::MCMCNum;
                         nsamples=200,
                         verbose=pmcmc.verbose, use_progress_meter=pmcmc.use_progress_meter,
                         misc=Dict{Symbol,Any}())

    @unpack theta0 = sol_in
    thetas_in = sol_in.thetas
    p = pmcmc
    demsize = size(theta0.gb.gl.dem.v)

    out = _retrace_samples(logdensity, thetas_in, pmcmc, demsize,
                           nsamples, verbose, use_progress_meter,
                           sol_in.logdensity_vals, sol_in.solverstats)

    return MCMCSol(out..., theta0, pmcmc, logdensity, misc)
end
# Split out logic such that it can be used in retrace_samples_condensed_MCMCSol as well:
function _retrace_samples(logdensity, thetas_orig, pmcmc::MCMCNum, demsize,
                          nsamples, verbose, use_progress_meter,
                          logdensity_vals,  # we want the density of the original pdf not the resampled one
                          solverstats)
    tic = time()

    # Those we want of the original thetas as that is more accurate:
    theta_expect = squeeze(mean(thetas_orig,2),2)
    theta_mode = thetas_orig[:,findmax(logdensity_vals)[2]]

    # Now thin the samples:
    inds = 1:(max(1, size(thetas_orig,2)÷nsamples)):size(thetas_orig,2)
    thetas_in = thetas_orig[:,inds]

    @assert KissMCMC.hasblob(logdensity, thetas_in[:,1])

    thetas, accept_ratio, blobs, logdensity_vals_ = KissMCMC.retrace_samples(logdensity, thetas_in,
                                                                            blob_reduce! = _blob_reduce!,
                                                                            use_progress_meter=use_progress_meter)

    @assert thetas_in == thetas # otherwise something went wrong

    (mc_vol, mc_vol_above, mc_hs1d, mc_ivs1d,
     mh2d, sh2d, miv2d, siv2d,
     min_h2d, max_h2d, min_iv2d, max_iv2d #=,hquant, ivquant=#) = unpack_blob(blobs, demsize)

    solverstats[:wall_time_resample] = time() - tic

    # return as the first arguments to MCMCSol
    return (thetas, nothing,
            theta_expect, theta_mode, std(thetas_orig,2)[:], hcat([[e...] for e in extrema(thetas_orig,2)]...)',
            blobs,
            logdensity_vals[inds], # logdensity_vals needs to be original density
            mc_vol, mc_vol_above,
            mc_hs1d, mc_ivs1d,
            mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d,
            solverstats)
end

############
# Blobs
###########
# see https://github.com/mauro3/KissMCMC.jl/blob/master/examples/notebooks/bayesian-ex-predictions-v2.ipynb

# Two methods needed to intitialize the storage-blob:
function _blob_reduce!(new_blob, niter::Int, nburnin::Int, nthin::Int)
    vol, vol_above, hs1d, ivs1d, hs2d, ivs2d = new_blob
    n = length(hs2d)
    # store 0D results
    vol_out = Array{eltype(hs1d)}((niter-nburnin)÷nthin)
    vol_above_out = Array{eltype(hs1d)}((niter-nburnin)÷nthin)
    # store 1D results
    h1d_out = Array{eltype(hs1d)}(length(hs1d), (niter-nburnin)÷nthin)
    iv1d_out = Array{eltype(ivs1d)}(length(ivs1d), (niter-nburnin)÷nthin)
    # store online-stats: mean, variance, extremas and quantiles
    h2d_out = Series(Group([Mean() for i=1:n]),
                     Group([Variance() for i=1:n]),
                     Group([Extrema() for i=1:n]),
                     #Group([Quantile([0.25,0.5,0.75]) for i=1:n]) # uses too much RAM
                     )
    iv2d_out = Series(Group([Mean() for i=1:n]),
                      Group([Variance() for i=1:n]),
                      Group([Extrema() for i=1:n]),
                      #Group([Quantile([0.25,0.5,0.75]) for i=1:n])
                      )
    return vol_out, vol_above_out, h1d_out, iv1d_out, h2d_out, iv2d_out
end
_blob_reduce!(new_blob, niter::Int, nburnin::Int, nthin::Int, nchains::Int) =
             [_blob_reduce!(new_blob,niter,nburnin,nthin) for i=1:nchains]

# Two methods to update it
function _blob_reduce!(stored_blob, new_blob, ni::Int)
    vol_out, vol_above_out, h1d_out, iv1d_out, h2d_out, iv2d_out = stored_blob
    vol, vol_above, hs1d, ivs1d, hs2d, ivs2d = new_blob
    vol_out[ni] = vol
    vol_above_out[ni] = vol_above
    h1d_out[:,ni] = hs1d
    iv1d_out[:,ni] = ivs1d
    n = length(hs2d)
    fit!(h2d_out, reshape(hs2d, n))
    fit!(iv2d_out, reshape(ivs2d, n))
    nothing
end
_blob_reduce!(stored_blob, new_blob, ni::Int, nc::Int) = _blob_reduce!(stored_blob[nc], new_blob, ni);

# One method to merge several chains
function _blob_reduce!(stored_blobs, chains2keep)
    # chain 1 has one stored_blobs[1], etc.
    vol_out, vol_above_out, h1d_out, iv1d_out, h2d_out, iv2d_out = stored_blobs[chains2keep[1]]
    n = size(h1d_out,1)
    h1d_out = reshape(h1d_out, length(h1d_out))
    iv1d_out = reshape(iv1d_out, length(iv1d_out))
    for i in chains2keep[2:end]
        append!(vol_out, stored_blobs[i][1])
        append!(vol_above_out, stored_blobs[i][2])
        append!(h1d_out, stored_blobs[i][3])
        append!(iv1d_out, stored_blobs[i][4])
        merge!(h2d_out, stored_blobs[i][5])
        merge!(iv2d_out, stored_blobs[i][6])
    end
    h1d_out = reshape(h1d_out, n, length(h1d_out)÷n)
    iv1d_out = reshape(iv1d_out, n, length(iv1d_out)÷n)
    return vol_out, vol_above_out, h1d_out, iv1d_out, h2d_out, iv2d_out
end

"""
Returns
- all samples of hs1d and iv1d
- mean, std, min, and max of the 2D fields
- h2d_stats, iv2d_stats: OnlineStats Series

hs1d, iv1d, mh2d, sh2d, miv2d, siv2d,
min_h2d, max_h2d, min_iv2d, max_iv2d,
# hquant, ivquant,
h2d_stats, iv2d_stats  = unpack_blob(blob, sz)
"""
function unpack_blob(blob, sz)
    vol, vol_above, h1d, iv1d, h2d, iv2d = blob
    hmean, hvar, hext#=, hquant=# = [value.(u) for u in value(h2d)]
    ivmean, ivvar, ivext#=, ivquant=# = [value.(u) for u in value(iv2d)]
    hmin = reshape([h[1] for h in hext], sz)
    hmax = reshape([h[2] for h in hext], sz)
    ivmin = reshape([iv[1] for iv in ivext], sz)
    ivmax = reshape([iv[2] for iv in ivext], sz)
    return tuple(vol, vol_above, h1d, iv1d,
                 map(x->reshape(x,sz), (hmean, sqrt.(hvar), ivmean, sqrt.(ivvar)))...,
                 hmin,hmax,ivmin,ivmax, #hquant, ivquant,
                 h2d, iv2d)
end
# nice:
unpack_blob(blob::Void, sz) = (nothing, nothing,
                               nothing,nothing,nothing,nothing,nothing,nothing,
                               nothing,nothing,nothing,nothing,nothing,nothing)





##########
# Parameters
##########
