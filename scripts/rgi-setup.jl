# The glaciers for ESA:

println("Time to load BITEModel.jl")
@time using BITEModel
const BM=BITEModel
using Parameters
import PyPlot
using Base.Test


"""
    get_glaciers(region)

Return the glacier numbers to be run.  Sorted with
large glaciers first.

Return:
- glnrs1 -- with radar
- glnrs2 -- without radar
- glnrs_all -- both combined.
"""
function get_glaciers(region)

    glnrs, lon, lat, subregion, area  = BM.read_attrib_file(region)
    # purge some
    inds = 1:length(glnrs)
    # if region==7
    #     inds = 1:length(glnrs)
    #     # inds = (.!( (lon.>18.5) .& (lat.>79.7) ) .& # rm Westfonna
    #     #         .!( (lon.>20) .& (lat.>79.2) )) # rm Austfonna
    # elseif region==11
    #     # Alps only
    #     inds = subregion.==1  # only Alps: 3892 glaciers
    if region==14
        inds = subregion.==2  # only Karakoram: 13757 glaciers
    end
    # else
    #     inds = 1:length(glnrs)
    # end;
    glnrs_all = glnrs[inds]

    # now get radar vs non-radar
    glnrs1 = BM.get_glaciers_with_radar(region)
    # purge glnrs1 which are not in glnrs_all
    glnrs1 = [n for n in glnrs1 if (n in glnrs_all)]
    glnrs2 = [n for n in glnrs_all if !(n in glnrs1)]

    # sort with big glaciers first
    area = BM.read_attrib_file(region, glnrs1)[5]
    glnrs1 = glnrs1[reverse(sortperm(area))]

    area = BM.read_attrib_file(region, glnrs2)[5]
    glnrs2 = glnrs2[reverse(sortperm(area))]

    area = BM.read_attrib_file(region, glnrs_all)[5]
    glnrs_all = glnrs_all[reverse(sortperm(area))]

    return glnrs1, glnrs2, glnrs_all
end

function fit_it(nr, region,
                sigma_of_model, fit_vars, fit_target, runtyp, posterior_2D_or_1D::Val{I},
                expect=nothing, stdev=nothing;
                store_thetas=false,
                update_cache=true, prior_only_resamples=300,
                retsol=false,
                dir_output=error("Provide dir_output as kwarg"),
                save_geotiff=true,
                save_sol=false,
                flroot="",
                fit_sigma=false,
                debugging=true,
                plkws...
                ) where I
    (gid, gl,gb,pp,pm,pn,pl,th0d,
     theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) = BM.init(nr, region, fit_vars, fit_target,
                                           sigma_of_model, runtyp;
                                           fit_sigma=fit_sigma,
                                           update_cache=update_cache,
                                           plkws...)

    # Update prior (if needed) to take values of expect and stdev into account.
    # This assumes that the previous posteriors had independent normal distributions.
    if expect!=nothing
        @assert length(expect)==length(theta0.th0)
        logprior_ = logprior
        # update priors to take expect and stdev into account
        logprior = theta -> logprior_(theta) +
            sum( - (theta[i]-expect[i])^2/(2*stdev[i]^2) for i=1:length(expect))
        # update the logposterior and logposterior1d to take the new logprior into accout:
        logposterior, logposterior1d = make_logposterior(loglikelihood, loglikelihood1d, logprior)
    end
    if I==2 && !(runtyp=="test" || runtyp=="testmid")
        merge!(pmcmc_defaults, Dict(:niter=>2*10^4,
                                    :nthin =>10,
                                    :nburnin=>5*10^3,
                                    :nchains=>50))
    end

    tosample = (logposterior, logposterior1d, logprior)[I]
    pmcmc = BM.MCMCNum(;pmcmc_defaults...)

        # some tests:
    if debugging
        fwdm_fn(theta0.th0);
        print("    Time to run forward model:")
        @time fwdsol = fwdm_fn(theta0.th0);
        print("    Posterior value:")
        @show logposterior(theta0.th0)[1]
        print("    Time to calculate posterior:")
        @time logposterior(theta0.th0)[1]
        print("    Starting MCMC:")
    end

    if logprior(theta0.th0)==-Inf
        logprior_debug(theta0.th0)
        error("Zero density of prior of $gl at $(theta0.th0)")
    end

    sol = BM.mcmc(tosample, theta0, pmcmc;
                  verbose=debugging,
                  use_progress_meter=progress_bar,
                  misc=Dict(:logposterior=>logposterior,
                            :logposterior1d=>logposterior1d,
                            :logprior=>logprior,
                            :logprior_debug=>logprior_debug,
                            :loglikelihood=>loglikelihood,
                            :loglikelihood1d=>loglikelihood1d,
                            :fwdm_fn=>fwdm_fn,
                            :fwdm1d_fn=>fwdm1d_fn,
                            :fit_target=>fit_target,
                            :runtyp=>runtyp,
                            :posterior_or_prior=>I,
                            :expect=>expect,
                            :stdev=>stdev,
                            # pp, pm, pn, etc are part of FWDSol
                            ))

    # if we sampled the logprior or logposterior1d, then run the forward model for
    # prior_only_resamples to get ice thickness samples
    if I==2
        tosample = logposterior
        sol.misc[:fit_target] = BM.FitTarget.prior
        sol = BM.retrace_samples(tosample, sol, pmcmc, nsamples=prior_only_resamples, misc=sol.misc)
    end

    # postproc


    debugging && println("MCMC of experiment $exp of $(gl.gid) done in $(signif(sol.solverstats[:wall_time]/60,4)) minutes")
    debugging && println("   sampling $(["log-posterior", "log-prior"][I])")
    #debugging && I==1 && println("   Rhat=$(mean(sol.solverstats[:Rhat])), $(maximum(sol.solverstats[:Rhat]))")
    debugging && I==1 && println("   eff_sample_size=$(sol.solverstats[:eff_sample_size]), $(maximum(sol.solverstats[:eff_sample_sizes]))")
    debugging && println("\n")

    if dir_output!="" && save_geotiff
        BM.make_geotiff(sol, dir_output*"/h-$flroot-$(BM.getrgi(gid)).tif")
    end
    if dir_output!="" && save_sol
        open(dir_output*"/sol-$flroot-$(BM.getrgi(gid)).jls", "w") do io
            serialize(io, sol)
        end
    end
    if retsol
        return BM.condense_MCMCSol(sol, store_thetas=store_thetas), sol
    else
        return BM.condense_MCMCSol(sol, store_thetas=store_thetas)
    end
end

var2filename(gid, sigma_h_model, sigma_iv_model) =
    dir_output*"$(BM.getrgi(gid))_sh=$(sigma_h_model)_siv=$(sigma_iv_model)"

import OnlineStats
function calculate_prior_update(out1::Vector, drop_sigma)
    out1_ = filter(o-> o isa Dict, out1)
    # if there is no out1 data, return nothing:
    length(out1_)==0 && return nothing, nothing, nothing, nothing

    # find out length without sigma
    if drop_sigma
        d = Dict(out1_[1][:theta0][:names]...)
        # remove all sigmas
        haskey(d, :sigma_iv_model) && delete!(d, :sigma_iv_model)
        haskey(d, :sigma_h_model) && delete!(d, :sigma_h_model)
        # and see what's the last index
        len = maximum(map(maximum, values(d)))
    else
        len = length(out1_[1][:thetas_expect])
        @assert all([length(o[:thetas_expect]) for o in out1_].==len) "Not all o[:thetas_expect] the same length"
    end

    # calculate expectation.
    expect_mean = mean([o[:thetas_expect][1:len] for o in out1_])
    expect_median = median(hcat([o[:thetas_expect][1:len] for o in out1_]...),2)[:] # might be more robust because of outliers
    expect = expect_median

    # use OnlineStats to get the total variance (and quantiles (todo))
    v = OnlineStats.Group([OnlineStats.Variance() for i=1:len])
    q = OnlineStats.Group([OnlineStats.Quantile() for i=1:len])
    for o in out1_
        ns = o[:solverstats][:nsamples]
        v1 = OnlineStats.Group([OnlineStats.Variance(o[:thetas_std][i]^2, o[:thetas_expect][i], OnlineStats.EqualWeight, ns)
                                for i=1:len])
        merge!(v, v1)
    end
    stdev = sqrt.(map(x->OnlineStats.value(x), v))

    return expect, expect_median, expect_mean, stdev
end
