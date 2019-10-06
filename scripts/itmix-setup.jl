# needs a variable debugging defined to load

error("Outdated: Remove Rhat")

println("Time to load BITEModel.jl")
@time import BITEModel;
const BM=BITEModel
const I2=BM.ITMIX2
using VAWTools, KissMCMC
using FileIO, JLD
import DataStructures: OrderedDict

# update_cache = false
# verbose = false
# plotyes = false
# itmixversion = 2

const progress_bar = true
const save_results = true

"""
    setup_one_exp(exp, gl, gb, pp, pm, pn,
                      n_1d_btilde, n_1d_fsl, n_1d_temp)

Essentially calls `BM.init_inverse` but also sets up:
- theta0
- removes all unused radar-lines from gb

Also has the extra feature to make more experiments:
- exp0: no radar lines
- exp 17-32: no IV, otherwise as 1.16

Above two only make sense (in the sense that there is something to fit)
for glaciers with IV data.
"""
function setup_one_exp(exp, gl, gb, pp, pm, pn,
                       n_1d_btilde, n_1d_fsl, n_1d_temp,
                       runtyp, fit_target, rm_iv;
                       th0_with_sigma=false)
    gln = BM.ITMIX2.make_ITMIX2_glacier(gl, exp, rm_iv)
    gbn = BM.Bands(gb, gl=gln)

    th0d = BM.theta0_dict_defaults(pm, n_1d_btilde, n_1d_fsl, n_1d_temp,
                                   with_sigma=th0_with_sigma)

    # update IC for some glaciers with otherwise zero density of prior
    # because of terminus flux or continuous glacier prior.
    if gl.gid==BITEModel.ITMIXGlacier(:Mocho, 2)
        # Mocho has zero thickness at the top with the given btilde
        # (which leads to prior==-Inf)
        th0d[:btilde][1] = 1.5
    elseif gl.gid==BITEModel.ITMIXGlacier(:Austfonna, 2)
        th0d[:btilde][1] = 1.0
    elseif gl.gid==BITEModel.ITMIXGlacier(:Freya, 2)
        th0d[:btilde][1] = 0.5
    elseif gl.gid==BITEModel.ITMIXGlacier(:Kesselwandferner, 2)
        th0d[:btilde][1] = 0.5
    end
    return gln, gbn, th0d, BM.init_inverse(gbn, pp, pm, pn,
                                           runtyp=runtyp,
                                           theta0_dict=th0d,
                                           fit_target=fit_target
                                           )...
end

"""
    run_one_exp(exp, gl, gb, pp, pm, pn, pl;
                     n_1d_btilde = 5,
                     n_1d_fsl = 3,
                     n_1d_temp = 1)

Runs one ITMIX2 experiment.  Useful to use in `map, `pmap`, etc.

Also has the extra feature to make more experiments:
- exp0: no radar lines
- exp 17-32: no IV, otherwise as 1.16

Above two only make sense (in the sense that there is something to fit)
for glaciers with IV data.
"""
function run_one_exp(exp, gl, gb, pp, pm, pn, pl, runtyp, fit_target;
                     rm_iv=false,
                     n_1d_btilde=3,
                     n_1d_fsl=3,
                     n_1d_temp=1,
                     dir="results_Werder"
                     )
    # I get an odd bug in pmap with ERROR: LoadError: On worker 2: UndefVarError: h_iv not defined
    # when using symbols for runtyp and fit_target, thus the clutch of passing them as strings:
    runtyp = Symbol(runtyp)

    println("\n --- Running experiment $exp of $(BM.getname(gl))")

    (gln, gbn, th0d, theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) =
         setup_one_exp(exp, gl, gb, pp, pm, pn,
                       n_1d_btilde, n_1d_fsl, n_1d_temp,
                       runtyp, fit_target, rm_iv)
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

    # sample posterior
    sol = BM.mcmc(logposterior, theta0, pmcmc;
                  verbose=debugging,
                  use_progress_meter=progress_bar,
                  misc=Dict(:logposterior=>logposterior,
                            :logposterior1d=>logposterior1d,
                            :logprior=>logprior,
                            :logprior_debug=>logprior_debug,
                            :loglikelihood=>loglikelihood,
                            :loglikelihood1d=>loglikelihood1d,
                            :fwdm_fn=>fwdm_fn,
                            :fit_target=>fit_target))
    debugging && println("MCMC of experiment $exp of $(BM.getname(gl)) done in $(signif(wall_time/60,4)) minutes")
    debugging && println("   Rhat=$(mean(sol.Rhat)), $(maximum(sol.Rhat))")
    #debugging && print_results(thetas, accept_ratio, names=BM.get_varnames(theta0))
    if save_results
        fln_jld = I2.save_itmix2_run(sol.mh2d, sol.sh2d, gln, pl, sol.miv2d, sol.siv2d, dir=dir)[end]
        FileIO.save(fln_jld,
                    "Rhat", sol.Rhat,
                    "sample_size", sol.sample_size,
                    "nthin", sol.nthin,
                    "thetas", sol.thetas,
                    "reshape_revert", sol.reshape_revert,
                    "runtyp", runtyp,
                    "fit_target", fit_target,
                    "theta0_names", sol.theta0.names,
                    "theta0_th0", sol.theta0.th0,
                    "th0_dict", th0d,
                    "wall_time", wall_time)
    end

    return res
end

println("itmix-setup.jl loaded")
