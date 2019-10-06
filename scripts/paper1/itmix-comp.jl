# This fits Unteraargletscher (in its ITMIX incarnation) and compares the
# model on some radar lines (done with the f3.jl and f4.jl scripts)

# copied and adapted from ../itmix2-runs.jl

##########################################
## Fitting
using Parameters
include("../itmix-setup.jl")

######
# Choose glacier & options
######
plotyes = false
itmixversion = 2
r_tracks = [5,34,42]


glacier = :Unteraar
gid = BM.ITMIXGlacier(glacier, 2)
gl,gb,pp,pm,pn,pl = BM.init_forward(gid)

pm = BM.MPara(pm,
              sigma_h_model=20,
              error_model=[BM.ErrorModel.sqrerr, BM.ErrorModel.abserr][1],
              nthin_h=1)

####################
# Do it
####################
runtyp = [:test,:testmid,:prodlow,:prodmid,:prod][3]
fit_target = [BM.FitTarget.h, BM.FitTarget.iv, BM.FitTarget.h_iv, BM.FitTarget.length][1]
n_1d_btilde = 3
n_1d_fsl = 3
n_1d_temp = 1
rm_iv = false

(gln, gbn, th0d,
 theta0, logposterior, logposterior1d,
 logprior, logprior_debug, loglikelihood, loglikelihood1d,
 fwdm_fn, fwdm1d_fn, pmcmc_defaults, fit_target) =
     setup_one_exp(r_tracks, gl, gb, pp, pm, pn,
                   n_1d_btilde, n_1d_fsl, n_1d_temp,
                   runtyp, fit_target, rm_iv)

# some tests:
print("Time to run forward model:")
@time fwdsol = fwdm_fn(theta0.th0);
# display(BM.plot_run(fwdsol, reuse=false, plot2d=[:h,:iv,:ivlog][1]))
print("Prior value:")
@show logprior(theta0.th0)[1]
#logprior_debug(theta0.th0)[1]
print("Posterior value:")
@show logposterior(theta0.th0)[1]
print("Time to calculate posterior:")
@time logposterior(theta0.th0)[1]

########
# max posterior estimation
#######
if false
    using Optim
    println("Optimizing for maximum posterior...")
    theta_max2 = optimize(theta_vec -> -logposterior(theta_vec)[1], float(theta0.th0), iterations=200, f_tol=1.0)
    # these don't work so well:
    #   theta_max2 = optimize(theta_vec -> -logposterior(theta_vec)[1], float(theta0.th0), GradientDescent())
    #theta_max2 = optimize(theta_vec -> -logposterior(theta_vec)[1], float(theta0.th0),  LBFGS())

    @show theta_max2
    error("asdf")
    theta0.th0[:] = theta_max2.minimizer
end

###############
## MCMC fitting
###############
println("Starting MCMC...")
pmcmc = BM.MCMCNum(;pmcmc_defaults...)

# sample posterior
sol = BM.mcmc(logposterior, theta0, pmcmc; verbose=true,
              misc=Dict(:logposterior=>logposterior,
                        :logposterior1d=>logposterior1d,
                        :logprior=>logprior,
                        :logprior_debug=>logprior_debug,
                        :loglikelihood=>loglikelihood,
                        :loglikelihood1d=>loglikelihood1d,
                        :fwdm_fn=>fwdm_fn,
                        :fit_target=>fit_target))
varnames = BM.get_varnames(theta0)
using KissMCMC
print_results(sol.thetas, sol.solverstats[:accept_ratio], sol.solverstats[:effective_sample_sizes], names=varnames)

# sample the prior to get that distribution too:
println("\n\nStarting MCMC of the prior to get that distribution too...")

pmcmc_prior = BM.MCMCNum(pmcmc, niter=10^5,
                         nthin = 50,
                         nburnin = 5*10^4,
                         nchains=100)

sol_prior = BM.mcmc(logprior, theta0, pmcmc_prior; verbose=true)
print_results(sol_prior.thetas, sol_prior.solverstats[:accept_ratio], sol_prior.solverstats[:effective_sample_sizes], names=varnames)

date = "$(lowercase(Dates.monthabbr(Dates.today())))$(VAWTools.int2str2(Dates.day( Dates.today())))"
fl = "figs/f2-$date-$(join(r_tracks, "-"))-$fit_target-$(pm.error_model)-sigma$(round(Int,pm.sigma_h_model))-thin$(pm.nthin_h).jls"
println("\n $fl")
open(fl, "w") do io
    serialize(io, (sol, sol_prior))
end
