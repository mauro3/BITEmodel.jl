
# ## Run forward model
# ####################
# priority_cutoff = 1
# for gl in instances(I2Gls.ITMIX2Glaciers)
#     if getpriority(tab, gl)>priority_cutoff
#         continue
#     end
#     gid = ITMIX2Glacier(gl)
#     for (_, exp, radar) in experiments(tab, gl)
#         # note that the radar-lines do not matter here as we don't fit the model
#         @time gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
#                                                   update_cache=update_cache,
#                                                   use_glogem=use_glogem,
#                                                   radarlines=getradar(radar))
#         @time out = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
#         save_itmix2(gl, exp, out)
#     end
# end

##########################################
## Fitting
using Parameters
include("itmix-setup.jl")

######
# Choose glacier & options
######
plotyes = false
itmixversion = 2
itmix_exp = 10

glacier = :Unteraar
gid = BM.ITMIXGlacier(glacier, 2)
gl,gb,pp,pm,pn,pl = BM.init_forward(gid)

####################
# Do it
####################
runtyp = [:test,:mid,:prodlow,:prodmid,:prod][3]
fit_target = [BM.FitTarget.h, BM.FitTarget.iv, BM.FitTarget.h_iv, BM.FitTarget.length][3]
n_1d_btilde = 3
n_1d_fsl = 3
n_1d_temp = 1
rm_iv = false

(gln, gbn, th0d,
 theta0, logposterior, logposterior1d,
 logprior, logprior_debug, loglikelihood, loglikelihood1d,
 fwdm_fn, fwdm1d_fn, pmcmc_defaults, fit_target) =
     setup_one_exp(itmix_exp, gl, gb, pp, pm, pn,
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
sol = BM.mcmc(logposterior, theta0, pmcmc; verbose=true);
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

# BM.savemcmc(thetas, blobs, gb,
#                   pp, pm, pn, pmcmc,
#                   theta0, run, sigmas, error_on_dirty=false)

#end
# ### plots

if plotyes
    display(BM.plotinv1d(sol, reuse=false))
    display(BM.plotinv1d_err(sol, reuse=false))

    display(BM.plotinv2d(sol, reuse=false))
    display(BM.plotinv2d_h(sol, reuse=false))
    display(BM.plotinv2d_iverr(sol, reuse=false))
    @unpack thetas, theta0 = sol
    display(BM.plottheta(thetas, theta0))
    display(BM.plottheta(thetas, theta0, toplot=:btilde, reuse=false))
    display(BM.plottheta(thetas, theta0, toplot=:fsl, reuse=false))
    display(BM.plottheta(thetas, theta0, toplot=:temp, reuse=false))

    display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :btilde, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :fsl, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, sol_prior.thetas), theta0, :temp, reuse=false, width=1))
end
