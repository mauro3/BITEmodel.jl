using Optim
using BITEModel
const BM=BITEModel

######
# Choose glacier & options
######
@def plotyes = false

# ITMIX glacier with IV data:
glaciers_iv = [:Brewster, :Hellstugubreen, :NorthGlacier,
               :SouthGlacier, :Tasman, :Unteraar, :Synthetic1,
               :Synthetic2, :Synthetic3]
glacier = glaciers_iv[1]
glacier = :Unteraar
gid = [ITMIXGlacier(glacier), SyntheticGlacier(:bench)][1]
pl_kws = Dict{Symbol,Any}()
if isa(gid, ITMIXGlacier)
    pl_kws[:use_glogem] = true
end

gl,gb,pp,pm,pn,pl = BM.init_forward(gid; pl_kws...)

####################
# Do it
####################

(theta0, logposterior, logposterior1d, logprior, logprior_debug,
 loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
 pmcmc_defaults, fit_target) =
     BM.init_inverse(gb, pp, pm, pn, n_1d_theta=n_1d_theta, run=run,
                     theta0_dict = theta0_dict     )
# some tests:
fwdm_fn(theta0.th0);
print("Time to run forward model:")
@time fwdsol = fwdm_fn(theta0.th0);
print("Posterior value:")
@show logposterior(theta0.th0)[1]
print("Time to calculate posterior:")
@time logposterior(theta0.th0)[1]

########
# max posterior estimation
#######
if false
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
pmcmc = BM.MCMCNum(pmcmc,
                   niter=0.1*10^6,
                   nthin = 50)

# sample posterior
res = mcmc(logposterior, theta0, pmcmc; verbose=true)
varnames = BM.get_varnames(theta0)
using KissMCMC
print_results(res.thetas, res.accept_ratio, names=varnames)

# sample the prior to get that distribution too:
pmcmc_prior = BM.MCMCNum(pmcmc, niter=0.1*10^6,
                         nthin = 50)
res_prior = mcmc(logprior, theta0, pmcmc_prior; verbose=true)



# BM.savemcmc(thetas, blobs, gb,
#                   pp, pm, pn, pmcmc,
#                   theta0, run, sigmas, error_on_dirty=false)

#end
# ### plots

if plotyes
    display(BM.plotinv1d(gb, blobs, reuse=false))
    display(BM.plotinv1d_err(gb, blobs, reuse=false))

    display(BM.plotinv2d(gl, blobs, reuse=false))
    display(BM.plotinv2d_h(gl, blobs, reuse=false))
    display(BM.plotinv2d_iverr(gl, blobs, reuse=false))

    display(BM.plottheta(thetas, theta0))
    display(BM.plottheta(thetas, theta0, toplot=:btilde, reuse=false))
    display(BM.plottheta(thetas, theta0, toplot=:fsl, reuse=false))
    display(BM.plottheta(thetas, theta0, toplot=:temp, reuse=false))

    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :btilde, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :fsl, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :temp, reuse=false, width=1))
end
