# Invert for the two glacier which have radar

println("Loading BITEModel.jl")
@time using BITEModel
const BM=BITEModel
import KissMCMC
using Base.Test
import JLD2
import Plots

i = 4
niter = [10^3, 3*10^4, 10^5, 3*10^5][i]
nthin = [1, 10, 50, 100][i]

with_good_radar = [409, 428]
region_svalbard = 7
gid = BM.RGIGlacier.(with_good_radar, region_svalbard)[2]

plotyes = false
verbose = false

gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose)
# Inversion
###########
n_1d_theta = 3

run=[:h, :iv, :h_iv][3]
(theta0, logposterior, logposterior1d, logprior, logprior_debug,
 loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
 pmcmc_defaults, fit_target) =
     BM.init_inverse(gb, pp, pm, pn, n_1d_theta=n_1d_theta, run=run,
                     theta0_dict = Dict(:btilde=> zeros(BM.F,n_1d_theta), # NOTE in m/y (error)
                                        :fsl=> zeros(BM.F,n_1d_theta), # error
                                        :temp=> zeros(BM.F,n_1d_theta), # error
                                        :dist_exp=> pm.dist_exp,
                                        :iv_h_exp=> pm.iv_h_exp,
                                        :iv_dist_exp1=> pm.iv_dist_exp1,
                                        #:iv_dist_exp2=> pm.iv_dist_exp2,
                                        )
                     )
theta0_prior = deepcopy(theta0)

# some tests:
print("Time to run forward model:")
@inferred fwdm_fn(theta0.th0);
@time hs1d, ivs1d, hs2d, ivs2d, pm_, taus1d, taus1d_l, taus2d = fwdm_fn(theta0.th0);

print("Time to calculate posterior:")
@inferred logposterior(theta0.th0)
@time logposterior(theta0.th0)[1];


print("Prior value:")
@show logprior(theta0.th0)
print("Likelihood value:")
@show loglikelihood(theta0.th0)[1]
print("Posterior value:")
@show logposterior(theta0.th0)[1]


# import NLopt

# # using NLopt
# #
# # Much faster than Optim.jl. LN_SBPLX seems better than LN_NELDERMEAD

# ic = copy(theta0.th0)

# od_ll = (th, grad) -> (logprior(th)[1]>-Inf ? loglikelihood(th)[1] : -Inf)
# od_post = (th, grad) -> logposterior(th)[1]

# # https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-derivative-free-optimization
# alg = [:LN_COBYLA, :LN_NELDERMEAD, :LN_BOBYQA, :LN_PRAXIS, :LN_SBPLX][end]
# opt = NLopt.Opt(alg, length(ic))

# NLopt.max_objective!(opt, od_post)
# tol = 1e-6
# NLopt.ftol_rel!(opt, tol)
# NLopt.ftol_abs!(opt, tol)
# NLopt.xtol_abs!(opt, tol)
# NLopt.xtol_rel!(opt, tol)
# NLopt.maxeval!(opt, 5000)

# @time maxf, theta_start, ret = NLopt.optimize(opt, ic)
# good_returns = [:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED]
# if !(ret in good_returns)
#     println("Optimization not successful: $ret\n   Best max found=$maxf)")
# else
#     println("Optimization successful.\n     Max found=$maxf")
# end
# theta0.th0[:] = theta_start


#############
# MCMC
############
pmcmc = BM.MCMCNum(niter = niter,
                   nthin = nthin,
                   sigma_ppdf = ones(length(theta0.th0))*3e-2,
                   alg=[:emcee, :metro][1]) # metro needs many more iterations:
                                            # - emcee good with 10^5
                                            # - metro needs 10^?

# sample posterior
thetas, theta_expect, theta_mode, accept_ratio, blobs, lpost, mc_hs1d, mc_ivs1d,
  mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = mcmc(logposterior, theta0, pmcmc; verbose=true)

# sample the prior to get that distribution too:
pmcmc_prior = BM.MCMCNum(niter=5*10^5, nthin=1000, ball_radius=0.45, drop_low_accept_ratio=false)
thetas_prior, theta_expect_prior, theta_mode_prior, accept_ratio_prior, _, _, mc_hs1d_prior, mc_ivs1d_prior =
    mcmc(logprior, theta0_prior, pmcmc_prior; verbose=true)

varnames = BM.get_varnames(theta0)
KissMCMC.print_results(thetas, accept_ratio, names=varnames)
KissMCMC.print_results(mc_hs1d, 0.0, names=["h" for i=1:size(mc_hs1d,1)])

JLD2.@save("out-$(BM.getrgi_nr(gid)).jld2",
                  gid, pmcmc, pm, pn, run, n_1d_theta,
           #theta0, gb, don't work because of functions,
           # blobs, # does't work in JLD2
                  thetas, theta_expect, accept_ratio, mc_hs1d, mc_ivs1d,
                  mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d,
                  thetas_prior, theta_expect_prior, accept_ratio_prior)

JLD.@save("out-$(BM.getrgi_nr(gid)).jld",
                  # gid, pm, pn, run, n_1d_theta, #pmcmc does not work in JLD
           #theta0, gb, don't work because of functions,
           blobs,
                  # thetas, theta_expect, accept_ratio, mc_hs1d, mc_ivs1d,
                  # mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d,
          # thetas_prior, theta_expect_prior, accept_ratio_prior
          )

# JLD2.@load("out-409.jld2",
#            gb, pmcmc, pp, pm, pn, run,
#            theta0,
#            thetas, theta_expect, accept_ratio, blobs, mc_hs1d, mc_ivs1d,
#            mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d,
#            thetas_prior, theta_expect_prior, accept_ratio_prior)

pyplot()
fig = BM.plottheta(thetas, theta0, toplot=:all)
Plots.savefig("thetas-$(BM.getrgi_nr(gid)).png", fig)

if plotyes

    display(BM.plotinv1d(gb, blobs, reuse=false))
    display(BM.plotinv1d_iverr(gb, blobs, reuse=false))

    display(BM.plotinv2d(gl, blobs, reuse=false))
    display(BM.plotinv2d_h(gl, blobs, reuse=false))
    display(BM.plotinv2d_iverr(gl, blobs, reuse=false))

    display(BM.plottheta(thetas, theta0, toplot=:all))
    display(BM.plottheta(thetas_prior, theta0, toplot=:all))
    # display(BM.plottheta(thetas, theta0, toplot=:btilde, reuse=false))
    # display(BM.plottheta(thetas, theta0, toplot=:fsl, reuse=false))
    # display(BM.plottheta(thetas, theta0, toplot=:temp, reuse=false))

    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :btilde, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :fsl, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :temp, reuse=false, width=1))

    display(BM.plot2d_h(gl, mh2d))

end


nothing
