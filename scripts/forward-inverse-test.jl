# This runs the forward model to make an ice thickness and then inverts for it.

# TODO:
# - run several tests
# - make into a test suite
# - figure out whether dropping temp or fsl is better

error("not working")

println("Loading BITEModel.jl")
using Optim
@time using BITEModel
const BM=BITEModel
import NLopt

plotyes = false
verbose = false
gid = SyntheticGlacier(:bench)

# Run it
########
println("Time to to make load-paras")

# parameters to play with
@show meas_dens = [:low, :normal, :high][2]
@show fsl0 = [0.1, 0.5, 0.9][1]
@show temp0 = [-10, -0.5][2]
@show sigma_small = 0.5

# TODO: This is bad as gl,pp,pm,pn,pl get loaded again below!
# I am not sure how the parameters flow...
# However, I also need gl in the
pl = BM.LoadPara(gid)
dt, pl = BM.make_datatable(gid,pl)
# Change defaults for :temp and/or :fsl
dt.fsl.opts[:fsl] = fsl0
dt.temp.opts[:temp] = temp0
gl,pp,pm,pn,pl = load_glacier(dt,pl)

# add the parameter to the options so the
# forward model is used to produce the thickness:
error_fac = 0

h_opts = pl.dataset_opts[BM.ThicknessData]
h_opts[:use_forward_model] = true
if meas_dens==:low
    h_opts[:profile_xlocs] = [2000.0]
    h_opts[:profile_ylocs] = [0.0]
elseif meas_dens==:high
    h_opts[:profile_xlocs] = [100:100:5900;]
    h_opts[:profile_ylocs] = [-250:250:250;]
end
h_opts[:fwd_para] = (gl, pp, pm, pn)
h_opts[:err_sigma] = 10.0 * error_fac
h_opts[:err_cor_length] = 600.0
h_opts[:profile_margin_dist] = 10.0

iv_opts = pl.dataset_opts[BM.IVData]
iv_opts[:use_forward_model] = true
iv_opts[:fwd_para] = (gl, pp, pm, pn)
iv_opts[:err_sigma] = 1.0 * error_fac
iv_opts[:err_cor_length] = 500.0
if meas_dens==:high
    iv_opts[:iv_grid] = (gl.dem.x, gl.dem.y)
elseif meas_dens==:low
    iv_opts[:iv_grid] = (gl.dem.x[1:10:end], gl.dem.y[1:10:end])
end
# for convenience:
h_opts_dt = dt.thickness.opts
iv_opts_dt = dt.iv.opts

# Now update all
dt, pl = BM.make_datatable(gid,pl)
# Again change defaults for :temp and/or :fsl
dt.fsl.opts[:fsl] = fsl0
dt.temp.opts[:temp] = temp0
gl,pp,pm,pn,pl = load_glacier(dt,pl)
gb,pm = BM.make_bands(gl, pp, pm, pn);

# Now, finally update pm (only benign things!!!)
## set model uncertainty to almost 0
pm = BM.MPara(pm, sigma_h_model=sigma_small, sigma_iv_model=sigma_small)

print("Time to run forward model:"); gc()
@time fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

qtot1d, qd1d = BM.postproc1d(gb, pp, pm)

if plotyes
    print("Time to plot:")
    @time display(BM.plot_run(fwdsol,
                              reuse=false, plot2d=[:h,:iv,:ivlog][1]))
    # @time display(BM.plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d,
    #                           ivs2d, taus1d_l, pp, pm,
    #                           reuse=false, plot2d=[:h,:iv,:ivlog][2]))
end

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

th0_true = deepcopy(theta0.th0)

if h_opts[:err_sigma]==0 && iv_opts[:err_sigma]==0
    # check that th0_true is indeed minimal
    good = true
    maxll = loglikelihood(th0_true)[1]
    for delta = [1e-8, 1e-4, 1e-3,1e-2,1,10,100,1000]
        for i=1:length(th0_true)
            th = copy(th0_true)
            th[i] += delta
            if logprior(th)>-Inf && loglikelihood(th)[1]>maxll
                good = false
                warn("loglikelihood < for $th, delta=$delta, i=$i")
                println(loglikelihood(th)[1]-maxll)
            end
            th[i] -= 2delta
            if logprior(th)>-Inf && loglikelihood(th)[1]>maxll
                good = false
                warn("loglikelihood < for $th, delta=-$delta, i=$i")
                println(loglikelihood(th)[1]-maxll)
            end
        end
    end
    if good
        println("Log-likelihood indeed minimal at th0_true = $maxll")
    else
        println("Log-likelihood NOT minimal at th0_true = $maxll")
    end
else
    println("Log-likelihood at th0_true = $(loglikelihood(th0_true)[1])")
end

# some tests:
fwdm_fn(theta0.th0);
print("Time to run forward model:")
@time hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l = fwdm_fn(theta0.th0);
print("Time to calculate posterior:")
@time logposterior(theta0.th0)[1]


print("Prior value:")
@show logprior(th0_true)
print("Likelihood value:")
@show loglikelihood(th0_true)[1]
print("Posterior value:")
@show logposterior(th0_true)[1]

########
# max posterior estimation
#######
println("Optimizing for maximum posterior...")

# # Using Optim.jl
# #
# # Work only reliably when close to the maximum.  Also the tolerances need to be reasonably low.
# #
# # - all converge if started at the minimum.  However Nealder-Mead takes 999 objective calls!
# # - Nealder-Mead claims to not be converged even though it gets closest to the max
# #
# opts = Optim.Options(x_tol=1e-8,
#                      f_tol=1e-8,
#                      g_tol=1e-7,
#                      iterations=1000)

# srand(1)
# ic = th0_true + rand(length(th0_true)) * 0.
# #ic = theta0.th0
# od = theta_vec -> logprior(theta_vec)[1]>-Inf ? -loglikelihood(theta_vec)[1] : Inf
# @show @time theta_max_likelihood1 = optimize(od, ic, opts)
# @show @time theta_max_likelihood2 = optimize(od, ic, GradientDescent(), opts)
# @show @time theta_max_likelihood2 = optimize(od, ic, AcceleratedGradientDescent(), opts)
# @show @time theta_max_likelihood3 = optimize(od, ic,  BFGS(), opts)

# println("Optimizing for maximum posterior...")
# @show @time theta_max_posterior = optimize(theta_vec -> -logposterior(theta_vec)[1], ic, BFGS(), opts)
# # these don't work so well:
# #   theta_max2 = optimize(theta_vec -> -logposterior(theta_vec)[1], float(theta0.th0), GradientDescent())
# #theta_max2 = optimize(theta_vec -> -logposterior(theta_vec)[1], float(theta0.th0),  LBFGS())


# using NLopt
#
# Much faster than Optim.jl. LN_SBPLX seems better than LN_NELDERMEAD

ic = th0_true + rand(length(th0_true)) * 0.01

od_ll = (th, grad) -> (logprior(th)[1]>-Inf ? loglikelihood(th)[1] : -Inf)
od_post = (th, grad) -> logposterior(th)[1]

# https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-derivative-free-optimization
alg = [:LN_COBYLA, :LN_NELDERMEAD, :LN_BOBYQA, :LN_PRAXIS, :LN_SBPLX][end]
opt = NLopt.Opt(alg, length(ic))

NLopt.max_objective!(opt, od_post)
tol = 1e-6
NLopt.ftol_rel!(opt, tol)
NLopt.ftol_abs!(opt, tol)
NLopt.xtol_abs!(opt, tol)
NLopt.xtol_rel!(opt, tol)
NLopt.maxeval!(opt, 5000)

@time maxf, theta_start, ret = NLopt.optimize(opt, ic)
good_returns = [:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED]
if !(ret in good_returns)
    println("Optimization not successful: $ret\n   Best max found=$maxf vs true $(logposterior(th0_true)[1])")
else
    println("Optimization successful.\n     Max found=$maxf vs true $(logposterior(th0_true)[1])")
end
theta0.th0[:] = theta_start

#############
# MCMC
############
pmcmc = BM.MCMCNum(niter = 10^3,
                   nthin = 20,
                   sigma_ppdf = ones(length(theta0.th0))*3e-2,
                   alg=[:emcee, :metro][1]) # metro needs many more iterations:
                                            # - emcee good with 10^5
                                            # - metro needs 10^?

# sample posterior
thetas, theta_expect, theta_mode, accept_ratio, blobs, lpost, mc_hs1d, mc_ivs1d,
  mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = mcmc(logposterior,
                                                                        theta0, pmcmc; verbose=true)

# sample the prior to get that distribution too:
# pmcmc_prior = BM.MCMCNum(niter=10^4,
#                          nthin = 10)
# thetas_prior, theta_expect_prior, accept_ratio_prior = mcmc(logprior, theta0, pmcmc_prior; verbose=true)


# pmcmc_prior2 = BM.MCMCNum(niter=10^6,
#                           nthin = 1, alg=:metro,
#                           sigma_ppdf=ones(length(theta0.th0))/5  )
# thetas_prior2, theta_expect_prior, accept_ratio_prior = mcmc(logprior, theta0, pmcmc_prior2; verbose=true)


varnames = BM.get_varnames(theta0)
using KissMCMC
print_results(thetas, accept_ratio, names=varnames, theta_true=th0_true)

if plotyes
    display(BM.plotinv1d(gb, blobs, reuse=false))
    display(BM.plotinv1d_iverr(gb, blobs, reuse=false))

    display(BM.plotinv2d(gl, blobs, reuse=false))
    display(BM.plotinv2d_h(gl, blobs, reuse=false))
    display(BM.plotinv2d_iverr(gl, blobs, reuse=false))

    display(BM.plottheta(thetas, theta0, toplot=:all))
    # display(BM.plottheta(thetas, theta0, toplot=:btilde, reuse=false))
    # display(BM.plottheta(thetas, theta0, toplot=:fsl, reuse=false))
    # display(BM.plottheta(thetas, theta0, toplot=:temp, reuse=false))

    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :btilde, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :fsl, reuse=false, width=1))
    display(BM.plottheta_violin((thetas, thetas_prior), theta0, :temp, reuse=false, width=1))
end


nothing
