using KissMCMC, BITEModel, VAWTools; const BM=BITEModel
const I2 = BM.ITMIX2
#include("itmix2-runs.jl")

pmcmc = BM.MCMCNum(;pmcmc_defaults...,
                   niter=400,
                   nthin = 1,
                   nburnin=10,
                   nchains=4)

# number of points to fit 1D quantities on
n_1d_theta = 3

good_exps = Dict(gl=>[] for gl in BM.ITMIXv2_glaciers)
bad_gl = Dict(gl=>[] for gl in BM.ITMIXv2_glaciers)
bb = map(x->BM.ITMIXGlacier(x,2), [:Starbuck, :Columbia])
bb = map(x->BM.ITMIXGlacier(x,2), [:Aqqutikitsoq])

for priority=1:4
    println("\n\n++++++++++++++++ PRIORITY $priority ++++++++++++++++++++++++++\n\n")
    for gid in I2.priorities[priority]
        exp = 1 # to make it available
        println("\n\n============================\nRunning Glacier $(BM.getname(gid))")
        #    try
        gl,gb,pp,pm,pn,pl = BM.init_forward(gid; pl_kws...);
        for exp=1:I2.nexps
            println("\n\n---------------\nRunning experiment $exp of $(BM.getname(gid))")
            gln = BM.ITMIX2.make_ITMIX2_glacier(gl, exp)
            gbn = BM.Bands(gb, gl=gln)
            th0d = BM.theta0_dict_defaults(n_1d_theta, pm)
            (theta0, logposterior, logposterior1d, logprior, logprior_debug,
             loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
             pmcmc_defaults, fit_target) =
                 BM.init_inverse(gbn, pp, pm, pn, n_1d_theta = n_1d_theta,
                                                               theta0_dict=th0d)
            # some tests:
            fwdm_fn(theta0.th0);
            print("    Time to run forward model:")
            @time hs1d, ivs1d, hs2d, ivs2d, pm_, taus1d, taus1d_l, taus2d =
                fwdm_fn(theta0.th0);
            print("    Posterior value:")
            @show logposterior(theta0.th0)[1]
            print("    Time to calculate posterior:")
            @time logposterior(theta0.th0)[1]

            # sample posterior
            thetas, theta_expect, accept_ratio, blobs,
            Rhat, sample_size, nthin, reshape_revert,
            mc_hs1d, mc_ivs1d,
            mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = mcmc(logposterior, theta0, pmcmc; verbose=true);
            print_results(thetas, accept_ratio, names=BM.get_varnames(theta0))
            push!(good_exps[gid], exp)
            I2.save_itmix2_run(mh2d, gln, pl)
        end
        # catch e
            #     println("\nFailed with $e\n\n")
        #     push!(bad_gl[gid], (exp, e))
        # end
    end
end
