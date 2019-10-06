# This runs the forward model to make an ice thickness and then inverts for it.

# TODO:
# - run several tests
# - make into a test suite
# - figure out whether dropping temp or fsl is better


println("Loading BITEModel.jl")
@everywhere begin
    using Plots, StatPlots
    if strip(readstring(`hostname`))=="x230"
        pyplot()
    else
        ENV["GKSwstype"]="100" # for headless/no display
        gr()
    end
    using Optim
    using BITEModel
    const BM=BITEModel
    import JLD2,JLD
    using FileIO
#    import NLopt
    import KissMCMC
    using Base.Iterators

    plotyes = false
    verbose = false
    gid = [BM.SyntheticGlacier(:bench), BM.ITMIXGlacier(:Unteraar)][1]

    # Run it
    ########
    mcmc_niter = 2*10^5
    # parameters to play with
    @show error_fac = [0,1][1]
    seed = 5

    @show meas_dens_s = [:low, :normal, :high]
    @show fsl0s = [0.1, 0.5, 0.9]
    @show temp0s = [-10, -0.5]
    @show sigma_smalls = [5.0, 0.5]
    n_1d_theta = 3
    pm0 = BM.MPara()
    theta0_dict_def = Dict(:btilde=> zeros(BM.F,n_1d_theta), # NOTE in m/y (error)
                           :fsl=> zeros(BM.F,n_1d_theta), # error
                           :temp=> zeros(BM.F,n_1d_theta), # error
                           :dist_exp=> pm0.dist_exp,
                           :iv_h_exp=> pm0.iv_h_exp,
                           :iv_dist_exp1=> pm0.iv_dist_exp1)
    #:iv_dist_exp2=> pm.iv_dist_exp2, this does not do much!?
    vars2pop = [:na] #[:temp, :fsl, :na]

    runs=[:h, :iv, :h_iv]


    function var2filename(run, meas_dens, fsl0, temp0, sigma_small, popit)
        if error_fac==0
            out = "output/fwd-inv/$(BM.getname(gid))_run-$(run)_md-$(meas_dens)_fsl-$(fsl0)_temp-$(temp0)_sigma-$(sigma_small)_sans-$(popit)"
        else
            out = "output/fwd-inv/$(BM.getname(gid))_with-err_run-$(run)_md-$(meas_dens)_fsl-$(fsl0)_temp-$(temp0)_sigma-$(sigma_small)_sans-$(popit)"
        end
    end
end

# Use Core.println
# https://discourse.julialang.org/t/status-of-threads/3729
# Threads.@threads does not work due to plotting and probably printing

# to use up to 27 workers use product:
@sync @parallel for (run,meas_dens,fsl0) in collect(product(runs,meas_dens_s,fsl0s))
    for temp0 in temp0s, sigma_small in sigma_smalls, popit in vars2pop
        tic()
        fl = var2filename(run, meas_dens, fsl0, temp0, sigma_small, popit)
        Core.println(fl)
        theta0_dict = deepcopy(theta0_dict_def)
        pop!(theta0_dict, popit, nothing)

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

        # reseed to make the same error for each run:
        srand(seed)

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
        h_opts[:err_cor_length] = 300.0
        h_opts[:profile_margin_dist] = 30.0

        iv_opts = pl.dataset_opts[BM.IVData]
        iv_opts[:use_forward_model] = true
        iv_opts[:fwd_para] = (gl, pp, pm, pn)
        iv_opts[:err_sigma] = 1.0 * error_fac
        iv_opts[:err_cor_length] = 300.0
        if meas_dens==:high
            iv_opts[:iv_grid] = (gl.dem.x, gl.dem.y)
        elseif meas_dens==:low
            iv_opts[:iv_grid] = (gl.dem.x[1:10:end], gl.dem.y[1:10:end])
        end

        # if ITMIX swap out IV and thickness
        if gid isa BM.ITMIXGlacier
            pl.dataset_LOADERS[BM.ThicknessData] = BM.SyntheticBenchLoader
            pl.dataset_LOADERS[BM.IVData] = BM.SyntheticBenchLoader
        end

        # Now update all
        dt, pl = BM.make_datatable(gid,pl)
        # for convenience:
        h_opts_dt = dt.thickness.opts
        iv_opts_dt = dt.iv.opts

        # Again change defaults for :temp and/or :fsl
        dt.fsl.opts[:fsl] = fsl0
        dt.temp.opts[:temp] = temp0
        gl,pp,pm,pn,pl = load_glacier(dt,pl)
        gb,pm = BM.make_bands(gl, pp, pm, pn);

        # Now, finally update pm (only benign things!!!)
        ## set model uncertainty to almost 0
        pm = BM.MPara(pm, sigma_h_model=sigma_small, sigma_iv_model=sigma_small)

        #print("Time to run forward model:"); gc()
        hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
            BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

        qtot1d, qd1d = BM.postproc1d(gb, pp, pm)

        fig = BM.plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d,
                          ivs2d, taus1d_l, pp, pm,
                          reuse=false, plot2d=[:h,:iv,:ivlog][1])
        Plots.savefig(fig, fl*"--map.png")

        # Inversion
        ###########
        (theta0, logposterior, logposterior1d, logprior, logprior_debug,
         loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
         pmcmc_defaults, fit_target) =
            BM.init_inverse(gb, pp, pm, pn, n_1d_theta=n_1d_theta, run=run,
                            theta0_dict = theta0_dict     )

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
                        Core.println(loglikelihood(th)[1]-maxll)
                    end
                    th[i] -= 2delta
                    if logprior(th)>-Inf && loglikelihood(th)[1]>maxll
                        good = false
                        warn("loglikelihood < for $th, delta=-$delta, i=$i")
                        Core.println(loglikelihood(th)[1]-maxll)
                    end
                end
            end
        end

        # some tests:
        fwdm_fn(theta0.th0);
        # print("Time to run forward model:")
        hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l = fwdm_fn(theta0.th0);
        # print("Time to calculate posterior:")
        logposterior(theta0.th0)[1]

        #         ###########
        #         # using NLopt
        #         #
        #         # Much faster than Optim.jl. LN_SBPLX seems better than LN_NELDERMEAD

        #         ic = th0_true + rand(length(th0_true)) * 0.01

        #         od_ll = (th, grad) -> (logprior(th)[1]>-Inf) ? loglikelihood(th)[1] : -Inf
        #         od_post = (th, grad) -> logposterior(th)[1]

        #         # https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/#local-derivative-free-optimization
        #         alg = [:LN_COBYLA, :LN_NELDERMEAD, :LN_BOBYQA, :LN_PRAXIS, :LN_SBPLX][end]
        #         opt = NLopt.Opt(alg, length(ic))

        #         NLopt.max_objective!(opt, od_post)
        #         tol = 1e-6
        #         NLopt.ftol_rel!(opt, tol)
        #         NLopt.ftol_abs!(opt, tol)
        # NLopt.xtol_abs!(opt, tol)
        # NLopt.xtol_rel!(opt, tol)
        # NLopt.maxeval!(opt, 5000)

        # maxf, theta_start, ret = NLopt.optimize(opt, ic)
        # good_returns = [:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED]
        # if !(ret in good_returns)
        #     Core.println("Optimization not successful: $ret\n   Best max found=$maxf vs true $(logposterior(th0_true)[1])")
        #     # else
        #     #     println("Optimization successful.\n     Max found=$maxf vs true $(logposterior(th0_true)[1])")
        # end
        # theta0.th0[:] = theta_start

        theta0.th0[:] = th0_true + rand(length(th0_true)) * 0.01

        #############
        # MCMC
        ############
        pmcmc = BM.MCMCNum(niter = mcmc_niter,
                           nthin = 20,
                           sigma_ppdf = ones(length(theta0.th0))*3e-2,
                           alg=[:emcee, :metro][1]) # metro needs many more iterations:
        # - emcee good with 10^5
        # - metro needs 10^?

        # sample posterior
        thetas, theta_expect, accept_ratio, blobs, mc_hs1d, mc_ivs1d,
        mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = mcmc(logposterior,
                                                                              theta0, pmcmc;
                                                                              verbose=false,
                                                                              use_progress_meter=false)

        # Figure out errors of h and iv versus the truth
        # 1D
        true_h1d = h_opts_dt[:fwd_true_h1d]
        true_iv1d = iv_opts_dt[:fwd_true_iv1d]
        true_h_th, true_iv_th = BM.h_iv_at_thetas(theta0, true_h1d, true_iv1d)
        h_th, iv_th = BM.h_iv_at_thetas(theta0, mc_hs1d, mc_ivs1d)

        # do some scatter plots of these:
        err_h_th = h_th .- true_h_th
        err_iv_th = iv_th .- true_iv_th

        rel_err_h_th = (h_th .- true_h_th)./true_h_th
        rel_err_iv_th = (iv_th .- true_iv_th)./true_h_th

        # # 2D
        # hfn = h_opts_dt[:fwd_thick_fn]
        # ivfn = iv_opts_dt[:fwd_iv_fn]
        # true_h2d = Gridded(gl.dem.x, gl.dem.y, hfn.(gl.dem.x, gl.dem.y'));
        # true_iv2d = Gridded(gl.dem.x, gl.dem.y, ivfn.(gl.dem.x, gl.dem.y'));

        # # sample the prior to get that distribution too:
        # pmcmc_prior = BM.MCMCNum(niter=5*10^4,
        #                          nthin = 10)
        # thetas_prior, theta_expect_prior, accept_ratio_prior = mcmc(logprior, theta0, pmcmc_prior; verbose=true)


        # pmcmc_prior2 = BM.MCMCNum(niter=10^6,
        #                           nthin = 1, alg=:metro,
        #                           sigma_ppdf=ones(length(theta0.th0))/5  )
        # thetas_prior2, theta_expect_prior, accept_ratio_prior = mcmc(logprior, theta0, pmcmc_prior2; verbose=true)


        varnames = BM.get_varnames(theta0)
        #print_results(thetas, accept_ratio, names=varnames)

        fig = BM.plottheta(thetas, theta0, toplot=:all, size=(1200,1200))
        Plots.savefig(fig, fl*"--thetas.png")

        # plot of h-errs against other variables
        tmp = copy(thetas)
        tmp[theta0.names[:fsl],:] = err_h_th
        va = strip.(varnames)
        va[theta0.names[:fsl]] = ["h_top","h_mid","h_low"]
        fig = StatPlots.cornerplot(tmp'; size=(1200,1200), label=va)
        Plots.savefig(fig, fl*"--thetas-h-temp.png")

        tmp = copy(thetas)
        tmp[theta0.names[:temp],:] = err_h_th
        va = strip.(varnames)
        va[theta0.names[:temp]] = ["h_top","h_mid","h_low"]
        fig = StatPlots.cornerplot(tmp'; size=(1200,1200), label=va)
        Plots.savefig(fig, fl*"--thetas-h-fsl.png")

        tmp = copy(thetas)
        tmp[theta0.names[:fsl],:] = err_iv_th
        va = strip.(varnames)
        va[theta0.names[:fsl]] = ["iv_top","iv_mid","iv_low"]
        fig = StatPlots.cornerplot(tmp'; size=(1200,1200), label=va)
        Plots.savefig(fig, fl*"--thetas-iv-temp.png")

        tmp = copy(thetas)
        tmp[theta0.names[:temp],:] = err_iv_th
        va = strip.(varnames)
        va[theta0.names[:temp]] = ["iv_top","iv_mid","iv_low"]
        fig = StatPlots.cornerplot(tmp'; size=(1200,1200), label=va)
        Plots.savefig(fig, fl*"--thetas-iv-fsl.png")


        # summaries
        # as JLD
        res = KissMCMC.summarize_run(thetas, theta_true=th0_true, names=varnames)
        save(fl*".jld2", "res", res)

        # as text
        dir_,fll = splitdir(fl)
        ar = hcat([[strip(string(v)),e] for (v,e) in zip(res[:var],res[:err])]...)
        l1,l2 = split(reprmime("text/plain", ar),'\n')[2:3]
        write(fl*"--summary", "$fll\n    $l1\n    $l2\n")
        toc()
    end
end
