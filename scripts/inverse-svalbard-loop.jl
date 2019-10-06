# invert all the Svalbard glaciers


@everywhere debugging = false
debugging && println("\n\nDEBUG RUN !!!!!!!!\n\n")
include("svalbard-setup.jl")

@everywhere begin
    using BSON
    include("svalbard-setup.jl")
    # runs fast if true

    update_cache = true
    svalbard_nrs = debugging ? svalbard_nrs[1:20] : svalbard_nrs

    import Plots, StatPlots
    import Base.Test: @inferred
    import DataStructures: OrderedDict
    Plots.pyplot()
    #Plots.gr()
    # if strip(readstring(`hostname`))=="x230"
    #     pyplot()
    # else
    #     ENV["GKSwstype"]="100" # for headless/no display
    #     gr()
    # end
    import JLD2,JLD
    import FileIO
    import KissMCMC

    # Run it
    ########
    n_1d_theta = 3

    mcmc_niter = debugging ? 10^3 : 5*10^5
    # parameters to play with
    srand(5)

    sigma_h_model = [0.0, 10.0]
    sigma_iv_model = [0.0, 10.0]
    pm0 = BM.MPara()

    n_1d_btilde=3
    n_1d_fsl=3
    n_1d_temp=1

    fit_target = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.length][2]
    runtyp = ["test", "testmid", "prodlow", "prod"][4]

    # used for plotting:
    theta0_with_h_iv_dict = OrderedDict(:btilde=> zeros(BM.F,n_1d_theta), # NOTE in m/y (error)
                                        :fsl=> zeros(BM.F,n_1d_theta), # error
                                        :dist_exp=> pm0.dist_exp,
                                        :h => zeros(BM.F,n_1d_theta), # for plotting
                                        :iv => zeros(BM.F,n_1d_theta), # for plotting
                                        )
    theta0_with_h_dict = deepcopy(theta0_with_h_iv_dict)
    pop!(theta0_with_h_dict, :iv)
    theta0_with_iv_dict = deepcopy(theta0_with_h_iv_dict)
    pop!(theta0_with_iv_dict, :h)

    dir_output = "output/svalbard/04dec/"*runtyp
end

!isdir(dir_output) && mkpath(dir_output)

@everywhere begin
    function var2filename(gid, sigma_h_model, sigma_iv_model)
        out = dir_output*"$(BM.getrgi(gid))_sh=$(sigma_h_model)_siv=$(sigma_iv_model)"
    end
end

#@sync @parallel for nr in svalbard_nrs[1:end÷2]
#@sync @parallel for nr in svalbard_nrs[end÷2+1:end]
@sync @parallel for nr in svalbard_nrs
    for sh in sigma_h_model, siv in sigma_iv_model
        gid = BM.RGIGlacier(nr, region_svalbard)
        fl = var2filename(gid, sh, siv)
        if debugging
            println(fl)
            tic()
        end
        try
            dt, pl = init_svalbard_gl_dt(gid)
            dt.para.opts[:pm][:sigma_h_model] = sh
            dt.para.opts[:pm][:sigma_iv_model] = siv
            gl,gb,pp,pm,pn,pl,dt = init_svalbard_gl(dt, pl)

            # Inversion
            ###########
            th0d = BM.theta0_dict_defaults(pm, n_1d_btilde, n_1d_fsl, n_1d_temp)

            (theta0, logposterior, logposterior1d, logprior, logprior_debug,
             loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
             pmcmc_defaults, fit_target) =
                BM.init_inverse(gb, pp, pm, pn,
                                runtyp=Symbol(runtyp),
                                theta0_dict=th0d,
                                fit_target=fit_target)
            pmcmc = BM.MCMCNum(;pmcmc_defaults...)

            # # Make not continuous glacier continuous:
            # if isinf(logprior_debug(theta0.th0,false)[2][2])
            #     theta0.th0[theta0.names[:btilde]] += 1
            # end

            if debugging
                @inferred logposterior(theta0.th0)
                print("Time for logposterior:  ")
                @time logposterior(theta0.th0)
                println("Value logposterior:  $(logposterior(theta0.th0)[1])")
            end


            res = if debugging
                print("Time for $(pmcmc.niter) MCMC steps:  ")
                @time BM.mcmc(logposterior, theta0, pmcmc;
                              verbose=false, use_progress_meter=false)
            else
                BM.mcmc(logposterior, theta0, pmcmc;
                        verbose=false, use_progress_meter=false)
            end
            # # h and iv averaged to where parameters are
            # h_th, iv_th = BM.h_iv_at_thetas(theta0, mc_hs1d, mc_ivs1d)

            # theta0_h_iv = BM.Theta0(theta0_with_h_iv_dict, gb, pp, pm, pn)
            # theta0_h = BM.Theta0(theta0_with_h_dict, gb, pp, pm, pn)
            # theta0_iv = BM.Theta0(theta0_with_iv_dict, gb, pp, pm, pn)

            # varnames = BM.get_varnames(theta0)
            # varnames_h_iv = BM.get_varnames(theta0_h_iv)
            # varnames_h = BM.get_varnames(theta0_h)
            # varnames_iv = BM.get_varnames(theta0_iv)

            # thetas_h_iv = vcat(thetas, h_th, iv_th)
            # thetas_h = vcat(thetas, h_th)
            # thetas_iv = vcat(thetas, iv_th)

            # #print_results(thetas, accept_ratio, names=varnames)

            # fig = BM.plottheta(thetas, theta0, toplot=:all, size=(1400,1400))
            # Plots.savefig(fig, fl*"--thetas.png")

            # fig = BM.plottheta(thetas_h, theta0_h, toplot=:all, size=(1400,1400))
            # Plots.savefig(fig, fl*"--thetas_h.png")

            # fig = BM.plottheta(thetas_iv, theta0_iv, toplot=:all, size=(1400,1400))
            # Plots.savefig(fig, fl*"--thetas_iv.png")

            # # more plots
            # # 1D results with bands
            # fig = BM.plotinv1d(gb, blobs)
            # Plots.savefig(fig, fl*"--1d.png")

            # # 2D error binned back to 1D
            # fig = BM.plotinv1d_err(gb, blobs)
            # Plots.savefig(fig, fl*"--1d_err.png")

            # # lots of 2D plots
            # fig = BM.plotinv2d(gl, blobs, size=(1400,1400))
            # Plots.savefig(fig, fl*"--2d_overview.png")

            # # error of mean IV
            # fig = BM.plotinv2d_iverr(gl, blobs, reuse=false)
            # Plots.savefig(fig, fl*"--2d_iverr.png")

            # # summaries
            # # as JLD
            # res = KissMCMC.summarize_run(thetas_h_iv, names=varnames_h_iv)
            # FileIO.save(fl*".jld2", "res", res)

            # # as text
            # fll = splitdir(fl)[2]
            # ar = hcat([[strip(string(v)),m,s] for (v,m,s) in zip(res[:var],res[:mean],res[:std])]...)
            # l1,l2,l3 = split(reprmime("text/plain", ar),'\n')[2:4]
            # write(fl*"--summary", "$fll\n    $l1\n    $l2\n    $l3\n")

            # ## forward
            # ## run forward model with most likely paras
            # res = KissMCMC.summarize_run(thetas, names=varnames_h_iv)

            # hs1d, ivs1d, hs2d, ivs2d, pm_,
            # taus1d, taus1d_l, taus2d = fwdm_fn(res[:mode])
            # fig = BM.plot_run(gb, hs1d, taus1d, ivs1d, hs2d, taus2d,
            #                   ivs2d, taus1d_l, pp, pm,
            #                   plot2d=[:h,:iv,:ivlog][1])
            # Plots.savefig(fig, fl*"--map.png")

            # FileIO.save(fl*"-map2d.jld2",
            #             "mh2d", mh2d, "sh2d", sh2d,
            #             "miv2d", miv2d, "siv2d", siv2d,
            #             "min_h2d", min_h2d, "max_h2d", max_h2d,
            #             "min_iv2d", min_iv2d, "max_iv2d", max_iv2d)

        @time bson(fl*".bson", Dict(:res => res)) # takes 2 minutes!
        catch e
            e==InterruptException() && rethrow(e)
            fll = splitdir(fl)[2]
            write(fl*"--summary", "$fll\n    error $e\n")
            println("error: $e")
        end
        if debugging
            toc()
            println("\n")
        end
    end
end
