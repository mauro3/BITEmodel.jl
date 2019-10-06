## Loading of inverse-model stuff

"""
Returns a default theta0 initialization dict.
"""
function theta0_dict_defaults(pm, n_1d_btilde, n_1d_fsl, n_1d_temp; with_sigma=false)
    out = OrderedDict(:btilde=> zeros(F,n_1d_btilde), # NOTE in m/y (error)
                      :fsl=> zeros(F,n_1d_fsl), # error
                      :temp=> zeros(F,n_1d_temp), # error
                      :dist_exp=> pm.dist_exp,
                       :iv_h_exp=> pm.iv_h_exp,
                      :iv_dist_exp1=> pm.iv_dist_exp1,
                      :iv_dist_exp2=> pm.iv_dist_exp2,
                      )
    if with_sigma
        out[:sigma_h_model] = pm.sigma_h_model
        out[:sigma_iv_model] = pm.sigma_iv_model
    end
    return out
end

"""
The default ball radius to be used in the emcee sampler initial conditions.
The size is chosen to be in line with the range of allowed values in priors.jl.
"""
const ball_radius_default = Dict(:btilde=> 1.4, # corresponds to sigma of 1 for dhdt and 1 for bdot
                                 :fsl=> 0.5,
                                 :temp=> 15.0,
                                 :dist_exp=> 0.8,
                                 :iv_h_exp=> 15.0, # @assert iv_h_exp>=0
                                 :iv_dist_exp1=> 0.4, # @assert 0<=iv_dist_exp1<=1
                                 :iv_dist_exp2=> 15.0, # @assert iv_dist_exp2>0
                                 :sigma_h_model => 100.0,
                                 :sigma_iv_model => 100.0,
                                 )
"""
    pmcmc_dict_defaults(theta0, runtyp=[:test,:mid,:prod][1])

Return, hopefully, sensible defaults for the construction of MCMCNum type with

    MCMCNum(pmcmc_dict_defaults(theta0, runtyp)...; other_kws)

The second argument allows to choose between quick-test, slower-test and production run
settings.
"""
function pmcmc_dict_defaults(theta0, runtyp)
    out = Dict(:sigma_ppdf => 0.001 + zeros(theta0.th0),
               :ball_radius => [ball_radius_default[name] for (name, inds) in theta0.names for i in inds]
               )
    dof = length(theta0.th0)
    return if runtyp==:test
        nchains = dof+2
        merge(out, Dict(:niter=>11*nchains,
                        :nthin => 1,
                        :nburnin=>10,
                        :nchains=>nchains))
    elseif runtyp==:testmid
        nchains = max(2*(dof+2), 30)
        merge(out, Dict(:niter=>10^4,
                        :nburnin=>9000,
                        :nchains=>nchains))
    elseif runtyp==:prodlow
        nchains = max(2*(dof+2), 50)
        merge(out, Dict(:niter=>10^5,
                        :nthin => 20,
                        :nburnin => 7*10^4,
                        :nchains=>nchains))
    elseif runtyp==:prodmid
        nchains = max(2*(dof+2), 50)
        merge(out, Dict(:niter=>5*10^5,
                        :nthin => 20,
                        :nburnin => 4*10^5,
                        :nchains=>nchains))
    elseif runtyp==:prodhigh
        nchains = min(10*dof, max(dof+2, 50))
        merge(out, Dict(:niter=>10^6,
                        :nchains=>nchains,
                        :nthin => 50,
                        :nburnin => 7*10^5))
    elseif runtyp==:prod
        nchains = min(10*dof, max(dof+2, 50))
        merge(out, Dict(:niter=>3*10^5,
                        :nchains=>nchains,
                        :nthin => 20,
                        :nburnin => 2*10^5))
    else
        error("Unknown option `runtyp`==$runtyp")
    end
end


"""
$(SIGNATURES)

Initialize all inverse-model parameters.

Returns
- theta0
- logposterior
- logprior, logprior_debug
- loglikelihood
- fwdm_fn
- pmcmc-dict-defaults
- fit_target_eff: what run-type was effectively chosen (depending on data availability)
"""
init_inverse(gb::Bands, pp, pm, pn; kwargs...) =
    init_inverse(gb.gl.gid, gb, pp, pm, pn; kwargs...)

function init_inverse(gid::GlacierID, gb::Bands, pp, pm, pn;
                      runtyp=[:test,:mid,:prod][1],
                      n_1d_btilde = 3, n_1d_fsl = 3, n_1d_temp = 1,
                      fit_target=FitTarget.h_iv,
                      theta0_dict=theta0_dict_defaults(pm, n_1d_btilde, n_1d_fsl, n_1d_temp),
                      check_pdf=true, # if true check and try to update if -Inf
                      do_grid_serach=true
                      )
    theta0_dict = deepcopy(theta0_dict)
    if pm.iv_nye
        pop!(theta0_dict, :iv_h_exp, nothing)
        pop!(theta0_dict, :iv_dist_exp1, nothing)
    else
        pop!(theta0_dict, :iv_dist_exp2, nothing)
    end

    ## Modify fit_target:
    if fit_target==FitTarget.h_iv
        # - if there is no IV data, then drop iv-part
        if gb.gl.iv isa BITEModel.SomeData{Void}
            fit_target = FitTarget.h
        end
        # - if there is no h data, then drop h-part
        if gb.gl.h isa BITEModel.SomeData{Void}
            fit_target = FitTarget.iv
        end
    end
    #  - if there is no data to fit to, the only fit the glacier length
    if  (fit_target==FitTarget.h && gb.gl.h isa BITEModel.SomeData{Void}) ||
        (fit_target==FitTarget.iv && gb.gl.iv isa BITEModel.SomeData{Void})
        fit_target = FitTarget.length
    end

    if fit_target==FitTarget.length || fit_target==FitTarget.prior
        pop!(theta0_dict, :sigma_h_model, nothing)
        pop!(theta0_dict, :sigma_iv_model, nothing)
    elseif fit_target==FitTarget.iv
        pop!(theta0_dict, :sigma_h_model, nothing)
    elseif fit_target==FitTarget.h
        pop!(theta0_dict, :sigma_iv_model, nothing)
    elseif  fit_target==FitTarget.h_iv
        nothing
    else
        error("Unknow $(fit_target)")
    end

    theta0 = Theta0(theta0_dict, gb, pp, pm, pn)

    logprior, logprior_debug = make_logprior(theta0)
    pn.test && Base.Test.@inferred logprior(theta0.th0)

    # forward model
    loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn = make_loglikelihood(theta0, pp, pm, pn, fit_target)
    pn.test && Base.Test.@inferred fwdm_fn(theta0.th0)
    pn.test && Base.Test.@inferred loglikelihood(theta0.th0)
    pn.test && Base.Test.@inferred loglikelihood1d(theta0.th0)

    logposterior, logposterior1d = make_logposterior(loglikelihood, loglikelihood1d, logprior)
    pn.test && Base.Test.@inferred logposterior(theta0.th0)

    # Try hard to get a btilde which works decently
    if check_pdf && logposterior(theta0.th0)[1]==-Inf || do_grid_serach
        orig = theta0_dict[:btilde]

        # Note: only fit length (through prior for calving glaciers and through 1d likelihood otherwise)

        if length(orig)!=3
            # Do an abbreviated grid search
            best = -Inf
            delta_btilde_best = zeros(orig)
            for delta_btilde = logspace(0, 1.0, 10) - 0.95
                for plusminus = [-1,1] # add or subtract delta_btilde
                    for ind=0:2^(length(orig))-1
                        # index vector
                        vec = [parse.(Bool, c) for c in bin(ind, length(orig))]
                        theta0_dict[:btilde] = orig + plusminus * delta_btilde .* vec
                        theta0 = Theta0(theta0_dict, gb, pp, pm, pn)
                        lp = logposterior1d(theta0.th0)[1]
                        if lp>best
                            best = lp
                            delta_btilde_best .= plusminus * delta_btilde .* vec
                        end
                    end
                end
            end
        else
            # Do a grid search
            best = -Inf
            delta_btilde_best = [0.0, 0.0, 0.0]
            nn = 5
            delta_btildes = [-(logspace(0, 1.0, nn) - 0.95); 0 ; logspace(0, 1.0, nn) - 0.95]
            for di in delta_btildes
                for dj in delta_btildes
                    for dk in delta_btildes
                        theta0_dict[:btilde] = orig .+ [di, dj, dk]
                        theta0 = Theta0(theta0_dict, gb, pp, pm, pn)
                        lp = logposterior1d(theta0.th0)[1]
                        if lp>best
                            best = lp
                            delta_btilde_best = [di, dj, dk]
                        end
                    end
                end
            end
        end
        # use the delta_btilde_best
        theta0_dict[:btilde] = orig .+ delta_btilde_best
        theta0 = Theta0(theta0_dict, gb, pp, pm, pn)

        @assert logposterior(theta0.th0)[1]>-Inf "Could not find a theta0 such that logposterior>-Inf"
    end

    return (theta0, logposterior, logposterior1d, logprior, logprior_debug,
            loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
            pmcmc_dict_defaults(theta0, runtyp), fit_target)
end
