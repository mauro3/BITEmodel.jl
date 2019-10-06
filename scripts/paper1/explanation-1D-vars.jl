using Parameters, Plots, VAWTools
include("../itmix-setup.jl")

if !isdefined(:gl_all_radar)
    glacier = :Unteraar
    gid = BM.ITMIXGlacier(glacier, 2)
    gl_all_radar,gb,pp,pm,pn,pl = BM.init_forward(gid)
end

pm = BM.MPara(pm,
              sigma_h_model=20,
              error_model=[BM.ErrorModel.sqrerr, BM.ErrorModel.abserr][1],
              nthin_h=1)

runtyp = [:test,:testmid,:prodlow,:prodmid,:prod][3]
fit_target = [BM.FitTarget.h, BM.FitTarget.iv, BM.FitTarget.h_iv, BM.FitTarget.length][1]
n_1d_btilde = 3
n_1d_fsl = 3
n_1d_temp = 1
rm_iv = false
r_tracks = [5,34,42]

(gln, gbn, th0d,
 theta0, logposterior, logposterior1d,
 logprior, logprior_debug, loglikelihood, loglikelihood1d,
 fwdm_fn, fwdm1d_fn, pmcmc_defaults, fit_target) =
     setup_one_exp(r_tracks, gl_all_radar, gb, pp, pm, pn,
                   n_1d_btilde, n_1d_fsl, n_1d_temp,
                   runtyp, fit_target, rm_iv)

theta = zeros(theta0.th0)
theta[1:3] = [-2.5, 0.7, 0.0]

xtheta = VAWTools.piecewiselinear(gb.ele, gb.xmid).(theta0.th_ele[:btilde])


new_btilde = BM.make_fields1d(theta, theta0)[1]

plot(gb.xmid/1e3, gb.bdot-gb.dhdt, xlabel="χ (km)", ylabel="b̃ (ma⁻¹)", label="Measurements", size=(400,300))
plot!(xtheta/1e3, theta[1:3], marker=:circle, label="θ")
plot!(gb.xmid/1e3, new_btilde, label="Used")

#savefig("figs/explanation-1D-vars.png")
#savefig("figs/explanation-1D-vars.pdf")
