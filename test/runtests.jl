using BITEModel
using Base.Test
const BM = BITEModel

## Forward model
################

plotyes = false
verbose = true

global gid, hs1d, ivs1d, hs2d, ivs2d, glaciertype, # this is a bit ugly but
       gl, gb, pp, pm, pn                          # needed because of
                                                   # include("../scripts/forward.jl")

@testset "Forward model with update_cache=$update_cache" for update_cache in [false,true]
    # Synthetic (run by including scripts/forward.jl)
    ###########
    glaciertype = :synth
    plotyes = false
    include("../scripts/forward.jl")
    Base.Test.@inferred BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
    Base.Test.@inferred BM.make_bands(gl, pp, pm, pn)
    fwdsol = BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)

    hs1d_ = [0.0, 167.025, 182.466, 189.223, 188.814, 183.506, 175.65,
             166.279, 151.378, 141.485, 131.018, 119.681, 109.229, 98.3605,
             86.7481, 76.5545, 67.8539, 63.5029, 57.3474, 50.484, 1.0]
    ivs1d_ = [0.0, 3.3608, 5.31713, 6.91192, 8.26612, 9.03612, 9.67009,
              10.0828, 10.5786, 11.0246, 10.9542, 10.8031, 10.8764, 10.6826,
              9.68466, 8.8367, 7.85266, 6.36658, 4.69378, 2.97981, 0.0]
    @test norm(hs1d_-fwdsol.hs1d)<2e-3
    @test norm(ivs1d_-fwdsol.ivs1d)<5e-4

    # test 2D at some random indices
    @test length(fwdsol.hs2d)==15921
    inds = [8179, 15516, 3523, 5846, 3236, 13314, 4959, 11, 12505, 6528,
            11823, 13553, 1043, 1420, 4210, 1991, 34, 4056, 5848, 15447,
            13452, 390, 8071, 9023, 12655, 3717, 7865, 8238, 14551, 4782]
    hs2d_ = [176.314, 0.0, 129.79, 172.317, 108.296, 0.0, 0.0, 0.0,
             57.2049, 0.0, 104.003, 0.0, 0.0, 0.0, 55.4563, 0.0, 0.0, 161.232,
             174.862, 0.0, 109.54, 0.0, 95.0336, 222.902, 119.131, 82.9983,
             98.8966, 230.184, 0.0, 132.357]
    ivs2d_ = [8.22742, 0.0, 6.17785, 7.93149, 6.84787, 0.0, 0.0, 0.0, 0.0,
              0.0, 7.88523, 0.0, 0.0, 0.0, 7.28846, 0.0, 0.0, 5.7638,
              7.93149, 0.0, 4.89699, 0.0, 0.0, 6.58108, 6.01924,
              7.96476, 6.96876, 6.70248, 0.0, 8.22056]
    @test norm(hs2d_-fwdsol.hs2d[inds])<1e-3
    @test norm(ivs2d_-fwdsol.ivs2d[inds])<5e-4

    # # ITMIX (run by hand)
    # #####################
    # glacier=:Urumqi # pretty fast this one
    # pl_kws = Dict{Symbol,Any}(:use_glogem => false)
    # gid = ITMIXGlacier(glacier)
    # gl,gb,pp,pm,pn,pl = BM.init_forward(gid, verbose;
    #                                     update_cache=update_cache, pl_kws...)
    # hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb =
    #     BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
    # Base.Test.@inferred BM.fwdm(gb, pp, pm, pn, gb.bdot-gb.dhdt, gb.fsl, gb.temp)
    # Base.Test.@inferred BM.make_bands(gl, pp, pm, pn)

    # hs1d_ = [0.0, 25.4951, 29.2112, 32.2141, 32.3665, 30.4775, 29.6186,
    #          29.2802, 27.4351, 28.3908, 29.75, 37.2109, 39.2642,
    #          40.3066, 42.8661, 44.964 , 48.2977, 53.8874, 56.7994,
    #          57.8052, 51.7617, 40.3284, 31.8477, 5.63301]
    # ivs1d_ = [0.0, 1.47265, 2.2546, 3.01979, 3.92654, 4.90726,
    #           4.87954, 4.17233, 3.82581, 3.51174, 3.2171, 2.61505, 2.57666,
    #           2.68062, 2.84157, 3.30573, 3.93735, 4.08862, 4.08278, 4.08057,
    #           4.2603, 4.31557, 2.45511, 0.00525341]
    # @test norm(hs1d_-hs1d)<5e-4
    # @test norm(ivs1d-ivs1d_)<5e-4

    # # test 2D at some random indices
    # @test length(hs2d)==133331
    # inds = [85082, 118470, 7808, 31087, 78719, 96341, 104214, 9715,
    #         91721, 131435, 49463, 127799, 36643, 131306, 71573, 88897,
    #         110465, 58186, 93094, 29907, 105007, 65361, 118730, 67826,
    #         118921, 66973, 114737, 84487, 19245, 126314]
    # hs2d_ = [0.0, 54.1206, 0.0, 0.0, 39.7876, 0.0, 0.0, 39.0787,
    #          23.5233, 0.0, 0.0, 19.2427, 0.0, 0.0, 22.237, 50.1744,
    #          0.0, 71.186, 0.0, 35.3194, 27.0517, 77.9176, 0.0, 21.425,
    #          0.0, 0.0, 62.8376, 82.3543, 43.1886, 28.9674]
    # ivs2d_ = [0.0, 1.73392, 0.0, 0.0, 3.50925, 0.0, 0.0, 1.3902,
    #           5.06693, 0.0, 0.0, 3.04547, 0.0, 0.0, 2.57997, 3.94191,
    #           0.0, 2.76635, 0.0, 1.6151, 6.13513, 3.14324, 0.0,
    #           2.31482, 0.0, 0.0, 1.87241, 4.19182, 1.40529, 1.68439]
    # @test norm(hs2d_-hs2d[inds])<5e-4
    # @test norm(ivs2d_-ivs2d[inds])<5e-4
end

# make sure this runs
include("../scripts/forward-advanced.jl")

## Inverse model
################
@testset "Inverse model" begin
    gid = SyntheticGlacier(:bench)
    gl,gb,pp,pm,pn,pl = BM.init_forward(gid)
    (theta0, logposterior, logposterior1d, logprior, logprior_debug,
     loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
     pmcmc_defaults, fit_target) =
         BM.init_inverse(gb, pp, pm, pn)
    Base.Test.@inferred fwdm_fn(theta0.th0)
    Base.Test.@inferred logprior(theta0.th0)
    Base.Test.@inferred loglikelihood(theta0.th0)
    Base.Test.@inferred logposterior(theta0.th0)

    pmcmc = BM.MCMCNum(;pmcmc_defaults...,
                       niter=5*10^4,
                       nthin = 50)

    # sample posterior
    srand(0)
    thetas, theta_expect, theta_mode, accept_ratio, blobs, lpost, mc_hs1d, mc_ivs1d,
      mh2d, sh2d, miv2d, siv2d, min_h2d, max_h2d, min_iv2d, max_iv2d = mcmc(logposterior, theta0, pmcmc; verbose=true)
    @test round(accept_ratio,5) == 0.11801

    @test norm(theta_expect- [1.81352, 1.40765, -0.29361, -0.170199,
                              -0.0292185, -0.15876, -0.3106, -0.165316, -0.0345219, 0.172756,
                              0.406927, 0.162265, -0.144996, 0.315479, 0.16641, -0.00986009,
                              0.0872772, -0.109316, 0.711055, -0.239801, -0.121702, -0.0558777,
                              -0.219639]) < 1e-5
    # test 2D at some random indices
    inds = [8179, 15516, 3523, 5846, 3236, 13314, 4959, 11, 12505, 6528,
            11823, 13553, 1043, 1420, 4210, 1991, 34, 4056, 5848, 15447,
            13452, 390, 8071, 9023, 12655, 3717, 7865, 8238, 14551, 4782]

    @test norm(mh2d[inds]-[195.185, 0.0, 112.593, 175.143, 110.134, 0.0, 0.0, 0.0, 48.5617,
                           0.0, 115.334, 0.0, 0.0, 0.0, 48.79, 0.0, 0.0, 127.909, 173.234,
                           0.0, 86.9006, 0.0, 80.6748, 169.704, 104.964, 85.447, 86.3577,
                           180.154, 0.0, 147.823])<2e-3

    @test norm(sh2d[inds]-[14.5854, 0.0, 21.0239, 18.2294, 11.4231, 0.0, 0.0, 0.0, 9.05472,
                           0.0, 10.1193, 0.0, 0.0, 0.0, 19.4267, 0.0, 0.0, 24.128, 21.1017,
                           0.0, 16.3925, 0.0, 15.0425, 27.2756, 19.6065, 9.78661, 33.134,
                           32.304, 0.0, 11.914])<2e-4

    @test norm(miv2d[inds]-[11.8817, 0.0, 8.42797, 10.3403, 10.106,
                            0.0, 0.0, 0.0, 0.0, 0.0, 12.5615, 0.0, 0.0, 0.0, 16.4126, 0.0,
                            0.0, 7.73834, 10.3403, 0.0, 7.61675, 0.0, 0.0, 7.85299, 8.40403,
                            14.5605, 16.3398, 8.09658, 0.0, 12.651])<1e-4
end
