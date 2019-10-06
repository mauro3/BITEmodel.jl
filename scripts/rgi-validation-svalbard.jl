# RGI-based runs for validation

using BSON, VAWTools, DataStructures
region = [3, 4, 14, 11, 7][end-1]
runtyp = ["test", "testmid", "prodlow", "prod"][3]
runtyp == :test && println("\n\nTEST RUN !!!!!!!!\n\n")
parallel = true
all_glaciers = parallel
save_geotiff = true
save_sols = false # makes huge files and did not work with NamedTuples...
fit_sigma = false # not implemented yet

today = Dates.today()
region_str = VAWTools.int2str2(region)

if runtyp in ["test", "testmid"]
    dir_output = "output/validation-test/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))/"
else
    if fit_sigma
        dir_output = "output/validation-$(region_str)_fit-sigma/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))/"
    else
        dir_output = "output/validation-$(region_str)/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))/"
    end
end
!isdir(dir_output) && mkpath(dir_output)


@show region, runtyp, parallel, dir_output

include("rgi-setup.jl") # makes sure all is pre-compiles

tmp = BM.get_git_shas(error_on_dirty=false)
description = """
    bugfix with expectation value

    fsl uniform and 0.5

    Landmask for adjacent glacier now activated

    Validation on Svalbard

    Now after some major bug-fix in expectation-prior update calculation.

    BITEModel commit $(tmp["BITEModel"][1])
              branch $(tmp["BITEModel"][2])
    """

if parallel && nprocs()<2
    # if Sys.CPU_CORES>=17

    #     addprocs(17) # /2 to get to physical cores
    # else
    #     addprocs(Sys.CPU_CORES÷2 + 1) # /2 to get to physical cores
    # end
    # addprocs(Sys.CPU_CORES÷5) # five regions
    addprocs(Sys.CPU_CORES)
end

@everywhere begin
    const debugging = true
    include("rgi-setup.jl")
    using NamedTuples
    const progress_bar = true
end

update_cache = true
glnrs_with_radar = get_glaciers(region)[1] #[8:end] # drop the first few, too big
# delete the ones with a bad DEM
if region==7
    for x in [00025, 00028, 00037, 00042, 00228, 00246, 00279, 00409, 00474, 00529, 00531, 00570, 01425, 01466, 01481, 01523, 01559]
        i = findfirst(glnrs_with_radar, x)
        if i>0
            deleteat!(glnrs_with_radar, i)
        end
    end

    # also delete the ones which error (none found currently)

    # Svalbard glaciers with dhdt data as hand-selected in qgis
    svalbard_gls_with_dhdt = [133, 168, 170, 174, 175, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196,
                              197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 210, 213, 225, 266, 267, 268, 269, 270, 271, 272, 273, 274,
                              275, 276, 277, 278, 279, 280, 281, 284, 293, 294, 296, 313, 314, 315, 316, 325, 326, 327, 328, 329, 330, 331, 332, 333,
                              334, 336, 337, 338, 339, 340, 341, 342, 343, 344, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358,
                              359, 360, 361, 362, 363, 364, 365, 366, 367, 368, 369, 370, 371, 372, 373, 374, 375, 376, 377, 378, 379, 380, 381,
                              382, 383, 384, 385, 386, 388, 389, 390, 391, 392, 393, 394, 395, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406,
                              407, 408, 409, 410, 411, 412, 413, 414, 415, 416, 417, 418, 419, 420, 421, 426, 427, 428, 429, 430, 431, 432, 433,
                              434, 435, 436, 437, 1403, 1414, 1415, 1461, 1467, 1468, 1469, 1470, 1475, 1476, 1478, 1480, 2, 3, 4, 6, 8, 9, 10, 11,
                              12, 13, 14, 15, 16, 17, 18, 19, 21, 22, 24, 86, 87, 88, 89, 90, 91, 93, 95, 96, 97, 100, 101, 102, 103, 104, 105, 106,
                              107, 108, 109, 112, 113, 114, 115, 116, 118, 119, 120, 121, 124, 125, 126, 128, 129, 130, 131, 132, 134, 135, 136, 137,
                              138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
                              161, 162, 163, 164, 165, 166, 167, 169, 172, 173, 176, 177, 179, 208, 209, 226, 227, 228, 229, 230, 231, 232, 233, 234
                              , 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 252, 253, 254, 255, 256, 257, 258, 259, 260,
                              261, 262, 263, 284, 295, 297, 301, 441, 442, 443, 446, 1402, 1406, 1419, 1420, 1421, 1422, 1423, 1438, 1439]
    # 39 have both:
    svalbard_gls_with_dhdt_and_radar = [s for s in svalbard_gls_with_dhdt if s in glnrs_with_radar]
elseif region==11
    # in the Alps not all have IV
    alps_gl_with_iv_n_radar = [1450, 1478, 1328, 2766, 1238, 2739, 1198, 2773, 2630, 1275, 2801, 2709, 1698, 2704, 2715, 2624, 1144, 2596,
                               2793, 2890, 2740, 2507, 2787, 2774, 2746, 1876, 2490, 1791, 1509, 2584, 2006, 2558, 2249, 2755, 1576, 1296, 2884, 2679, 1987, 1986, 1067,
                               2869, 2909, 1786, 1857, 2600, 1928, 2448, 2673, 1492, 2583, 2244, 1894, 2549, 2027, 1199, 2745, 1367, 1344, 1376]
end

glnrs_fit = glnrs_with_radar[1:2:end]
glnrs_test = glnrs_with_radar[2:2:end]
if length(glnrs_fit)>length(glnrs_test)
    # set to same length
    glnrs_fit = glnrs_fit[1:end-1]
end

# To reduce to just a few glacier
if !all_glaciers
    st = [17, 34, 100][2] # 5, 3, 1
    glnrs_fit = glnrs_fit[1:st:end]
    glnrs_test = glnrs_test[1:st:end]
else
    # st = 3
    # glnrs_fit = glnrs_fit[1:st:end][1:20]
    # glnrs_test = glnrs_test[1:st:end][1:20]
end

if region==7
    # 22:
    glnrs_fit_with_dhdt = [s for s in glnrs_fit if s in svalbard_gls_with_dhdt_and_radar]
    # 17:
    glnrs_test_with_dhdt = [s for s in glnrs_test if s in svalbard_gls_with_dhdt_and_radar]
elseif region==11
    # 30 each
    glnrs_fit_with_iv = [s for s in glnrs_fit if s in alps_gl_with_iv_n_radar]
    glnrs_test_with_iv = [s for s in glnrs_test if s in alps_gl_with_iv_n_radar]
end

## MCMC

# parameters to play with
sigma_of_model = @NT(sigma_h_model = 20.0,
                     sigma_iv_model = 20.0)


# number of divisions in elevation band variables
fit_vars = @NT(n_1d_btilde=3,
               n_1d_fsl=3,
               n_1d_temp=1)

@show region, debugging, runtyp, length(glnrs_fit)

# Method errors contain their arguments, which in our case can be huge.
# Remove them in case of error in pmap
@everywhere strip_MethodErrors(err) = err isa MethodError ? MethodError(0, err.f, err.world) : err

# Write description of model run
write(dir_output*"description.txt", description)

# Full fitting Mode 1
#####################

tic = time()

# runtyp=:prodlow
# nr = 2819
# out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(1), dir_output="", retsol=true)
# error()

# out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(1), dir_output="", retsol=true)
# out_fit = map(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
#                         fit_target, runtyp, Val(1)), glnrs_fit)

expect_mean = Dict()
expect_median = Dict()
stdev = Dict()

out_ = 1 # to get "debug" access to below calculation in case of errors
for fit_target = [BM.FitTarget.h, BM.FitTarget.iv, BM.FitTarget.h_iv, BM.FitTarget.length]

    glnrs_fit_ = if region==11 && fit_target==BM.FitTarget.iv
        #  for the Alps only use the ones with IV
        glnrs_fit_with_iv
    else
        glnrs_fit
    end

    out_ = pmap(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
                                fit_target, runtyp, Val(1),
                                dir_output=dir_output, save_geotiff=save_geotiff,
                                save_sol=false, retsol=save_sols, flroot="FIT_$(fit_target)",
                                fit_sigma=fit_sigma,
                                ),
                glnrs_fit_, on_error=strip_MethodErrors, retry_delays = zeros(3));
    out_fit = save_sols ? [ (o isa Exception) ? o : o[1] for o in out_] : out_;

    @assert length(out_fit)>0 "All fitting glaciers errored!"

    toc1 = time() - tic
    println("Iteration over fit glaciers $fit_target completed in $(toc1/60)min")

    # Make the priors of the fitted parameters
    ##########################################
    drop_sigma = true
    _, expect_median_, expect_mean_, stdev_ = calculate_prior_update(out_fit, drop_sigma)
    expect_mean[Symbol(fit_target)] = expect_mean_
    expect_median[Symbol(fit_target)] = expect_median_
    stdev[Symbol(fit_target)] = stdev_

    ## store
    fl = joinpath(dir_output, "reg$(region_str)_FIT_$(fit_target)")
    println("Storing $fl.bson ...")
    @time bson(fl*".bson", Dict(:out_fit => deepcopy(out_fit),
                                :description => description,
                                :commit => BM.get_git_shas(error_on_dirty=false),
                                :glnrs_fit => glnrs_fit_,
                                :glnrs_test => glnrs_test,
                                :runtyp => runtyp,
                                :expect_median => expect_median_,
                                :expect_mean => expect_mean_,
                                :stdev => stdev_,
                                ))

    println("... done.")
    if save_sols
        sols = [o[2] for o in out_ if !(o isa Exception)]

        println("Storing $fl.jls ...")
        sols_fit = OrderedDict(BM.getrgi_nr(o[:rgi])=>s for (o,s) in zip(out_fit,sols))
        open(fl*".jls", "w") do io
            serialize(io, sols_fit)
        end
        println("... done.")
    end

end
#  Base case: no change in prior
ln = length(first(expect_median)[2])
expect_median[:base] = zeros(ln)
expect_mean[:base] = zeros(ln)
stdev[:base] = Inf.*ones(ln)

fl = joinpath(dir_output, "reg$(region_str)_priors")
bson(fl*".bson", Dict(:expect_median => deepcopy(expect_median),
                      :expect_mean => deepcopy(expect_mean),
                      :stdev => deepcopy(stdev),
                      ))

# keys(expect_mean) == [:base, :h, :h_iv, :iv, :length]

# Validation: posterior1d only fitting
############################


tic2 = time()
# validation using updated priors both against test-set and fit-set (the latter is for double-checking)
for (fl_str, glnrs_) in zip(["testgls", "fitgls"], [glnrs_test, glnrs_fit])
    for k in keys(expect_median)
        fit_target = BM.FitTarget.length
        expect = expect_median[k]
        stdev_ = stdev[k]
        if k!=:base
            @assert expect!=nothing
            @assert stdev_!=nothing
        end

        #out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(2), expect, stdev_, dir_output=dir_output, retsol=true)
        out_ = pmap(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
                                 fit_target, runtyp, Val(2),
                                 expect, stdev_,
                                 dir_output=dir_output, save_geotiff=save_geotiff,
                                 flroot="PRIOR_$k-$fl_str",
                                 fit_sigma=false,
                                 save_sol=false, retsol=save_sols),
                    glnrs_, on_error=strip_MethodErrors, retry_delays = zeros(3))
        out_test = save_sols ? [ (o isa Exception) ? o : o[1] for o in out_] : out_
        @assert length(out_test)>0 "All testing glaciers errored!"

        toc2 = time() - tic2

        println("Iteration over validation glaciers $k completed in $(toc2/60)min")


        fl = joinpath(dir_output, "reg$(region_str)_PRIOR_$k-$fl_str")
        println("Storing $fl.bson ...")
        @time bson(fl*".bson", Dict(:out_test => deepcopy(out_test),
                                    :description => description,
                                    :commit => BM.get_git_shas(error_on_dirty=false),
                                    :expect => expect,
                                    :stdev => stdev_,
                                    ))
        println("... done.")
        if save_sols
            println("Storing $fl.jls ...")
            sols = [o[2] for o in out_ if !(o isa Exception)]
            sols_test = OrderedDict(BM.getrgi_nr(o[:rgi])=>s for (o,s) in zip(out_fit,sols))
            open(fl*".jls", "w") do io
                serialize(io, sols_test)
            end
            println("... done.")
        end
    end
end

#####
# Last IV fitting with priors from h_iv and base (all else doesn't make sense)
fit_target = BM.FitTarget.iv
for k in [:base, :h_iv]
    # Here we need to add back the stripped sigma_iv term.  With stdev==Inf it
    # does not modify the prior.
    expect = [expect_median[k]..., 0]
    stdev_ = [stdev[k]..., Inf]

    glnrs_test_ = if region==11 && fit_target==BM.FitTarget.iv
        #  for the Alps only use the ones with IV
        glnrs_test_with_iv
    else
        glnrs_test
    end

    out_ = pmap(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
                             fit_target, runtyp, Val(1),
                             expect, stdev_,
                             dir_output=dir_output, save_geotiff=save_geotiff,
                             flroot="PRIOR_$(k)_FIT_iv",
                             fit_sigma=fit_sigma,
                             save_sol=false, retsol=save_sols),
                glnrs_test_, on_error=strip_MethodErrors, retry_delays = zeros(3))
    out_test = save_sols ? [ (o isa Exception) ? o : o[1] for o in out_] : out_
    @assert length(out_test)>0 "All IV-testing glaciers errored!"

    toc1 = time() - tic

    fl = joinpath(dir_output, "reg$(region_str)_PRIOR_$(k)_FIT_$(fit_target)")
    println("Storing $fl.bson ...")
    @time bson(fl*".bson", Dict(:out_test => deepcopy(out_test),
                                :description => description,
                                :commit => BM.get_git_shas(error_on_dirty=false),
                                :glnrs_fit => glnrs_fit,
                                :glnrs_test => glnrs_test_,
                                :runtyp => runtyp,
                                :expect => expect,
                                :stdev => stdev_,
                                ))

    println("... done.")
    if save_sols
        sols = [o[2] for o in out_ if !(o isa Exception)]
        println("Storing $fl.jls ...")
        sols_test_iv = OrderedDict(BM.getrgi_nr(o[:rgi])=>s for (o,s) in zip(out_fit,sols))
        open(fl*".jls", "w") do io
            serialize(io, sols_test_iv)
        end
        println("... done.")
    end
end

toc = time() - tic
println("Total time $(toc/60/60) hours")
