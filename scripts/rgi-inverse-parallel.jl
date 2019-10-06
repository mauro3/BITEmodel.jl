# RGI-based runs for ESA

using BSON, VAWTools
region = [3, 4, 14, 11, 7][3]
runtyp = ["test", "testmid", "prodlow", "prod"][3]
runtyp == :test && println("\n\nTEST RUN !!!!!!!!\n\n")
parallel = true
all_glaciers = parallel
repeat_mode_1 = 1 # number of mode-1 repetitions to be able to calculate Rhat
repeat_mode_2 = 1 # number of mode-2 repetitions to be able to calculate Rhat

fit_sigma = true
use_synthetic_dhdt = true

today = Dates.today()
if runtyp in ["test", "testmid"]
    dir_output = "output/rasters-test/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))"
else
    if fit_sigma
        dir_output = "output/rasters-fit_sigma/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))"
    else
        dir_output = "output/rasters/$(lowercase(Dates.monthabbr(today)))-$(VAWTools.int2str2(Dates.day(today)))"
    end
end
if use_synthetic_dhdt
    dir_output *= "-without-dhdt/"
else
    dir_output *= "/"
end

!isdir(dir_output) && mkpath(dir_output)

region_str = VAWTools.int2str2(region)

@show region, runtyp, parallel, dir_output

include("rgi-setup.jl") # makes sure all is pre-compiles

tmp = BM.get_git_shas(error_on_dirty=false)
description = """
    Now looking into getting Rhat convergence stats after all.

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
glnrs1, glnrs2, glnrs_all = get_glaciers(region)

# To reduce to just a few glacier both of mode1 & 2
if !all_glaciers
    glnrs1 = length(glnrs1)>1 ? glnrs1[1:2] : glnrs1
    glnrs2 = glnrs2[1:5]
end

## MCMC

# parameters to play with
sigma_of_model = @NT(sigma_h_model = 10.0,
                     sigma_iv_model = 10.0)


# number of divisions in elevation band variables
fit_vars = @NT(n_1d_btilde=3,
               n_1d_fsl=3,
               n_1d_temp=1)

fit_target = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.length][1]

@show region, debugging, runtyp, length(glnrs1)

# Full fitting Mode 1
#####################
start_tic = time()

# Method errors contain their arguments, which in our case can be huge.
# Remove them in case of error in pmap
@everywhere strip_MethodErrors(err) = err isa MethodError ? MethodError(0, err.f, err.world) : err
for rm1 in 1:repeat_mode_1
    tic = time()

    # runtyp=:prodlow
    # nr = 2819
    # out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(1), dir_output="", retsol=true, fit_sigma=fit_sigma)
    # error()

    # out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(1), dir_output="", retsol=true, fit_sigma=fit_sigma)
    global out1 = pmap(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
                                    fit_target, runtyp, Val(1), dir_output=dir_output,
                                    fit_sigma=fit_sigma, store_thetas=true, use_synthetic_dhdt=use_synthetic_dhdt),
                       glnrs1, on_error=strip_MethodErrors, retry_delays = zeros(3))

    toc1 = time() - tic
    println("Iteration over mode-1 glaciers completed in $(toc1/60)min")

    ext = repeat_mode_1>1 ? "-$rm1" : ""
    println("Storing $dir_output/$(region_str)_out1$ext.bson ...")
    write(dir_output *"/$(region_str)_out$ext.txt", description)
    @time bson(dir_output *"/$(region_str)_out1$ext.bson", Dict(:out1 => deepcopy(out1),
                                                                :description => description,
                                                                :commit => BM.get_git_shas(error_on_dirty=false),
                                                                :glnrs1 => glnrs1,
                                                                :glnrs2 => glnrs2,
                                                                :runtyp => runtyp,
                                                                ))
    #out1 = BSON.load(dir_output *"/$(region_str)_out1.bson")[:out1];
    println("... done.")

    if length(out1)>0 && length([o[:thetas_expect] for o in out1 if o isa Dict])==0
        error("All Mode-1 runs errored")
    end
end

# Make the priors of the fitted parameters
##########################################
#
# I'm just doing the easy thing here and using the expectation and std of
# the thetas as component-wise prior uppdate.

function calculate_prior_update_reg14()
    # specify files
    fls = ["output/rasters-fit_sigma/production-with-dhdt/07_out1.bson",
           "output/rasters-fit_sigma/production-with-dhdt/11_out1.bson",
           "output/rasters-fit_sigma/production-without-dhdt/03_out1.bson",
           "output/rasters-fit_sigma/production-without-dhdt/04_out1.bson"]

    ## use a mean
    # exs = []
    # sts = []
    # for f in fls
    #     out1 = BSON.load(f)[:out1];
    #     expect, expect_median, expect_mean, stdev = calculate_prior_update(out1, true)
    #     push!(exs, expect_median)
    #     push!(sts, stdev)
    # end
    # return exs, sts, mean(hcat(exs...),2)[:], mean(hcat(sts...),2)[:]

    # merge all runs into one.  I.e. treat it as one big region.
    out1 = []
    for f in fls
        o = BSON.load(f)[:out1]
        # @show calculate_prior_update(o, true)[1]
        append!(out1, o)
    end
    expect, expect_median, expect_mean, stdev = calculate_prior_update(out1, true)
    return expect, stdev

end
if region==14
    expect, stdev = calculate_prior_update_reg14()
    expect_mean = nothing
    expect_median = expect
else
    expect, expect_median, expect_mean, stdev = calculate_prior_update(out1, true)
end

# Priors-only fitting Mode 2
############################
for rm1 in 1:repeat_mode_2
    tic2 = time()
    #out,sol = fit_it(nr, region, sigma_of_model, fit_vars, fit_target, runtyp, Val(2), expect, stdev, dir_output="", retsol=true, fit_sigma=false)
    global out2 = pmap(nr -> fit_it(nr, region, sigma_of_model, fit_vars,
                             fit_target, runtyp, Val(2),
                             expect, stdev, dir_output=dir_output,
                             fit_sigma=false, store_thetas=true, use_synthetic_dhdt=use_synthetic_dhdt),
                glnrs2, on_error=strip_MethodErrors, retry_delays = zeros(3))

    toc2 = time() - tic2

    println("Iteration over mode-2 glacier completed in $(toc2/60)min")
    ext = repeat_mode_2>1 ? "-$rm1" : ""
    println("Storing $dir_output/$(region_str)_out2$ext.bson ...")
    @time bson(dir_output *"/$(region_str)_out2$ext.bson", Dict(:out2 => deepcopy(out2),
                                                                :description => description,
                                                                :commit => BM.get_git_shas(error_on_dirty=false),
                                                                :expect => expect,
                                                                :expect_mean => expect_mean,
                                                                :expect_median => expect_median,
                                                                :stdev => stdev,
                                                                ))
    println("... done.")
end

toc = time() - start_tic
println("Total time $(toc/60/60) hours")

if length(out2)>0 && length([o[:thetas_expect] for o in out2 if o isa Dict])==0
    error("All Mode-2 runs errored")
end
