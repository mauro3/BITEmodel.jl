# invert the ITMIX 2 with IV for ESA

include("itmix-setup.jl") # run this first on single process to make sure all precompilation this through
                          # before the parallel run starts

if nprocs()<2
    if Sys.CPU_CORES>=17
        addprocs(17) # /2 to get to physical cores
    else
        addprocs(Sys.CPU_CORES÷2 + 1) # /2 to get to physical cores
    end
end

@everywhere begin
    const debugging = true
    include("itmix-setup.jl")
end

# whether doing a testing, semi-testing or production run
runtyp = ["test", "testmid", "prodlow", "prod"][3]
runtyp == :test && println("\n\nTEST RUN !!!!!!!!\n\n")

fit_target = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.length][2]  # making this a global to avoid below odd bug.
pl_kws = Dict{Symbol,Any}()

dir = "results_Werder/esa-14nov"
skip_iv_point_measurements = false

for gid in BM.ITMIXGlacier.(BM.itmix_glaciers_iv, 2)
    if skip_iv_point_measurements && gid in BM.ITMIXGlacier.(BM.itmix_glaciers_iv_point, 2)
        continue
    end
    println("\n\n============================  Running Glacier $(BM.getname(gid))")
    for rm_iv in [true, false]
        println("\n\n============================  rm_iv = $(rm_iv)")
        println("\n\n============================")

        gl,gb,pp,pm,pn,pl = BM.init_forward(gid; pl_kws...);

        if nprocs()>1
            pmap(expnr->run_one_exp(expnr, gl, gb, pp, pm, pn, pl, runtyp, fit_target; rm_iv=rm_iv, dir=dir), 0:16)
        else
            map(expnr->run_one_exp(expnr, gl, gb, pp, pm, pn, pl, runtyp, fit_target; rm_iv=rm_iv, dir=dir), 0:16)
        endI==1 ? true : false
    end
end
