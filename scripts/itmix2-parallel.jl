# invert all the ITMIX 2

include("itmix-setup.jl") # run this first to make sure all precompilation this through
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
runtyp = ["test", "testmid", "prodlow", "prod"][1]
runtyp == :test && println("\n\nTEST RUN !!!!!!!!\n\n")

fit_target = [BM.FitTarget.h, BM.FitTarget.h_iv, BM.FitTarget.length][2]  # making this a global to avoid below odd bug.
pl_kws = Dict{Symbol,Any}()

dir = "results_Werder/tmp"

for priority=1:4
    println("\n\n++++++++++++++++ PRIORITY $priority ++++++++++++++++++++++++++\n\n")
    for gid in I2.priorities[priority] # Note: parallel execution will block here until all processes have finished
                                       # the pmap.
        println("\n\n============================  Running Glacier $(BM.getname(gid))")
        #    try
        gl,gb,pp,pm,pn,pl = BM.init_forward(gid; pl_kws...);
        # I get an odd bug in pmap with ERROR: LoadError: On worker 2: UndefVarError: h_iv not defined
        if nprocs()>1
            pmap(exp->run_one_exp(exp, gl, gb, pp, pm, pn, pl, runtyp, fit_target; dir=dir), 0:I2.nexps)
        else
            map(exp->run_one_exp(exp, gl, gb, pp, pm, pn, pl, runtyp, fit_target; dir=dir), 0:I2.nexps)
        end
    end
end
