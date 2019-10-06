# Evaluates the ITMIX runs by creating a `data` DataFrame containing summary info and, optionally, plots.
#
# NOTE: the run-plots show the mean thickness field not the one produced by using mtheta in the forward model.

error("remove Rhat")

# run this first to make sure all precompilation this through
include("itmix-setup.jl")
using FileIO, DataFrames, StatsBase
import Plots. PyPlot, Interpolations
global pl

const make_plots = false
const make_plots_thetas = make_plots
parallel = false

n_1d_btilde = 3 # number of spatial parameters
n_1d_fsl = 3 # number of spatial parameters
n_1d_temp = 1 # number of spatial parameters
subdir = ["",
          "res-7oct-2",
          "res-23-oct",
          "res-28-oct",
          "res-30-oct",
          "res-31-oct-rmiv",
          "res-03-nov-rmiv",
          "esa-14nov"
          ][end]

if !isdefined(:runit)
    runit = true
end


"""Construct an empty DataFrame"""
function make_data()
    error("remove Rhat")
    DataFrame(person = Symbol[],
              gl = Symbol[],
              prio = Int[],
              expnr = Int[],
              rm_iv = Bool[],
              nl = Int[],
              mean_h = Float64[],
              volume_area = Float64[],
              rmse = Float64[],
              rmse_fitting = Float64[],
              max = Float64[],
              max_fitting = Float64[],
              Rhat = Float64[],
              Rhat_max = Float64[],
              mtheta = Vector{Float64}[],
              err_h = VAWTools.Traj{Float64}[]
              )
end

function load_one(gid, expnr, rm_iv;
                  fit_target=BM.FitTarget.h_iv)
    gl,gb,pp,pm,pn,pl = gls[gid]
    priority = find([gid in pr for pr in I2.priorities])[1]
    fl, fls, fliv, flivs, flj = I2.filepaths(gl, pl, expnr, rm_iv, "results_Werder/"*subdir)
    dir_,fl_ = splitdir(fl)
    fl_ = splitext(fl_)[1]
    if !isfile(fl)
        println("File $fl not found")
        return nothing
    end

    jldloaded, jldd = try
        true, FileIO.load(flj)
    catch e
        println("File $flj not loading with $e")
        false, nothing
    end

    if jldloaded
        thetas, Rhat, sample_size = map(k->jldd[k], ("thetas","Rhat","sample_size"))
        mtheta = squeeze(mean(thetas,2),2)

        gln = BM.ITMIX2.make_ITMIX2_glacier(gl, expnr)
        gbn = BM.Bands(gb, gl=gln)
        th0d = BM.theta0_dict_defaults(pm, n_1d_btilde, n_1d_fsl, n_1d_temp)
        if length(mtheta)==8
            # make it compatible with older runs
            pop!(th0d, :iv_dist_exp2)
        end
        (theta0, logposterior, logposterior1d, logprior, logprior_debug,
         loglikelihood, loglikelihood1d, fwdm_fn, fwdm1d_fn,
         pmcmc_defaults, fit_target) =
             BM.init_inverse(gbn, pp, pm, pn,
                             theta0_dict=th0d,
                             fit_target=fit_target)
        fwdsol = fwdm_fn(theta0.th0)
    end

    r_lines_fitting = I2.get_rlines(gid, expnr)
    hs2d = VAWTools.read_agr(fl)
    hs2ds = VAWTools.read_agr(fls)
    iv2d = try VAWTools.read_agr(fliv) end
    iv2ds = try VAWTools.read_agr(flivs) end
    if gid==BITEModel.ITMIXGlacier(:Columbia, 2) || gid==BITEModel.ITMIXGlacier(:Elbrus, 2)
        # downsample agr
        hs2d = downsample(hs2d, 2)
        hs2ds = downsample(hs2ds, 2)
    end

    # tuples (control, fitting, all)
    RMSE, max_err, quantiles, smaller_20, smaller_std, smaller_2std, err_h = BM.compare_to_radar(gl, hs2d, hs2ds, r_lines_fitting)

    mean_ice_thickness, volume_area = VAWTools.meannan(hs2d.v), BM.volume_area_mean_h(gl)

    data_out = if jldloaded
        [:Werder, gid.id, priority, expnr, rm_iv, length(r_lines_fitting),
         map(x->signif(x,3), (mean_ice_thickness, volume_area, RMSE[1], RMSE[2], (max_err[1]), (max_err[2]), (mean(Rhat)), (maximum(Rhat))))...,
         mtheta, err_h]
    else
        [:Werder, gid.id, priority, expnr, rm_iv, length(r_lines_fitting),
         map(x->signif(x,3), (mean_ice_thickness, volume_area, RMSE[1], RMSE[2], (max_err[1]), (max_err[2]), 0.0, 0.0))...,
         Float64[], err_h]
    end

    fwdsol = fwdm_fn(mtheta);

    return data_out, hs2d, hs2ds, iv2d, iv2ds, r_lines_fitting, dir_, fl_, fwdsol, jldd, mtheta, thetas, gbn
end

# """
# Evaluate one ITMIX2 experiment and store/save relevant input.

# Returns the data to be appended to data-DataFrame
# """
@everywhere function evaluate_one(expnr, rm_iv, gb, pp, pm, pn, pl;
                                  n_1d_btilde=3, n_1d_fsl=3, n_1d_temp=1,
                                  make_plots = false,
                                  make_plots_thetas = make_plots,
                                  subdir="",
                                  fit_target=BM.FitTarget.h_iv )
    gl = gb.gl
    gid = gl.gid
    (data_out, hs2d, hs2ds, iv2d, iv2ds, r_lines_fitting, dir_,
     fl_, fwdsol, jldd, mtheta, thetas, gbn) = load_one(gid, expnr, rm_iv, fit_target=fit_target)

    # plot it
    if make_plots
        BM.plot2d_h(gl, hs2d.v; clims=:radar, plotline=I2.get_rlines(gid, expnr), reuse=true,
                    outlines=true)
        Plots.savefig(joinpath(dir_, "plots", "$(fl_)-hs2d.png"))
        BM.plot2d_h(gl, hs2ds.v; clims=:model, plotline=I2.get_rlines(gid, expnr), plottraj=false,
                    colorbar_title="h std (m)", reuse=true, outlines=true)
        Plots.savefig(joinpath(dir_, "plots", "$(fl_)-hs2d-std.png"))
        if iv2d!=nothing
            BM.plot2d_h(gl, iv2d.v; plotiv=true, reuse=true, clims=gl.iv.vals isa VAWTools.Traj ? :radar : :model,
                        outlines=true)
            Plots.savefig(joinpath(dir_, "plots", "$(fl_)-iv2d.png"))
            if !(gl.iv isa BM.SomeData{Void})
                if gl.iv.vals isa VAWTools.Gridded
                    err = iv2d.v-gl.iv.vals.v
                    err[.!gl.ivmask] = NaN
                    mi, ma = VAWTools.extremanan(err)
                    clim = max(abs(mi), abs(ma)) # set clims symmetric because of colorscale
                    BM.plot2d_h(gl, err; plotiv=true, colorbar_title="iv model-data (m)",
                                seriescolor=:pu_or, clims=(-clim,clim), reuse=true, outlines=true)
                    Plots.savefig(joinpath(dir_, "plots", "$(fl_)-iv2d-err.png"))
                elseif gl.iv.vals isa VAWTools.Traj
                    #BM.plot2d_h(gl, iv2d.v; clims=:radar, plotiv=true, colorbar_title="iv model-data (m)")
                else
                    error()
                end
            end
            BM.plot2d_h(gl, iv2ds.v; plotiv=true, colorbar_title="iv std (m)", reuse=true,
                        plottraj=false, clims=:model, outlines=true)
            Plots.savefig(joinpath(dir_, "plots", "$(fl_)-iv2d-std.png"))
        end
        # this isn't quite correct as, e.g., hs1d of mtheta is not equal to expected hs1d
        #
        # NOTE this plots the mean hs2d, but the others are from running the model with mean paras
        BM.plot_run(fwdsol,
                    reuse=true, plot2d=[:h,:iv,:ivlog][1], clims=[:model,:radar][2],
                    plotline=I2.get_rlines(gid, expnr),
                    radarline_gl=gl)
        Plots.savefig(joinpath(dir_, "plots", "$(fl_)-run.png"))
    end
    if make_plots_thetas
        names = jldd["theta0_names"]
        theta0 = BM.Theta0{Float64}(mtheta, names,
                                    OrderedDict{Symbol,LinSpace{Float64}}(), gbn, pp, pm, pn)
        BM.plottheta(thetas, theta0, fsize=(1000,1000), reuse=true)
        Plots.savefig(joinpath(dir_, "plots", "$(fl_)-thetas.png"))
    end
    # PyPlot.close("all") does not help
    return data_out
end


# Good queries
# data[data[:gl].==:Unteraar,:]
# by(data, :gl, d->DataFrame( mean_rmse = mean(d[:rmse]), max_err = maximum(d[:max]), runs=length(d[:expnr])))

function summarize_data(data; expnrs=0:16, person=false)
    d = data[ in.(data[:expnr], (expnrs,)), :]
    bycols = person ? [:person, :gl] : [:gl]
    by(d, bycols, d->DataFrame(prio = round(Int, mean(d[:prio])),
                               median_mean_h =  median(d[:mean_h]),
                               quantiles_mean_h =  (summarystats(d[:mean_h]).q25, summarystats(d[:mean_h]).q75),
                               extrema_mean_h =  extrema(d[:mean_h]),
                               mean_rmse = mean(d[:rmse]),
                               rel_mean_rmse = mean(d[:rmse])/median(d[:mean_h]),
                               mean_max_err = mean(d[:max]),
                               maximum_max_err = maximum(d[:max]),
                               runs=length(d[:expnr]),
                               Rhat_median=median(d[:Rhat_max]),
                               Rhat_max=maximum(d[:Rhat_max])
                               ),
       sort = true)

    # Easier:
    # flds = [:gl, :rmse, :rmse_fitting, :max, :max_fitting, :Rhat, :Rhat_max]
    # sort(aggregate(d[flds], :gl, [mean, median, maximum]))
end

@everywhere begin
    include("itmix-setup.jl")
    using FileIO, DataFrames
    import Plots. PyPlot, Interpolations
    pl_kws = Dict{Symbol,Any}()
    PyPlot.matplotlib[:interactive](false) # memory leak https://github.com/JuliaPy/PyPlot.jl/issues/111

# global priority,gid,gl,gb,pp,pm,pn,pl,fl,fls,flj,hs2d,hs2ds,RMSE,max_err,smaller_20,thetas,theta0
# global fwdm_fn, jldd, gbn, mtheta, names


    gids = [BM.ITMIXv2_glaciers, BM.ITMIXGlacier.(BM.itmix_glaciers_iv,2)][2]
    if !isdefined(:gls) || length(gls)<length(gids)
        gls = try
            open("data/ITMIX2/results_Werder/$subdir/gls.jls", "r") do io
                deserialize(io)
            end
        catch
            gls = Dict()
            for gid in gids
                gls[gid] = BM.init_forward(gid; pl_kws...);
            end
            gls
        end
    end
end
# # load and save gls with
# open("data/ITMIX2/results_Werder/$subdir/gls.jls", "w") do io
#     serialize(io, gls)
# end
# gls = open("data/ITMIX2/results_Werder/$subdir/gls.jls", "r") do io
#     deserialize(io)
# end


if runit
    println("""
    ==================
    Running over directory $(joinpath(first(gls)[2][end].data_root,"results_Werder", subdir))
    ==================
    """)

    data = make_data()

    #gid = BITEModel.ITMIXGlacier(:Austfonna, 2)
    for priority=1:4, rm_iv in [false,true]
        println("\n\n++++++++++++++++ PRIORITY $priority ; rm_iv = $rm_iv ++++++++++++++++++++++++++\n\n")
        for gid in I2.priorities[priority] #[BM.ITMIXGlacier(:SouthGlacier,2)] #
            if !haskey(gls, gid)
                continue
            end
            gl,gb,pp,pm,pn,pl = gls[gid]

            dir_ = splitdir(I2.filepaths(gl, pl, 0, rm_iv, "results_Werder/"*subdir)[1])[1]
            !isdir(joinpath(dir_, "plots")) && mkdir(joinpath(dir_, "plots"))

            println("\n\n============================  Glacier $(BM.getname(gid))")
            #println("expnr, nl, rmse, rmse_fitting, max, Rhat, Rhat_max")
            if parallel
                out = pmap(expnr -> evaluate_one(expnr, rm_iv, gb, pp, pm, pn, pl;
                                                 make_plots=make_plots,
                                                 make_plots_thetas=make_plots,
                                                 subdir=subdir),
                           0:I2.nexps)
                [push!(data, o) for o in out if o!=nothing]
            else
                for expnr=0:I2.nexps

                println("\n\n                              Exp $expnr")
                    #try
                    data_out = evaluate_one(expnr, rm_iv, gb, pp, pm, pn, pl;
                                            make_plots=make_plots,
                                            make_plots_thetas=make_plots,
                                            subdir=subdir)
                    data_out!=nothing && push!(data, data_out)
                    # catch e
                    #     e isa InterruptException && rethrow(e)
                    #     print(e)
                # end
                end
            end
        end
    end
    # save data:
    FileIO.save(joinpath(pl.data_root,"results_Werder", subdir, "data.jld"), "data", data,
                "data_summarized", summarize_data(data)) #, "gls", gls) saving gls takes 3.5GiB!!!

    ## Load with
    # using DataFrames, FileIO, BITEModel
    # data, data_summarized = FileIO.load("data/ITMIX2/results_Werder/res-30-oct/data.jld", "data", "data_summarized")

end

# read others results

function read_others_results(gid, expnr, dir, person; prefix="thickness", thin=1, start=1)
    gl,gb,pp,pm,pn,pl = gls[gid]
    r_lines_fitting = I2.get_rlines(gid, expnr)

    expnr_ = expnr>9 ? Symbol("exp$expnr") : Symbol("exp0$expnr")
    fl = joinpath(dir, "$(prefix)_$(BM.getname(gid))_$(person)_$(expnr_).asc")
    hs2d = VAWTools.read_agr(fl)
    hs2d.v[hs2d.v.==0] = NaN
    if thin>1
        hs2d = VAWTools.downsample(hs2d, thin, start)
    end

    RMSE, max_err, quantiles, smaller_20, smaller_std, smaller_2std, err_h = BM.compare_to_radar(gl, hs2d, nothing, r_lines_fitting)
    mean_ice_thickness, volume_area = VAWTools.meannan(hs2d.v), BM.volume_area_mean_h(gl)

    data_out = [person, gid.id, -1, expnr, true, length(r_lines_fitting),
                map(x->signif(x,3), (mean_ice_thickness, volume_area, RMSE[1],
                                     RMSE[2], (max_err[1]), (max_err[2]), NaN, NaN))..., Float64[], err_h]
    return data_out
end
