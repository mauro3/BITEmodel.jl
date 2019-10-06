"""
Runs the forward model for a given theta:

hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb, bdot1d, dhdt1d = run_theta(theta, gl, gb, pm, pn, pp);
"""
function run_theta(theta, gl, gb, pm, pn, pp)
    pmm = MPara(pm, todict(theta))

    # this was missing:
    dhdt1d = [dhdtfn(b, pmm.dhdt_cutoff, pmm.dhdt_base) for b in gb.bands]
    bdot1d = map_onto_bands(gb, gl.bdot, mean, FILL)

    @time hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb, bdot1d, dhdt1d =
        fwdm(gb, gl, pp, pmm, pn, bdot1d, dhdt1d)
    return hs1d, taus1d, ivs1d, hs2d, taus2d, ivs2d, taus1d_l, vol_ratio, gb, bdot1d, dhdt1d
end

plan = """

- MPara gives the parameters needed.
- gl+pm give an "instance" of a glacier to run the model over.
  - mostly pm would have information on how to change some of the 1D input variables
  - this would need some checking whether gb needs updating

"""

"""
Returns a dict of SHAs and branches.  Errors if one of the repos is dirty,
unless error_on_dirty=false.
"""
function get_git_shas(;error_on_dirty=true)
    LG = Base.LibGit2
    out = Dict{String,Any}()

    for rep in ["KissMCMC", "VAWTools", "BITEModel"]
        gr = LG.GitRepo(Pkg.dir(rep))
        error_on_dirty && LG.isdirty(gr) && error("KissMCMC is dirty.  Either commit or force with error_on_dirty=false")
        out[rep] = [string(LG.head_oid(gr)), string(LG.headname(gr))]
    end
    return out
end


# # patch JLD: Not needed for JLD2
# JLD.readas(::Val{:FunctionSerializer}) = (x...) -> error("Storing functions with JLD is not supported")
# @suppress JLD.writeas(fn::Function) = Val{:FunctionSerializer}()


"""
Saves an MCMC model run.
"""
function savemcmc(thetas, blobs, gb::Bands,
                  pp::Phys, pm::MPara, pn::Num, pmcmc, #::MCMCNum,
                  theta0, run::Symbol, sigmas::Tuple;
                  flname="", force_overwrite=false, error_on_dirty=true)
    gl = gb.gl
    sh, siv1d, siv2d = sigmas
    sh = round(Int,sh)
    siv2d = round(Int,siv2d)
    if flname==""
        region = getregion(gl)
        lab = getname(gl)
        flname = "output/$region-$lab-$run-$sh-$siv2d.jld"
    end
    if ~force_overwrite && isfile(flname)
        error("File $flname exists!  Aborting.  Use 'force_overwrite=true'.")
    end
    shas = get_git_shas(error_on_dirty=error_on_dirty)

    file = JLD.jldopen(flname, "w")
    try
        file["thetas"] = thetas
        file["blobs"] = blobs
        file["gb"] = gb
        file["pp"] = pp
        file["pm"] = pm
        file["pn"] = pn
        file["pmcmc"] = pmcmc
        file["theta0"] = theta0
        file["run"] = run
        file["sigmas"] = sigmas
        file["git_shas"] = shas
        file["timestamp"] = Dates.now()
    finally
        close(file)
    end
    return flname
end

"""
Load results.

better use @load

gl, gb, pp, pm, pn, res, git_shas, timestamp = loadmcmc(flname)

thetas_e, accept_ratio_e, blobs_e = squash_chains(res... , drop_low_accept_ratio=true, blob_reduce! = blob_reduce!)
mc_vol, mc_vol_above, mc_hs1d, mc_ivs1d, mh2d, sh2d, miv2d, siv2d, pms = unpack_blob(blobs_e, size(gl.dem.v))
"""
function loadmcmc(flname)
    @unpack thetas, blobs, gb, pp, pm, pn, pmcmc,
    theta0, run, sigmas, git_shas, timestamp = JLD.load(flname)
    return thetas, blobs, gb, pp, pm, pn, pmcmc,
    theta0, run, sigmas, git_shas, timestamp
end
