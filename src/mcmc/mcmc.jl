# using MCMC to fit a Bayesian model

# inverse
using KissMCMC
using StatsBase
using OnlineStats
import JLD2
const JLD=JLD2

export thick_err, iv2d_err, iv1d_err, make_loglikelihood, make_logposterior,
    MCMCPara, mcmc, Theta0

# misc tools
include("tools.jl")

include("probabilistic-model.jl")
include("priors.jl")

# numerical fitting
include("numerics.jl")

include("plotting-inverse.jl")
include("postproc.jl")
include("evaluation.jl")
