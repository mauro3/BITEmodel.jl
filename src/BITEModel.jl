module BITEModel
using VAWTools
import Interpolations
const Interp = Interpolations
import Proj4
using Suppressor
using DocStringExtensions
import DataStructures: OrderedDict
import Roots

const SRangeL{F} = StepRangeLen{F,Base.TwicePrecision{F},Base.TwicePrecision{F}}

# a few short-hands
const F = Float64 # to switch between Float32 and Float64

# Base types
abstract type APara end #{Name} # Name is a symbol to encode the used model (both forward and inverse)

###############
# Misc helpers
###############

include("misc.jl")

################
# Pre-processing
################

# Various data structures for model
include("data-structs.jl")

# All related to data loading
# (probably should be its own package at some stage)
include("data-loading.jl")
include("data-loading-misc.jl")
include("loaders/data-loaders.jl")

###############
# Actual model
###############
# (F&H 2009)

## Forward model implementation, a variation of the Farinotti et al. (2009) model.
include("TwoD.jl")
include("OneD.jl")

# Forward model post-proc:
# plotting, volume,

###############
# Inverse model
###############
# Do some other time...
#
# - Bayesian Fitting
# - Huss ad-hoc fitting

include("mcmc/mcmc.jl")
include("data-loading-inverse.jl")

# """
# TODO

# The inverse model, i.e. the model fitting function.  It needs to
# specify fitting parameters and measurements.  All 2D input/output
# needs to be on the same gird (really?).

# Inputs:
# -------

# Model input with no errors:
# - dem: glacier DEM with ice mask

# Model input with errors (random fields):

# - bdot: mass balance rate
# - dhdt: elevation change rate

# Fitting parameters:

# - fq: function with fitting parameter
# - A: function with fitting parameter



# Measurements of forward model output:

# - us_obs: surface flow speed (usually 2D)
# - h_obs: ice thickness measurements (usually along radar tracks)

# Outputs:
# -------

# - h: ice thickness
# - h_err: error estimate
# - us + us_err: surface flow speed
# - fq + fq_err
# - A + A_err

# """
# function invm end

######
# Post-processing tools
#####
include("postproc.jl")

######
# Plotting
#####
include("plotting.jl")


end # module
