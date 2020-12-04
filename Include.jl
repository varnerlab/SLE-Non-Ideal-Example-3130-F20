# activate my project -
import Pkg
Pkg.activate()
Pkg.instantiate()

# load packages -
using Plots
using Optim

# load my codes -
include("./src/Compute.jl")
include("./src/Estimate.jl")