maindir = "../../"
src     = maindir * "src/"

import Pkg
Pkg.activate(maindir)

include(src * "problem.jl")
include(src * "risk_models.jl")
include(src * "algo.jl")
include(src * "ub.jl")

include("hydro_hist.jl")
beta = 0.4
niters = [300, 300, 300]

nstages = Hydro_Hist.nstages
inivol  = Hydro_Hist.inivol

risk      = mk_primal_avar(beta)
risk_dual = mk_copersp_avar(beta)

import Gurobi
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)
# import GLPK
# solver = GLPK.Optimizer
# 

# Solution algorithms
# Pure primal
primal_pb, primal_trajs, primal_lbs = primalsolve(Hydro_Hist.M, nstages, risk, solver, inivol, niters[1]; verbose=true)

# Pure dual
# Currently not working: with more scenarios does not show progress
using Random: seed!
seed!(1)
dual_pb, dual_ubs = dualsolve(Hydro_Hist.M, nstages, risk_dual, solver, inivol, niters[2]; verbose=true)

# Primal with interior bounds
# Currently not working: sometimes the state escapes the convex hull
seed!(2)
_, _, _, primal_aux, Ubs = primalsolve(Hydro_Hist.M, nstages, risk, solver, inivol, niters[3]; verbose=true, ub=true)
nothing
