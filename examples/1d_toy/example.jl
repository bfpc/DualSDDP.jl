maindir = "../../"
src     = maindir * "src/"

import Pkg
Pkg.activate(maindir)

include(src * "problem.jl")
include(src * "risk_models.jl")
include(src * "algo.jl")
include(src * "ub.jl")

include("hydro_conf.jl")
include("hydro.jl")


risk      = mk_primal_avar(beta)
risk_dual = mk_copersp_avar(beta)


# import Gurobi
# env = Gurobi.Env()
# solver = () -> Gurobi.Optimizer(env)
import GLPK
solver = GLPK.Optimizer

# Solution algorithms
# Pure primal
primalsolve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true)

# Pure dual
using Random: seed!
seed!(1)
dual_pb, dual_ubs = dualsolve(Hydro1d.M, nstages, risk_dual, solver, [inivol], niters; verbose=true)

# Primal with interior bounds
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_aux, Ubs = primalsolve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true, ub=true)
nothing
