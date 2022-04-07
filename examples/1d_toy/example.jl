import Pkg
Pkg.activate("../")

using Random: seed!
using DualSDDP

lip_factor = 100

include("hydro_conf.jl")
include("hydro.jl")


risk      = mk_primal_avar(alpha)
risk_dual = mk_copersp_avar(alpha)


# import Gurobi
# env = Gurobi.Env()
# solver = () -> Gurobi.Optimizer(env)
import GLPK
solver = GLPK.Optimizer

# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs = primalsolve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true)

# Pure dual
seed!(1)
dual_pb, dual_ubs = dualsolve(Hydro1d.M, nstages, risk_dual, solver, [inivol], niters; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs = primalub(Hydro1d.M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs = problem_child_solve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true)

include("../ex_plot_save.jl")
