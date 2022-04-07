import Pkg
Pkg.activate("../")

using Random: seed!
import JuMP
using DualSDDP

###################
# Parameters & Data

lip_factor = 1

include("hydro_hist.jl")
alpha = 0.2
beta = 0.3
niters = 2000
ub_step = 50

nstages = Hydro_Hist.nstages
inivol  = Hydro_Hist.inivol

risk      = mk_primal_avar(alpha; beta=beta)
risk_dual = mk_copersp_avar(alpha; beta=beta)

import Gurobi
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)

#####################
# Solution algorithms

# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs = primalsolve(Hydro_Hist.M, nstages, risk, solver, inivol, niters; verbose=true)

# Pure dual
# Currently not working: with more scenarios does not show progress
seed!(3)
dual_pb, dual_ubs = dualsolve(Hydro_Hist.M, nstages, risk_dual, solver, inivol, niters; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs = primalub(Hydro_Hist.M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs = problem_child_solve(Hydro_Hist.M, nstages, risk, solver, inivol, niters; verbose=true)

include("../ex_plot_save.jl")
