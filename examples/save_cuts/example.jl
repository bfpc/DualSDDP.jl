import Pkg
Pkg.activate("../")

using Random: seed!
using DualSDDP

# Parameters
nstages = 4
alpha   = 0.4
niters  = 40
ub_step = 10

nscen = 10

lip_factor = 100

include("hydro.jl")

risk      = mk_primal_avar(alpha)
risk_dual = mk_copersp_avar(alpha)

import GLPK
solver = GLPK.Optimizer

M = Hydro1d.M
state0 = Hydro1d.inivol

include("../../src/save_cuts.jl")

seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, niters; verbose=false)

write_cuts_to_file(primal_pb, "save_cuts/cuts_primal.json")

seed!(2)
dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=false)

write_cuts_to_file(dual_pb, "save_cuts/cuts_dual.json")


