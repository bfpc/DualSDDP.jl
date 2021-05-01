maindir = "../../"
src     = maindir * "src/"

import Pkg
Pkg.activate(maindir)

include(src * "problem.jl")
include(src * "risk_models.jl")
include(src * "algo.jl")
include(src * "ub.jl")

include("hydro_ar1.jl")
beta = 0.4
niters = 40

nstages = HydroAR1.nstages
inivol  = HydroAR1.inivol
iniflow = HydroAR1.u0

risk      = mk_primal_avar(beta)
risk_dual = mk_copersp_avar(beta)

primal_pb = mk_primal_decomp(HydroAR1.M, nstages, risk)
dual_pb   = mk_dual_decomp(HydroAR1.M, nstages, risk_dual)


import Gurobi
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)
# import GLPK
# solver = GLPK.Optimizer
# 
for m in primal_pb
  JuMP.set_optimizer(m, solver)
end
for m in dual_pb
  JuMP.set_optimizer(m, solver)
end

println("********")
println(" PRIMAL ")
println("********")
println("Forward-backward Iterations")
for i = 1:niters
  forward(primal_pb, [inivol; iniflow])
  backward(primal_pb)
  println("Iteration $i: LB = ", JuMP.objective_value(primal_pb[1]))
end
println()

# Currently not working, probably needs "Lipschitz" bounds on dual variables
println("********")
println("  DUAL  ")
println("********")
println("Forward-backward Iterations")
init_dual(dual_pb, [inivol; iniflow])
for i = 1:niters
  forward_dual(dual_pb)
  backward_dual(dual_pb)
  println("Iteration $i: UB = ", -JuMP.objective_value(dual_pb[1]))
end

# Currently not working, probably needs "Lipschitz" penalties
primal_pb = mk_primal_decomp(HydroAR1.M, nstages, risk);
stages    = mk_primal_decomp(HydroAR1.M, nstages, risk);
for m in primal_pb
  JuMP.set_optimizer(m, solver)
end
for m in stages
  JuMP.set_optimizer(m, solver)
end
traj = forward_backward(primal_pb, niters, [inivol; iniflow]; return_traj=true);
Ubs = convex_ub(stages,traj)
println("Recursive upper bounds on $niters trajectories")
Ubs

