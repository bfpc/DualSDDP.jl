import Pkg
Pkg.activate(".")

include("problem.jl")
include("risk_models.jl")

include("hydro_conf.jl")
include("hydro.jl")

risk      = mk_primal_avar(beta)
risk_dual = mk_copersp_avar(beta)

primal_pb = mk_primal_decomp(Hydro1d.M, nstages, risk)
dual_pb   = mk_dual_decomp(Hydro1d.M, nstages, risk_dual)


# Demo forward, needs solver
include("algo.jl")

# import Gurobi
# env = Gurobi.Env()
# solver = () -> Gurobi.Optimizer(env)
import GLPK
solver = GLPK.Optimizer

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
for i = 1:10
  forward(primal_pb, [inivol])
  backward(primal_pb)
  println("Iteration $i: LB = ", JuMP.objective_value(primal_pb[1]))
end

println()

println("********")
println("  DUAL  ")
println("********")
println("Forward-backward Iterations")
init_dual(dual_pb, [inivol])
for i = 1:10
  forward_dual(dual_pb)
  backward_dual(dual_pb)
  println("Iteration $i: UB = ", -JuMP.objective_value(dual_pb[1]))
end


