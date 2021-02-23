import Pkg
Pkg.activate(".")

include("problem.jl")
include("risk_models.jl")

include("hydro.jl")

f3      = mk_primal_avar(0.3)
f3_dual = mk_copersp_avar(0.3)

primal_pb = mk_primal_decomp(Hydro1d.M, 5, f3)
dual_pb   = mk_dual_decomp(Hydro1d.M, 5, f3_dual)


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

println("Forward on the primal")
forward(primal_pb, [53.222])
println()

println("Forward on the dual")
forward_dual(dual_pb, [0.0])

