import Pkg
Pkg.activate("../")

using Random: seed!
using DualSDDP

lip_factor = 100

include("evil_conf.jl")
include("evil.jl")


risk      = mk_primal_avar(alpha)
risk_dual = mk_copersp_avar(alpha)


# import Gurobi
# env = Gurobi.Env()
# solver = () -> Gurobi.Optimizer(env)
import GLPK
solver = GLPK.Optimizer

M = Evil1d.M
state0 = [inivol]

include("../ex_calc_problems.jl")

include("../ex_plot_save.jl")
