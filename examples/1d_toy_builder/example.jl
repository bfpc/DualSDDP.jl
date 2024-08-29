import Pkg
Pkg.activate("../")

using Random: seed!
using DualSDDP
using ..Build

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

M = Hydro1d.M
state0 = Hydro1d.inivol

include("../ex_calc_problems.jl")

include("../ex_plot_save.jl")
