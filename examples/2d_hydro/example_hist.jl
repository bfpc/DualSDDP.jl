import Pkg
Pkg.activate("../")

using Random: seed!
import JuMP
using DualSDDP

###################
# Parameters & Data

lip_factor = 1
nscen = 0

include("hydro_hist.jl")
alpha = 0.4
beta = 0.7
niters = 300
ub_step = 25

nstages = Hydro_Hist.nstages
inivol  = Hydro_Hist.inivol

risk      = mk_primal_avar(alpha; beta=beta)
risk_dual = mk_copersp_avar(alpha; beta=beta)

import Gurobi
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)

M = Hydro_Hist.M
state0 = inivol

include("../ex_calc_problems.jl")

include("../ex_plot_save.jl")
