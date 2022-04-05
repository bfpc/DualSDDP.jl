import Pkg
Pkg.activate("../")

import JuMP
using DualSDDP

###################
# Parameters & Data

lip_factor = 1

include("hydro_hist.jl")
alpha = 0.4
beta = 0.7
niters = [300, 300, 300]

nstages = Hydro_Hist.nstages
inivol  = Hydro_Hist.inivol

risk      = mk_primal_avar(alpha; beta=beta)
risk_dual = mk_copersp_avar(alpha; beta=beta)

import Gurobi
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)

#####################
# Solution algorithms

#problem child
problem_child_solve(Hydro_Hist.M, nstages, risk, solver, inivol, niters[1]; verbose=true)

# Pure dual
using Random: seed!
seed!(3)
dual_pb, dual_ubs = dualsolve(Hydro_Hist.M, nstages, risk_dual, solver, inivol, niters[2]; verbose=true)

# Primal with interior bounds
seed!(2)
primal_pb, primal_trajs, primal_lbs, Ubs = primalsolve(Hydro_Hist.M, nstages, risk, solver, inivol, niters[3]; verbose=true, ub=true)

###################
# Graphs and bounds

import PyPlot as plt
plt.figure(figsize=(6,4));
plt.plot(dual_ubs, label="Upper bounds");
plt.plot(primal_lbs, label="Lower bounds");
plt.xlabel("iteration #");
plt.legend();
plt.savefig("bounds.pdf");

plt.yscale("log");
plt.savefig("bounds_semilog.pdf");

plt.figure(figsize=(6,4));
plt.semilogy(dual_ubs./primal_lbs .- 1, label="Dual");
plt.axhline(y=Ubs[1,1]/primal_lbs[end] - 1, linestyle="--", color="black", label="inner")
plt.legend()
plt.xlabel("iteration #");
plt.title("Relative gap");
plt.savefig("gap.pdf");

import NPZ
NPZ.npzwrite("bounds.npz", lb=primal_lbs, ub=dual_ubs);
