import Pkg
Pkg.activate("../")

using DualSDDP

lip_factor = 100

include("evil_conf.jl")
include("evil.jl")


risk      = mk_primal_avar(beta)
risk_dual = mk_copersp_avar(beta)


# import Gurobi
# env = Gurobi.Env()
# solver = () -> Gurobi.Optimizer(env)
import GLPK
solver = GLPK.Optimizer

# Solution algorithms
problem_child_solve(Evil10d.M, nstages, risk, solver, inivol, niters; verbose=true)


# Pure dual
using Random: seed!
seed!(1)
dual_pb, dual_ubs = dualsolve(Evil10d.M, nstages, risk_dual, solver, inivol, niters; verbose=true)

# Primal with interior bounds
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_aux, Ubs = primalsolve(Evil10d.M, nstages, risk, solver, inivol, niters; verbose=true, ub=true)

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
plt.semilogy(dual_ubs./primal_lbs .- 1);
plt.xlabel("iteration #");
plt.title("Relative gap");
plt.savefig("gap.pdf");

import NPZ
NPZ.npzwrite("bounds.npz", lb=primal_lbs, ub=dual_ubs);
