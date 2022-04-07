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

# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs = primalsolve(Evil10d.M, nstages, risk, solver, inivol, niters; verbose=true)

# Pure dual
seed!(1)
dual_pb, dual_ubs = dualsolve(Evil10d.M, nstages, risk_dual, solver, inivol, niters; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs = primalub(Evil10d.M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs = problem_child_solve(Evil10d.M, nstages, risk, solver, inivol, niters; verbose=true)

###################
# Graphs and bounds

function plot_step_ub(ub; kwargs...)
  n = length(ub)
  xs = repeat(first.(ub), 1, 2)' |> (x -> reshape(x, 1,2n))
  ys = repeat(last.(ub), 1, 2)' |> (x -> reshape(x, 1,2n))
  plt.plot(xs[2:end], ys[1:end-1]; kwargs...)
end

import PyPlot as plt
plt.figure(figsize=(6,4));
plt.plot(1:niters, dual_ubs, label="Dual Upper bounds");
plot_step_ub(rec_ubs, label="Inner Upper bounds");
plt.plot(1:niters, io_ubs, label="IO Upper bounds");
plt.plot(1:niters, io_lbs, label="IO Lower bounds");
plt.plot(1:niters, primal_lbs, label="Lower bounds");
plt.xlabel("iteration #");
plt.legend();
plt.savefig("bounds.pdf");

plt.yscale("log");
plt.savefig("bounds_semilog.pdf");

gap_ub = [(n, u/primal_lbs[n] - 1) for (n,u) in rec_ubs];
plt.figure(figsize=(6,4));
plt.semilogy(1:niters, dual_ubs./primal_lbs .- 1; label="Dual UB / Primal LB");
plot_step_ub(gap_ub; label="Recursive Inner UB / Primal LB");
plt.semilogy(1:niters, io_ubs./io_lbs .- 1; label="IO UB/LB");
plt.legend();
plt.xlabel("iteration #");
plt.title("Relative gap");
plt.savefig("gap.pdf");

import NPZ
NPZ.npzwrite("bounds.npz", lb=primal_lbs, ub=dual_ubs, io_lbs=io_lbs, io_ubs=io_ubs, rec_ubs=last.(rec_ubs), rec_iters=first.(rec_ubs));
