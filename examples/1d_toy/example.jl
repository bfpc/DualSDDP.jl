import Pkg
Pkg.activate("../")

using Random: seed!
using DualSDDP

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

# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs = primalsolve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true)

# Pure dual
seed!(1)
dual_pb, dual_ubs = dualsolve(Hydro1d.M, nstages, risk_dual, solver, [inivol], niters; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs = primalub(Hydro1d.M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
io_pb, io_lbs, io_ubs = problem_child_solve(Hydro1d.M, nstages, risk, solver, [inivol], niters; verbose=true)

function plot_step_ub(ub; kwargs...)
  xs = repeat(first.(ub), 1, 2)' |> (x -> reshape(x, 1,8))
  ys = repeat(last.(ub), 1, 2)' |> (x -> reshape(x, 1,8))
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
NPZ.npzwrite("bounds.npz", lb=primal_lbs, ub=dual_ubs);
