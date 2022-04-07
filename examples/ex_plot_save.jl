import PyPlot as plt
import NPZ

###################
# Graphs and bounds

function plot_step_ub(ub; kwargs...)
  n = length(ub)
  xs = repeat(first.(ub), 1, 2)' |> (x -> reshape(x, 1,2n))
  ys = repeat(last.(ub), 1, 2)' |> (x -> reshape(x, 1,2n))
  plt.plot(xs[2:end], ys[1:end-1]; kwargs...)
end

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

########
# Saving

NPZ.npzwrite("bounds.npz", lb=primal_lbs, ub=dual_ubs, io_lbs=io_lbs, io_ubs=io_ubs, rec_ubs=last.(rec_ubs), rec_iters=first.(rec_ubs));
