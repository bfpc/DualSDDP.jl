import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load
using Printf: @printf

include("postproc_funs.jl")

function table_bounds(bounds_dict, primal_dict)
  alphas, betas, epsilons = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      @printf("Experiment α=%6.3f β=%6.3f:\n", alpha, beta)
      if (alpha, beta) in keys(primal_dict)
          @printf("    primal lb=%.17f\n", primal_dict[(alpha,beta)][1][end])
          @printf("    io     lb=%.17f\n", primal_dict[(alpha,beta)][5][end])
      end
      for epsilon in epsilons
        @printf("    ε=%.2e yields final upper bound at %.17f\n", epsilon, bounds_dict[(alpha, beta, epsilon)][end])
      end
    end
  end
end

# Small multiples on alpha/beta
function smallmultiples(bounds_dict, primal_dict; log::Bool=false)
  alphas, betas, epsilons = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  interact = plt.isinteractive()
  plt.ioff()
  fig, axs = plt.subplots(ncols=nalphas, nrows=nbetas, figsize=(14,14))
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      if (alpha, beta) in keys(primal_dict)
          pd = primal_dict[(alpha, beta)]
          lb = max(pd[1][end], pd[5][end])
      else
          println("  Warning: $((alpha,beta)): using lower bound from dual algorithm only!")
          lb = minimum([bounds_dict[(alpha, beta, e)][end] for e in epsilons])
      end
      for epsilon in epsilons
        axs[i,j].plot(bounds_dict[(alpha, beta, epsilon)] .- lb, label=string(epsilon))
      end
      axs[i,j].set_title("α = $alpha, β = $beta")
      if log
        axs[i,j].set_yscale("log")
      end
    end
  end
  axs[1,end].legend(title="Regularization ϵ")
  fig.suptitle("Upper bounds in dual algorithm - Lipschitz constant = $lip_factor")
  fig.tight_layout()

  figname = joinpath("data", "output", cfg["save_path"], "Epsilon-smallmultiples")
  if log
    figname *= "-log.pdf"
  else
    figname *= ".pdf"
  end
  plt.savefig(figname)

  if interact
    plt.ion()
  end
  plt.show()
end

#
# Script
#
# Dual bounds with varying epsilon
cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)
dual_dict = read_abe_dualub(cfg)

# Base parameters
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]

# Primal baselines
cfg_name = ARGS[2]
cfg = ConfigManager(cfg_name, @__DIR__)

primal_dict = read_ab_all(cfg)


table_bounds(dual_dict, primal_dict)
smallmultiples(bounds_dict, primal_dict; log=true)
