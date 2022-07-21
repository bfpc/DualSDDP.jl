import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load

include("postproc_funs.jl")

# Small multiples on alpha/beta
function smallmultiples(bounds_dict; log::Bool=false)
  alphas, betas, lips = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  interact = plt.isinteractive()
  plt.ioff()
  fig, axs = plt.subplots(ncols=nalphas, nrows=nbetas, figsize=(14,14); sharey=true, squeeze=false)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      lb = best_lb(bounds_dict, alpha, beta, lips)
      for lip in lips
        axs[i,j].plot(bounds_dict[(alpha, beta, lip)][2], label=string(lip))
      end
      axs[i,j].axhline(lb, linestyle="--", linewidth=1, color="b")
      axs[i,j].set_title("α = $alpha, β = $beta")
      if log
        axs[i,j].set_yscale("log")
      end
    end
  end
  axs[1,end].legend(title="Lipschitz parameter")
  fig.suptitle("Bounds for varying Lipschitz constants")
  fig.tight_layout()
  if interact
    plt.ion()
  end
  plt.savefig(joinpath("data", "output", cfg["save_path"], "Lip-smallmultiples.pdf"))
  plt.show()
end

#
# Script
#
cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

# Base parameters
parse!(cfg, 1)
dir = cfg["experiment"]["dir"]


# Table and graphs
bounds_dict = read_ablip_all(cfg)
table_relgaps(bounds_dict, "Lip", [1,5,10,20,100])
smallmultiples(bounds_dict; log=true)
