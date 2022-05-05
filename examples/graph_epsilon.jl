import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load
using Printf: @printf

# Currently unused
# plot_gap=1e-7

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

# These are "unparametrizable" yet
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]

function read_abe(cfg::ConfigManager)
  bounds_dict = Dict()
  N = total_combinations(cfg)

  for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]
    epsilon = cfg["parameters"]["epsilon"]

    bounds_dict[(alpha, beta, epsilon)] = (data["dual ub"])
  end
  return bounds_dict
end

function sets_in_dict(bounds_dict)
  alphas = Set()
  betas = Set()
  epsilons = Set()

  for (α, β, ϵ) in keys(bounds_dict)
    push!(alphas, α)
    push!(betas, β)
    push!(epsilons, ϵ)
  end
  alphas = sort([alphas...])
  betas = sort([betas...])
  epsilons = sort([epsilons...])

  return alphas, betas, epsilons
end

function table_bounds(bounds_dict)
  alphas, betas, epsilons = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      @printf("Experiment α=%6.3f β=%6.3f:\n", alpha, beta)
      for epsilon in epsilons
        @printf("    ε=%.2e yields final bound at %.17f\n", epsilon, bounds_dict[(alpha, beta, epsilon)][end])
      end
    end
  end
end

# Small multiples on alpha/beta
function smallmultiples(bounds_dict; log::Bool=false)
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
      lb = minimum([bounds_dict[(alpha, beta, e)][end] for e in epsilons])
      for epsilon in epsilons
        axs[i,j].plot(bounds_dict[(alpha, beta, epsilon)] .- lb, label=string(epsilon))
      end
      axs[i,j].set_title("α = $alpha, β = $beta")
      if log
        axs[i,j].yscale("log")
      end
    end
  end
  axs[1,end].legend(title="Regularization ϵ")
  fig.suptitle("Upper bounds in dual algorithm - Lipschitz constant = $lip_factor")
  fig.tight_layout()
  plt.show()
  if interact
    plt.ion()
  end
end

bounds_dict = read_abe(cfg)
table_bounds(bounds_dict)
smallmultiples(bounds_dict)
