import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load

include("postproc_funs.jl")

function gap_lines(values, ax, label, color)
  lb, ub, inner, inner_it, io_lb, io_ub = values
  gap       = ub./lb .- 1
  gap_inner = inner./lb[inner_it] .- 1
  gap_io    = io_ub./io_lb .- 1
  axs[i].semilogy(gap, label=label, color=color)
  axs[i].semilogy(inner_it, gap_inner, linestyle="--", color=color)
  axs[i].semilogy(gap_io, linestyle="dashdot", color=color)
end

function gaps_ab(bounds_dict, title, filename)
  # Small multiples: Varying β, keeping α fixed
  alphas, betas = sets_in_dict(bounds_dict)

  nalphas = length(alphas)

  fig, axs = plt.subplots(ncols=nalphas, figsize=(18,4), sharey=true)
  for i in 1:nalphas
    a = alphas[i]
    j = 0
    for b in betas
      if haskey(bounds_dict,(a,b))
        gap_lines(bounds_dict[(a,b)], axs[i], string(b), "C$j")
        j+=1
      end
    end
    axs[i].legend(title="β")
    axs[i].set_title("α = $a")
  end

  fig.suptitle(title)
  fig.tight_layout()
  plt.savefig(filename)
end

function gaps_ab_tr(bounds_dict, title, filename)
  # Small multiples: Varying α, keeping β fixed
  alphas, betas = sets_in_dict(bounds_dict)

  nbetas = length(betas)
  fig, axs = plt.subplots(ncols=nbetas, figsize=(18,4), sharey=true)
  for i in 1:nbetas
    b = betas[i]
    j = 0
    for a in alphas
      if haskey(bounds_dict,(a,b))
        gap_lines(bounds_dict[(a,b)], axs[i], string(a), "C$j")
        j+=1
      end
    end
    axs[i].legend(title="α")
    axs[i].set_title("β = $b")
  end

  fig.suptitle(title)
  fig.tight_layout()
  plt.savefig(filename)
end

#
# Script
#
# Upper and lower bounds for several algorithms
cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)
bounds_dict = read_ab_all(cfg)

# Base parameters
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]


title = "Model: $dir - Lipschitz factor: $lip_factor
         solid: primal-dual gap / dashed: primal-recursive gap / dashdot: outer-inner gap"
basename = joinpath("data", "output", cfg["save_path"], "Gaps")

gaps_ab(bounds_dict, title, basename * ".pdf")
gaps_ab_tr(bounds_dict, title, basename * "-tr.pdf")
