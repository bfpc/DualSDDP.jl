import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

N = total_combinations(cfg)

# These are "unparametrizable" yet
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]

bounds_dict = Dict()
alphas = Set()
betas = Set()

for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]

    bounds_dict[(alpha,beta)] = (data["lb"], data["ub"], data["inner"])
    push!(alphas, alpha)
    push!(betas, beta)
end

alphas = sort([alphas...])
betas = sort([betas...])

nalphas = length(alphas)
nbetas = length(betas)

# Small multiples 1: Varying β, keeping α fixed
fig, axs = plt.subplots(ncols=nalphas, figsize=(18,4), sharey=true)
for i in 1:nalphas
  a = alphas[i]
  j = 0
  for b in betas
    lb, ub, inner = bounds_dict[(a,b)]
    gap = ub./lb .- 1
    axs[i].semilogy(gap, label=string(b), color="C$j")
    axs[i].semilogy(first.(inner), last.(inner)./lb[first.(inner)] .- 1, linestyle="--", color="C$j")
    j+=1
  end
  axs[i].legend(title="β")
  axs[i].set_title("α = $a")
end

fig.suptitle(dir * " - Lipschitz factor: $lip_factor")
plt.savefig(joinpath("data", "output", cfg["save_path"], "Lip$(lip_factor).pdf"))


# Small multiples 2: Varying α, keeping β fixed
fig, axs = plt.subplots(ncols=nbetas, figsize=(18,4), sharey=true)
for i in 1:nbetas
  b = betas[i]
  j = 0
  for a in alphas
    lb, ub, inner = bounds_dict[(a,b)]
    gap = ub./lb .- 1
    axs[i].semilogy(gap, label=string(a), color="C$j")
    axs[i].semilogy(first.(inner), last.(inner)./lb[first.(inner)] .- 1, linestyle="--", color="C$j")
    j+=1
  end
  axs[i].legend(title="α")
  axs[i].set_title("β = $b")
end

fig.suptitle(dir * " - Lipschitz factor: $lip_factor")
plt.savefig(joinpath("data", "output", cfg["save_path"], "Lip$(lip_factor)_tr.pdf"))
