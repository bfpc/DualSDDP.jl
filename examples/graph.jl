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
betas = Set()
lambdas = Set()

for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    beta   = cfg["risk-aversion"]["beta"]
    lambda = cfg["risk-aversion"]["lambda"]

    bounds_dict[(beta,lambda)] = (data["lb"], data["ub"], data["inner"][1][1])
    push!(betas, beta)
    push!(lambdas, lambda)
end

betas = sort([betas...])
lambdas = sort([lambdas...])

nbetas = length(betas)
nlambdas = length(lambdas)

# Small multiples 1: Varying λ, keeping β fixed
fig, axs = plt.subplots(ncols=nbetas, figsize=(18,4), sharey=true)
for i in 1:nbetas
  b = betas[i]
  j = 0
  for l in lambdas
    lb, ub, inner = bounds_dict[(b,l)]
    gap = ub./lb .- 1
    axs[i].semilogy(gap, label=string(l), color="C$j")
    axs[i].axhline(y=inner/lb[end] - 1, linestyle="--", color="C$j")
    j+=1
  end
  axs[i].legend(title="λ")
  axs[i].set_title("β = $b")
end

fig.suptitle(dir * " - Lipschitz factor: $lip_factor")
plt.savefig(joinpath("data", "output", cfg["save_path"], "Lip$(lip_factor).pdf"))


# Small multiples 2: Varying β, keeping λ fixed
fig, axs = plt.subplots(ncols=nlambdas, figsize=(18,4), sharey=true)
for i in 1:nlambdas
  l = lambdas[i]
  j = 0
  for b in betas
    lb, ub, inner = bounds_dict[(b,l)]
    gap = ub./lb .- 1
    axs[i].semilogy(gap, label=string(b), color="C$j")
    axs[i].axhline(y=inner/lb[end] - 1, linestyle="--", color="C$j")
    j+=1
  end
  axs[i].legend(title="β")
  axs[i].set_title("λ = $l")
end

fig.suptitle(dir * " - Lipschitz factor: $lip_factor")
plt.savefig(joinpath("data", "output", cfg["save_path"], "Lip$(lip_factor)_tr.pdf"))
