import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load

ARGS=["config_4d.json"]

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

N = total_combinations(cfg)

# Read data from all config combinations
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]

bounds_dict = Dict()
alphas = Set()
betas = Set()

for idx=1:N
    println(idx)
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]

    bounds_dict[(alpha,beta)] = (data["primal lb"], data["dual ub"], data["inner recursive bound"],data["inner recursive iters"],data["io lb"],data["io ub"])
    push!(alphas, alpha)
    push!(betas, beta)
end

alphas = sort([alphas...])
betas = sort([betas...])

nalphas = length(alphas)
nbetas = length(betas)


function plot_step(ax, times, values;
                   range::Union{Nothing,Tuple{Int,Int}}=nothing, kwargs...)
  n = length(times)
  xs = repeat(times, 1, 2)'  |> (x -> reshape(x, 1,2n))
  ys = repeat(values, 1, 2)' |> (x -> reshape(x, 1,2n))
  if range == nothing
    ax.plot(xs[2:end], ys[1:end-1]; kwargs...)
  else
    idxs = (range[1] .<= xs .<= range[2])
    ax.plot(xs[idxs], ys[1:end-1][idxs[2:end]]; kwargs...)
  end
end

function plot_bounds(n_init, n_end, title)
    plot_it = collect(n_init:n_end)

    #fig, axs = plt.subplots(nrows=nbetas, ncols=nalphas, figsize=(11,11))
    fig, axs = plt.subplots(nrows=nbetas, ncols=1, figsize=(4,12))
    # horizontal
    #fig, axs = plt.subplots(nrows=1, ncols=nalphas, figsize=(11,11))
    for i in 1:nalphas
        a = alphas[i]
        #for j in 1:nbetas
            b = betas[i]
            axs[i,1].set_title("α=$a, β = $b")

            if haskey(bounds_dict,(a,b))
                lb, ub, inner, inner_it, io_lb, io_ub = bounds_dict[(a,b)]
                # '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
                axs[i,1].plot(plot_it,ub[n_init:n_end], label="Dual UB", color="C1")
                plot_step(axs[i,1], inner_it, inner, range=(n_init,n_end), label="Philpott UB", color="C2")
                axs[i,1].plot(plot_it,io_ub[n_init:n_end], label="Baucke UB", linestyle="--", color="C3")
                axs[i,1].plot(plot_it,io_lb[n_init:n_end], label="Baucke LB",  linestyle="--",color="C4")
                axs[i,1].plot(plot_it,lb[n_init:n_end], label="SDDP LB", color="C5")
                ## horizontal
                # axs[1,i].plot(plot_it,ub[n_init:n_end], label="Dual UB", color="C1")
                # plot_step(axs[1,i], inner_it, inner, range=(n_init,n_end), label="Philpott UB", color="C2")
                # axs[1,i].plot(plot_it,io_ub[n_init:n_end], label="Baucke UB", linestyle="--", color="C3")
                # axs[1,i].plot(plot_it,io_lb[n_init:n_end], label="Baucke LB",  linestyle="--",color="C4")
                # axs[1,i].plot(plot_it,lb[n_init:n_end], label="SDDP LB", color="C5")
            end
        #end
        #axs[1,nbetas].legend() #horizontal
        axs[1,1].legend() #vertical
    end

    fig.suptitle("Model: $dir - Lipschitz factor: $lip_factor")
    fig.tight_layout()
    plt.savefig(joinpath("data", "output", cfg["save_path"], "$title.pdf"))
end

plot_bounds(100, 300, "Bounds_end_diag")
