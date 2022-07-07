import Pkg
Pkg.activate(".")

ARGS=["config_4d_scen.json"]

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load



cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

N = total_combinations(cfg)

# Read data from all config combinations
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]


times_dict = Dict()
bounds_dict = Dict()
alphas = Set()
betas = Set()
branches = Set()

for idx=1:N
    println(idx)
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]

    times_dict[(alpha,beta)] = (data["primal t"], data["dual t"], data["inner recursive t"], data["io t"])
    bounds_dict[(alpha,beta)] = (data["primal lb"], data["dual ub"], data["inner recursive bound"],data["inner recursive iters"],data["io lb"],data["io ub"])

    push!(alphas, alpha)
    push!(betas, beta)
end

alphas = sort([alphas...])
betas = sort([betas...])

nalphas = length(alphas)
nbetas = length(betas)

function plot_bound_times(n_init, n_end, title)
# TODO
end

function moving_average(data;window=20)
    n = length(data)
    ma = zeros(n)
    for i in 1:window
        ma[i]= sum(data[1:window])/window
    end 
    for i in (window+1):(n-window) 
        ma[i]= sum(data[i-window:i+window])/(2*window)
    end
    for i in (n-window):n 
        ma[i]= sum(data[i-2*window:i])/(2*window)
    end
    return ma
end


function plot_times(n_init, n_end, title)
    plot_it = collect(n_init:n_end)

    fig, axs = plt.subplots(nrows=nbetas, ncols=nalphas, figsize=(11,11))
    for i in 1:nalphas
        a = alphas[i]
        for j in 1:nbetas
            b = betas[j]
            axs[i,j].set_title("α=$a, β = $b")

            if haskey(bounds_dict,(a,b))
                t_primal, t_dual, t_inner, t_io = times_dict[(a,b)]
                # '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
                axs[i,j].plot(plot_it, moving_average(t_dual[n_init:n_end])
                #./moving_average(t_primal[n_init:n_end])
                ,   label="Dual", color="C1")
                #axs[i,j].plot(plot_it, t_inner[n_init:n_end],  label="Philpott", color="C2")
                axs[i,j].plot(plot_it, moving_average(t_io[n_init:n_end])
                #./moving_average(t_primal[n_init:n_end])
                ,     label="Baucke", color="C3")
                #axs[i,j].plot(plot_it, moving_average(t_primal[n_init:n_end]), label="Primal", color="C5")
                #println( moving_average(t_io[n_init:n_end])./moving_average(t_primal[n_init:n_end])[end])
            end
        end
        axs[1,nbetas].legend()
    end

    fig.suptitle("Times\nModel: $dir - Lipschitz factor: $lip_factor")
    fig.tight_layout()
    plt.savefig(joinpath("data", "output", cfg["save_path"], "$title.pdf"))
end

plot_times(2, 300, "times")
