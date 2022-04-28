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


function plot_bounds(n_init, n_end, title)

    fig, axs = plt.subplots(nrows=nbetas, ncols=nalphas)
    for i in 1:nalphas
        a = alphas[i]
        
        for j in 1:nbetas
            b = betas[j]
            axs[i,j].set_title("α=$a, β = $b")
            
            plot_it = collect(n_init:n_end)
            if haskey(bounds_dict,(a,b))  
                lb, ub, inner, inner_it, io_lb, io_ub = bounds_dict[(a,b)]
                #gap = ub./lb .- 1
                # '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
                axs[i,j].plot(plot_it,ub[n_init:n_end], label="Dual UB", color="C1")
                axs[i,j].plot(inner_it, inner, label="Philpott UB", color="C2")
                axs[i,j].plot(plot_it,io_ub[n_init:n_end], label="Baucke UB", linestyle="--", color="C3")
                axs[i,j].plot(plot_it,io_lb[n_init:n_end], label="Baucke LB",  linestyle="--",color="C4")
                axs[i,j].plot(plot_it,lb[n_init:n_end], label="SDDP LB", color="C5")
                # axs[i,j].semilogy(inner_it, inner, linestyle="--", color="C$k")
                # axs[i,j].semilogy( io_ub./io_lb .- 1, linestyle="dashdot", color="C$k")
            end 
            #axs[i,j].legend()
        end
        axs[1,nbetas].legend()
        # axs[i].set_title("α = $a")
    end

    fig.suptitle(dir * " - Lipschitz factor: $lip_factor")
    println(joinpath("data", "output", cfg["save_path"], "$title.pdf"))
    plt.savefig(joinpath("data", "output", cfg["save_path"], "$title.pdf"))
end

plot_bounds(100,500,"end")