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
nscens = Set()

for idx=1:N
    println(idx)
    parse!(cfg, idx)
    data = load(cfg)

    nscen = cfg["parameters"]["nscen"]

    times_dict[(nscen)] = (data["primal t"], data["dual t"], data["inner recursive t"], data["io t"])
    bounds_dict[(nscen)] = (data["primal lb"], data["dual ub"], data["inner recursive bound"],data["inner recursive iters"],data["io lb"],data["io ub"])

    push!(nscens, nscen)
end

nscens = sort([nscens...])

nnscens = length(nscens)

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

    fig, axs = plt.subplots(nrows=nnscens, figsize=(4,11))
    for i in 1:nnscens
        n = nscens[i]
            
            axs[i].set_title("nscen=$n")

            if haskey(bounds_dict,(n))
                t_primal, t_dual, t_inner, t_io = times_dict[(n)]
                # '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
                axs[i].plot(plot_it, moving_average(t_dual[n_init:n_end])
                ./moving_average(t_primal[n_init:n_end])
                ,   label="Dual", color="C1")
                #axs[i,1].plot(plot_it, t_inner[n_init:n_end],  label="Philpott", color="C2")
                axs[i].plot(plot_it, moving_average(t_io[n_init:n_end])
                ./moving_average(t_primal[n_init:n_end])
                ,     label="Baucke", color="C3")
                #axs[i,1].plot(plot_it, moving_average(t_primal[n_init:n_end]), label="Primal", color="C5")
                #println( moving_average(t_io[n_init:n_end])./moving_average(t_primal[n_init:n_end])[end])
            end
        
        #axs[1,nbetas].legend()
    end

    fig.suptitle("Times\nModel: $dir - Lipschitz factor: $lip_factor")
    fig.tight_layout()
    plt.savefig(joinpath("data", "output", cfg["save_path"], "$title.pdf"))
end


## Generating table of iteration time
it = 100 
for n in [10,20,40,80]
           t_primal, t_dual, t_inner, t_io = times_dict[(n)]
           println(n," & ", round(moving_average(t_primal[1:300])[it],digits=3) 
           ," & ", round(moving_average(t_dual[1:300])[it],digits=3) 
           ," & ", round(moving_average(t_io[1:300])[it],digits=3) 
           , "\\\\" )
end

plot_times(2, 300, "times")
