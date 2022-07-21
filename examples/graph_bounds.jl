import Pkg
Pkg.activate(".")

import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, load

include("postproc_funs.jl")

# Linestyles
# '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'

function plot_step(ax, times, values;
                   range::Union{Nothing,Tuple{Int,Int}}=nothing, kwargs...)
  if range != nothing
    idxs = (range[1] .<= times .<= range[2])
    times = times[idxs]
    values = values[idxs]
  end
  n = length(times)
  xs = repeat(times, 1, 2)'  |> (x -> reshape(x, 1,2n))
  ys = repeat(values, 1, 2)' |> (x -> reshape(x, 1,2n))
  ax.plot(xs[2:end], ys[1:end-1]; kwargs...)
end

function bounds(ax, xs, values)
  lb, ub, inner, inner_it, io_lb, io_ub = values
  n_init, n_end = xs[1], xs[end]

  ax.plot(xs, ub[n_init:n_end], label="Dual UB", color="C1")
  plot_step(ax, inner_it, inner, range=(n_init,n_end), label="Philpott UB", color="C2")
  ax.plot(xs, io_ub[n_init:n_end], label="Baucke UB", linestyle="--", color="C3")
  ax.plot(xs, io_lb[n_init:n_end], label="Baucke LB", linestyle="--",color="C4")
  ax.plot(xs, lb[n_init:n_end], label="SDDP LB", color="C5")
end

function rel_gaps(ax, xs, values)
  lb, ub, inner, inner_it, io_lb, io_ub = values
  n_init, n_end = xs[1], xs[end]

  gap    = ub[n_init:n_end]./lb[n_init:n_end] .- 1
  gap_io = io_ub[n_init:n_end]./io_lb[n_init:n_end] .- 1

  ax.plot(xs, gap, label="Dual gap", color="C1")
  ax.plot(xs, gap_io, label="Baucke gap", color="C3")
  # ax[:set_ylim]([0,0.4])
end

function sm_alfabeta(n_init, n_end, bounds_dict, eachone, title, filename; extrapar=nothing)
    if extrapar != nothing
        alphas, betas, params = sets_in_dict(bounds_dict)
    else
        alphas, betas = sets_in_dict(bounds_dict)
    end
    nalphas = length(alphas)
    nbetas  = length(betas)

    plot_it = collect(n_init:n_end)
    fig, axs = plt.subplots(nrows=nbetas, ncols=nalphas, figsize=(11,11))
    for i in 1:nalphas
        a = alphas[i]
        for j in 1:nbetas
            b = betas[j]
            axs[i,j].set_title("α=$a, β=$b")

            key = extrapar == nothing ? (a,b) : (a,b,extrapar)
            if haskey(bounds_dict, key)
                eachone(axs[i,j], plot_it, bounds_dict[key])
            end
        end
    end
    axs[1,1].legend()

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
bounds_dict = read_ablip_all(cfg)

# Base parameters
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = 1

#
# Graphs
title = "Model: $dir - Lipschitz factor: $lip_factor"
outdir = joinpath("data", "output", cfg["save_path"])

# Gaps
filename = joinpath(outdir, "RelGap_beg_Lip$(lip_factor).pdf")
sm_alfabeta(1, 50, bounds_dict, rel_gaps, title, filename; extrapar=lip_factor)
filename = joinpath(outdir, "RelGap_end_Lip$(lip_factor).pdf")
sm_alfabeta(50, 200, bounds_dict, rel_gaps, title, filename; extrapar=lip_factor)

# Bounds
filename = joinpath(outdir, "Bounds_beg_Lip$(lip_factor).pdf")
sm_alfabeta(1, 50, bounds_dict, bounds, title, filename; extrapar=lip_factor)
filename = joinpath(outdir, "Bounds_end_Lip$(lip_factor).pdf")
sm_alfabeta(50, 200, bounds_dict, bounds, title, filename; extrapar=lip_factor)
