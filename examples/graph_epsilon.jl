import Pkg
Pkg.activate(".")

using Plots
using Config: parse!, ConfigManager, total_combinations, load

plot_gap=1e-7

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

N = total_combinations(cfg)

# These are "unparametrizable" yet
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]

bounds_dict = Dict()
epsilons = Set()

lb = 1e15

for idx=1:N
    println(idx)
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]
    epsilon = cfg["parameters"]["epsilon"]

    bounds_dict[epsilon] = (data["dual ub"])
    lb = min(data["dual ub"][end],lb)
    #lb = 1.883e8
    push!(epsilons, epsilon)
end
lb = (1-plot_gap)*lb

p = plot()
for epsilon in epsilons
    plot!(p,bounds_dict[epsilon].-lb, yaxis=:log, label=string(epsilon))
end
plot(p)

for k in 0:6
    epsilon = 1/10^k
    println("ε=$epsilon yields final gap with best lb= ",bounds_dict[epsilon][end]-lb/(1-plot_gap))
end

