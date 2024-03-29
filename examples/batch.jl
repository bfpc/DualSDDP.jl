import Pkg
Pkg.activate(".")

import JuMP
import Gurobi
import PyPlot as plt
using Config: parse!, ConfigManager, total_combinations, save
using Random: seed!

using DualSDDP
env = Gurobi.Env()
solver = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0)
# import GLPK
# solver = GLPK.Optimizer

nprint = 25 #printing of current results

function save_primal_vfs!(data, primal_pb)
    nstages = length(primal_pb)

    primal_VFs = [pvf_info(stage) for stage in primal_pb]
    primal_iters = length(primal_VFs[1])

    data["primal VF lbs"] = last.(primal_VFs)
    # (t,k) = k-th intercept at t-th stage
    intercepts = Array{Float64}(undef, nstages-1, primal_iters)
    for t = 1:nstages-1
      for k = 1:primal_iters
        intercepts[t,k] = primal_VFs[t][1][k][1]
      end
    end
    data["primal VF intercepts"] = intercepts

    # (t,k,j) = coefficient of j-th variable, in k-th cut, at t-th stage
    stagedim = length(primal_VFs[1][1][1][2])
    multipliers = Array{Float64}(undef, nstages-1, primal_iters, stagedim)
    for t = 1:nstages-1
      for k = 1:primal_iters
        multipliers[t,k,:] .= primal_VFs[t][1][k][2]
      end
    end
    data["primal VF multipliers"] = multipliers

    return
end

function save_dual_vfs!(data, dual_pb)
    nstages = length(dual_pb)

    dual_VFs = [dvf_info(stage) for stage in dual_pb]
    dual_iters = length(dual_VFs[1])

    # (t,k) = k-th intercept at t-th stage
    intercepts = Array{Float64}(undef, nstages-1, dual_iters)
    for t = 1:nstages-1
      for k = 1:dual_iters
        intercepts[t,k] = dual_VFs[t][k][1]
      end
    end
    data["dual VF intercepts"] = intercepts

    # (t,k) = k-th gamma coefficient at t-th stage
    mul_γ = Array{Float64}(undef, nstages-1, dual_iters)
    for t = 1:nstages-1
      for k = 1:dual_iters
        mul_γ[t,k] = dual_VFs[t][k][3]
      end
    end
    data["dual VF mul gamma"] = mul_γ

    # (t,k,j) = coefficient of j-th variable, in k-th cut, at t-th stage
    stagedim = length(dual_VFs[1][1][2])
    mul_π = Array{Float64}(undef, nstages-1, dual_iters, stagedim)
    for t = 1:nstages-1
      for k = 1:dual_iters
        mul_π[t,k,:] .= dual_VFs[t][k][2]
      end
    end
    data["dual VF multipliers"] = mul_π

    return
end

function experiment(cfg::ConfigManager, M::MSLBO, state0::Vector{Float64})
    # Risk Aversion
    ra     = cfg["risk-aversion"]
    alpha  = ra["alpha"]
    beta   = ra["beta"]

    risk      = mk_primal_avar(alpha; beta=beta)
    risk_dual = mk_copersp_avar(alpha; beta=beta)

    # Algorithms to run
    algo = cfg["algorithms"]
    primal      = algo["primal"]
    dual        = algo["dual"]
    philpott_ub = algo["philpott_ub"]
    pb_child    = algo["pb_child"]

    # Other parameters
    params  = cfg["parameters"]
    nstages = params["nstages"]
    # Regularization for probabilities in dual forward
    ϵ = params["epsilon"]

    # Dict of data to be saved
    data = Dict()

    # Solution algorithms
    # Pure dual
    if dual
      seed!(3)
      dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, params["dual_iters"]; verbose=true, ϵ, nprint)
      # Collect bounds, iteration times and value functions
      data["dual ub"] = dual_ubs
      data["dual t"]  = dual_times
      save_dual_vfs!(data, dual_pb)
    end

    # Primal with interior bounds
    if (primal || philpott_ub)
      seed!(2)
      primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true, nprint = nprint)

      # Collect bounds, iteration times and value functions
      data["primal lb"] = primal_lbs
      data["primal t"]  = primal_times
      save_primal_vfs!(data, primal_pb)

      if philpott_ub
        # Compute ub from primal trajectories "à la Philpott et al."
        ub_step = params["ub_step"]
        iters_ub = ub_step:ub_step:params["primal_iters"]
        ubs_p, ubs_times = primalub(M, nstages, risk, solver, primal_trajs, iters_ub; verbose=true)
        # Collect bounds, iterations and iteration times
        data["inner recursive iters"] = first.(ubs_p)
        data["inner recursive bound"] = last.(ubs_p)
        data["inner recursive t"] = ubs_times
      end
    end

    # Primal with inner and outer bounds
    if pb_child
      seed!(4)
      io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true,nprint = nprint)
      # Collect bounds and iteration times
      data["io lb"] = io_lbs
      data["io ub"] = io_ubs
      data["io t"]  = io_times
      # TODO
      # pb_child && save_io_vfs!(data, io_pb)
    end

    # Save collected data to disk
    save(cfg, data)
end

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

# Global configuration
parse!(cfg, 1)

# lip_factor and nscen should be in Main namespace for include() to work
dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]
nscen = cfg["parameters"]["nscen"]
MyModule = include(joinpath(dir, cfg["experiment"]["Model"]))

# If all went well, save the current config file to the output directory
outdir = joinpath("data", "output", cfg["save_path"])
mkpath(outdir)
cfg_save = joinpath(outdir, cfg_name)
cp(cfg_name, cfg_save)
chmod(cfg_save, 0o440)

ttime = 0

N = total_combinations(cfg)
for idx=1:N
    parse!(cfg, idx)

    # Hack to change lip_factor / nscen
    global lip_factor
    global nscen
    global MyModule
    new_lip_factor = cfg["parameters"]["Lip"]
    new_nscen = cfg["parameters"]["nscen"]
    if (new_lip_factor != lip_factor) || (new_nscen != nscen)
      lip_factor = new_lip_factor
      nscen = new_nscen
      MyModule = include(joinpath(dir, cfg["experiment"]["Model"]))
    end

    α = cfg["risk-aversion"]["alpha"]
    β = cfg["risk-aversion"]["beta"]
    nstages = cfg["parameters"]["nstages"]
    ε = cfg["parameters"]["epsilon"]
    println("##################################")
    println("$idx/$N experience $dir with horizon $nstages, $nscen branches/stage, for α=$α, β=$β, Lip_factor=$lip_factor and ε=$ε")
    println("##################################")

    global ttime += @elapsed experiment(cfg, MyModule.M, MyModule.inivol)
    println()
    println("Estimated time remaining: ", ttime*(N-idx)/idx )
    println()
    GC.gc() #freeing memory before next run
end

