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

function save_vfs!(data, primal_pb, dual_pb)
    nstages = length(primal_pb)

    # Primal
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

    # Dual
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

function experiment(cfg::ConfigManager, M::MSLBO, state0::Vector{Float64};
    primal = true, dual = true, philpott_ub = true, pb_child = true)
    # Risk Aversion
    ra     = cfg["risk-aversion"]
    alpha  = ra["alpha"]
    beta   = ra["beta"]
    epsilon = cfg["parameters"]["epsilon"]

    primal = cfg["parameters"]["primal"]
    dual = cfg["parameters"]["dual"]
    philpott_ub = cfg["parameters"]["philpott_ub"]
    pb_child = cfg["parameters"]["pb_child"]

    risk      = mk_primal_avar(alpha; beta=beta)
    risk_dual = mk_copersp_avar(alpha; beta=beta)

    # Other parameters
    params  = cfg["parameters"]
    nstages = params["nstages"]

    # Solution algorithms
    # Pure dual
    if dual
      println(epsilon)
      seed!(3)
      dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, params["dual_iters"]; verbose=true, epsilon = epsilon)
    else
      dual_pb, dual_ubs, dual_times = 0,0,0
    end

    # Primal with interior bounds

    if (primal || philpott_ub)
      seed!(2)
      primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true)

      if philpott_ub
        # Compute ub from primal trajectories "à la Philpott et al."
        ub_step = params["ub_step"]
        iters_ub = ub_step:ub_step:params["primal_iters"]
        ubs_p, ubs_times = primalub(M, nstages, risk, solver, primal_trajs, iters_ub; verbose=true)
      else
        ubs_p, ubs_times = 0,0
      end
    else
      primal_pb, primal_trajs, primal_lbs, primal_times = 0,0,0,0
      ubs_p, ubs_times = 0,0
    end

    # Primal with inner and outer bounds
    if pb_child
      seed!(4)
      io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true)
    else
      io_pb, io_lbs, io_ubs, io_times = 0,0,0,0
    end


    # Saving info
    data = Dict()
    #   Bounds
    data["primal lb"] = primal_lbs
    data["primal t"]  = primal_times
    data["dual ub"]   = dual_ubs
    data["dual t"]    = dual_times
    data["inner recursive iters"] = first.(ubs_p)
    data["inner recursive bound"] = last.(ubs_p)
    data["inner recursive t"] = ubs_times
    data["io lb"] = io_lbs
    data["io ub"] = io_ubs
    data["io t"]  = io_times

    save_vfs!(data, primal_pb, dual_pb)

    save(cfg, data)
end

cfg_name = ARGS[1]
cfg = ConfigManager(cfg_name, @__DIR__)

N = total_combinations(cfg)

# These are "unparametrizable" yet
parse!(cfg, 1)

dir = cfg["experiment"]["dir"]
lip_factor = cfg["parameters"]["Lip"]
include(joinpath(dir, cfg["experiment"]["Model"]))

ttime = 0

for idx=1:N
    parse!(cfg, idx)

    # Remove references to specific "module" name?
    #print(cfg)
    α = cfg["risk-aversion"]["alpha"]
    β = cfg["risk-aversion"]["beta"]
    L = cfg["parameters"]["Lip"]
    exp = cfg["experiment"]["dir"]
    nstages = cfg["parameters"]["nstages"]
    println("##################################")
    println("$idx/$N experience $exp with horizon $nstages for α=$α, β=$β and Lip_factor=$L")
    println("##################################")

    global ttime += @elapsed experiment(cfg, Hydro_Hist.M, Hydro_Hist.inivol)
    println()
    println("Estimated time remaining: ", ttime*(N-idx)/idx )
    println()
end

