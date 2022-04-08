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

function experiment(cfg::ConfigManager, M::MSLBO, state0::Vector{Float64})
    # Risk Aversion
    ra     = cfg["risk-aversion"]
    alpha  = ra["alpha"]
    beta   = ra["beta"]

    risk      = mk_primal_avar(alpha; beta=beta)
    risk_dual = mk_copersp_avar(alpha; beta=beta)

    # Other parameters
    params  = cfg["parameters"]
    nstages = params["nstages"]

    # Solution algorithms
    # Pure dual
    seed!(3)
    dual_pb, dual_ubs = dualsolve(M, nstages, risk_dual, solver, state0, params["dual_iters"]; verbose=true)

    # Primal with interior bounds
    seed!(2)
    primal_pb, primal_trajs, primal_lbs = primalsolve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true)

    # Compute ub from primal trajectories "à la Philpott et al."
    ub_step = params["ub_step"]
    iters_ub = ub_step:ub_step:params["primal_iters"]
    ubs_p = primalub(M, nstages, risk, solver, primal_trajs, iters_ub; verbose=true)

    # Primal with inner and outer bounds
    seed!(4)
    io_pb, io_lbs, io_ubs = problem_child_solve(M, nstages, risk, solver, state0, params["primal_iters"]; verbose=true)


    # Saving info
    data = Dict()
    #   Bounds
    data["primal lb"] = primal_lbs
    data["dual ub"]   = dual_ubs
    data["inner recursive iters"] = first.(ubs_p)
    data["inner recursive bound"] = last.(ubs_p)
    data["io lb"] = io_lbs
    data["io ub"] = io_ubs

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

for idx=1:N
    parse!(cfg, idx)

    # Remove references to specific "module" name?
    experiment(cfg, Hydro_Hist.M, Hydro_Hist.inivol)
end

