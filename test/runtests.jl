using Test
using Random: seed!
using MathOptInterface: OptimizerWithAttributes, Silent
using HiGHS: Optimizer

solver = OptimizerWithAttributes(Optimizer, Silent() => true)

using .DualSDDP

# A trivial one-stage problem
module OneStage
using ..DualSDDP: MSLBO

demand = 1.0:10.0
nscen = length(demand)

function A(t::Int, i::Int)
  return reshape([1.0, 0], 2, 1)
end

function B(t::Int, i::Int)
  return reshape([-1.0, 0], 2, 1)
end

function T(t::Int, i::Int)
  return [0.0 0; 1 -1]
end

function c(t::Int, i::Int)
  return [1.0, 0]
end

function d(t::Int, i::Int)
  return [0.0, demand[i]]
end

function Ux(t::Int)
  return [100.0]
end
function Uy(t::Int)
  return [10.0, 10]
end

function lb(t::Int, nstages::Int)
  return 0
end
function ub(t::Int, nstages::Int)
  return (nstages - t + 1)*10
end

function Lip(t::Int, nstages::Int)
  return (nstages - t + 1)
end

function prob(t::Int)
  return ones(nscen)/nscen
end

M = MSLBO(A,B,T,c,d,Ux,Uy,lb,ub,Lip,prob)

end

# Test
function test_risk1()
  alpha = 0.5
  beta = 0.5
  risk_primal = mk_primal_avar(alpha; beta)
  risk_dual = mk_copersp_avar(alpha; beta)
  risk_obj = 6.75
  target = [risk_obj, risk_obj]

  M = OneStage.M
  state0 = [0.0]

  # Pure primal
  seed!(1)
  nstages, niters = 1, 2
  primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk_primal, solver, state0, niters; verbose=false)
  @test primal_lbs ≈ target

  # Pure dual
  seed!(1)
  dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=false)
  @test dual_ubs ≈ target

  # Recursive upper bounds over primal trajectories
  rec_ubs, rec_times = primalub(M, nstages, risk_primal, solver, primal_trajs, [1,2]; verbose=false)
  @test last.(rec_ubs) ≈ target

  # Primal with outer and inner bounds
  seed!(1)
  io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk_primal, solver, state0, niters; verbose=false)
  @test io_lbs ≈ target
  @test io_ubs ≈ target

  nothing
end

test_1stage()

