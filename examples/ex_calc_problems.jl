# Need M::MSLBO, nstages::Int, risk, risk_dual, state0::Vector{Float64},
#      niters::Int, solver
# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, niters; verbose=true)

# Pure primal with only forward pass
seed!(2)
primal_pb_fw, primal_trajs_fw, primal_lbs_fw, primal_times_fw = primalsolve(M, nstages, risk, solver, state0, niters; verbose=true, backward_solve=false)

# Pure dual
seed!(1)
dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=true)

# Pure dual with only forward pass
seed!(1)
dual_pb_fw, dual_ubs_fw, dual_times_fw = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=true, backward_solve=false)

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, niters; verbose=true)

# Primal with outer and inner bounds, only forward pass
seed!(1)
io_pb_fw, io_lbs_fw, io_ubs_fw, io_times_fw = problem_child_solve(M, nstages, risk, solver, state0, niters; verbose=true, backward_solve=false)

