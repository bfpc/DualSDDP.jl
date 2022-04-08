# Need M::MSLBO, nstages::Int, risk, risk_dual, state0::Vector{Float64},
#      niters::Int, solver
# Solution algorithms
# Pure primal
seed!(2)
primal_pb, primal_trajs, primal_lbs, primal_times = primalsolve(M, nstages, risk, solver, state0, niters; verbose=true)

# Pure dual
seed!(1)
dual_pb, dual_ubs, dual_times = dualsolve(M, nstages, risk_dual, solver, state0, niters; verbose=true)

# Recursive upper bounds over primal trajectories
rec_ubs, rec_times = primalub(M, nstages, risk, solver, primal_trajs, ub_step:ub_step:niters; verbose=true)

# Primal with outer and inner bounds
seed!(1)
io_pb, io_lbs, io_ubs, io_times = problem_child_solve(M, nstages, risk, solver, state0, niters; verbose=true)

