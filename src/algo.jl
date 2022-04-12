" Randomly choose an index according to probability prob"
function choose(prob; norm=1.0)
  if norm == 0
    return rand(1:length(prob))
  end
  v = rand()*norm
  acc = 0
  for (j,p) in enumerate(prob)
    acc += p
    if acc > v
      return j
    end
  end
  return length(prob)
end

"Make a forward primal pass starting from state0"
function forward(stages, state0; debug=0,return_traj =false)
  if return_traj
    traj=[state0]
  end
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0)
    opt_recover(stage, "primal_fw", "Primal, forward: Failed to solve for $(i)-th stage.\nInitial state $(state0)")

    j = choose(stage.ext[:prob])
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    return_traj && push!(traj,state0)
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0")
      if debug > 1
        println("                         local decision: ", JuMP.value.(stage[:y]))
      end
    end
  end
  return_traj && return traj
end


"Modify objective function of first dual stage to account for initial state decision"
function init_dual(stages, x0)
  m1 = stages[1]
  o1 = JuMP.objective_function(m1)
  JuMP.set_objective_function(m1, o1 - m1[:π0]'*x0)
end

""" Make a foward pass in the dual

epsilon is a regularization parameter that guarantee selection of all nodes we positive probability
normalize the cut to γ=1 if γ > 0.
"""
function forward_dual(stages; debug=0, normalize=false, solvertol=1e-5, epsilon = 1e-2)
  ϵ = epsilon 

  gamma0 = 1.0
  state0 = Float64[]
  for (i,stage) in enumerate(stages)
    if i > 1
      set_initial_state!(stage, state0, gamma0)
    end
    opt_recover(stage, "dual_fw", "Dual, forward: Failed to solve for $(i)-th stage.\nInitial state $(state0), $(gamma0)")

    # Regularize probabilities, so that there's a chance to sample every scenario
    gammas = JuMP.value.(stage.ext[:vars][2])
    gammas_r = gammas .+ ϵ*gamma0
    j = choose(gammas_r, norm=sum(gammas_r))
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    gamma0 = gammas[j]
    if normalize && (gamma0 < solvertol)
      if debug > 0
        println("Setting γ0 = $(gamma0) to zero, below solver tolerance $(solvertol).")
        println("  Scenario variables:", JuMP.value.(stage[:λ][:,j]))
        println("                    :", JuMP.value.(stage[:ξ][:,j]))
        println("                    :", JuMP.value.(stage[:π][:,j]))
      end
      gamma0 = 0.0
    end
    if normalize && (gamma0 > solvertol)
      if debug > 1
        println("Normalizing γ0 = $(gamma0) to one")
      end
      state0 ./= gamma0
      gamma0   = 1.0
    end
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0, prob $gamma0")
      if debug > 1
        println("                  full   local decisions: ", JuMP.value.(stage[:λ]))
        println("                                        : ", JuMP.value.(stage[:ζ]))
        println("                                        : ", JuMP.value.(stage[:ξ]))
        println("  ", gammas)
      end
    end
  end
  if gamma0 == 0
    println("Finished at gamma = 0")
  end
end

""" Add a primal cut 

stage is the current stage problem, at which the cut is added
next is the next stage problem whose dual variable is obtained to define the cut 

the state at which the cut is computed is stored in next.ext
"""
function add_cut!(stage, next)
  ref = JuMP.FixRef.(next.ext[:vars][2])
  x0 = JuMP.value.(ref)
  multipliers = JuMP.dual.(ref)
  cst = JuMP.dual_objective_value(next)

  z = stage.ext[:vars][5]
  x = stage.ext[:vars][1]

  n_scen = length(z)
  for j = 1:n_scen
    JuMP.@constraint(stage, z[j] >= cst + multipliers'*(x[:,j] .- x0))
  end

  # Save cut coefficients
  push!(stage.ext[:cuts], PrimalCut(cst, multipliers, x0))
end


""" Debug function """
function compute_cut(stage, next, next_state::Vector{Float64})
  set_initial_state!(next, next_state)
  opt_recover(next, "primal_cut", "Primal, risk-cut: Failed to solve for initial state $(next_state)")
  ref = JuMP.FixRef.(next.ext[:vars][2])
  multipliers = JuMP.dual.(ref)
  cst = JuMP.dual_objective_value(next)

  return [cst, multipliers, next_state, cst - multipliers'*next_state]
end

""" Add a dual cut 

stage is the current stage problem, at which the cut is added
next is the next stage problem whose dual variable is obtained to define the cut 

the state at which the cut is computed is stored in next.ext
"""
function add_cut_dual!(stage, next)
  ref_π = JuMP.FixRef.(next.ext[:vars][3])
  ref_γ = JuMP.FixRef(next.ext[:vars][4])
  π0 = JuMP.value.(ref_π)
  γ0 = JuMP.value(ref_γ)
  mul_π = JuMP.dual.(ref_π)
  mul_γ = JuMP.dual(ref_γ)
  cst = JuMP.dual_objective_value(next)

  z = stage.ext[:vars][8]
  π = stage.ext[:vars][1]
  γ = stage.ext[:vars][2]

  n_scen = length(z)
  for j = 1:n_scen
    JuMP.@constraint(stage, z[j] >= mul_π'*π[:,j] + mul_γ*γ[j])
  end

  # Save cut coefficients
  push!(stage.ext[:cuts], DualCut(cst, mul_π, mul_γ, π0, γ0))
end

""" Perform a primal backward pass """
function backward(stages)
  for i in 1:(length(stages)-1)
    add_cut!(stages[i], stages[i+1])
  end
end

""" Perform a dual backward pass """
function backward_dual(stages)
  for i in 1:(length(stages)-1)
    add_cut_dual!(stages[i], stages[i+1])
  end
end

""" 
Solve the problem through primal SDDP

M is a MSLBO representing the problem
nstages is the horizon of the problem
risk is a function building the risk measure
state0 is the initial state
niters is the number of iterations ran before stopping

"""
function primalsolve(M::MSLBO, nstages, risk, solver, state0, niters;
                     verbose=false, nprint=10)
  pb = mk_primal_decomp(M, nstages, risk)
  for m in pb
    JuMP.set_optimizer(m, solver)
  end

  println("********")
  println(" PRIMAL ")
  println("********")
  trajs = []
  times = Float64[]
  lbs = Float64[]
  for i = 1:niters
    dt = @elapsed traj = forward(pb, state0; return_traj=true)
    push!(trajs, traj)
    lb = JuMP.objective_value(pb[1])
    push!(lbs, lb)
    
    if verbose && (i % nprint == 0)
      println("Iteration $i: LB = ", lb)
    end
    dt += @elapsed backward(pb)
    push!(times, dt)
  end
  if verbose
    println()
  else
    println("Lower bound: ", lbs[end])
  end

  return pb, trajs, lbs, times
end


""" Compute a primal upper bound "à la" Philpott et. al.
that is using backward propagation through given trajectories
at given niters iteration
"""
function primalub(M, nstages, risk,solver, trajs,niters::Int;verbose=false)
  println("******************************************")
  println(" PRIMAL Upper Bounds at $niters iterations")
  println("******************************************")
  stages = mk_primal_decomp(M, nstages, risk)
  for m in stages
    JuMP.set_optimizer(m, solver)
  end
  Ubs = convex_ub(stages, trajs)
  if verbose
    println("Recursive upper bounds on $niters trajectories")
    display(Ubs)
  else
    println("Upper bound: ", maximum(Ubs[:,0]))
  end
  return Ubs
end

""" Compute a primal upper bound "à la" Philpott et. al.
that is using backward propagation through given trajectories
at iteration defined by iterator niters"""
function primalub(M, nstages, risk,solver,trajs,niters;verbose = false)
  if verbose
    println("********")
    println(" PRIMAL UB")
    println("********")
  end
  ub = Tuple{Int,Float64}[]
  times = Float64[]
  for n in niters
    dt = @elapsed begin
    stages = mk_primal_decomp(M, nstages, risk)
    for m in stages
      JuMP.set_optimizer(m, solver)
    end
    Ubs = convex_ub(stages, trajs[1:n])
    end
    verbose && println("Primal (Philpott) UB at iteration $n: $(Ubs[1,1])")
    push!(ub,(n,Ubs[1,1]))
    push!(times, dt)
  end
  return ub, times
end


""" 
Solve the problem through dual SDDP

M is a MSLBO
nstages is the horizon of the problem
risk is a function building the risk measure
state0 is the initial state
niters is the number of iterations ran before stopping
"""
function dualsolve(M::MSLBO, nstages, risk, solver, state0, niters; verbose=false, nprint=10, epsilon = 1e-2)
  pb = mk_dual_decomp(M, nstages, risk)
  for m in pb
    JuMP.set_optimizer(m, solver)
  end
  println("********")
  println("  DUAL  ")
  println("********")
  init_dual(pb, state0)
  ubs = Float64[]
  times = Float64[]
  for i = 1:niters
    dt = @elapsed forward_dual(pb; normalize=true, epsilon = epsilon)
    ub = -JuMP.objective_value(pb[1])
    push!(ubs, ub)
    if verbose && (i % nprint == 0)
      println("Iteration $i: D-UB = ", ub)
    end
    dt += @elapsed backward_dual(pb)
    push!(times, dt)
  end
  if verbose
    println()
  else
    println("Upper bound: ", ubs[end])
  end

  return pb, ubs, times
end
