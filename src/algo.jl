import JuMP

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

function forward(stages, state0; debug=0,return_traj =false)
  if return_traj
    traj=[state0]
  end
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0)
    JuMP.optimize!(stage)

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



function init_dual(stages, x0)
  m1 = stages[1]
  o1 = JuMP.objective_function(m1)
  JuMP.set_objective_function(m1, o1 - m1[:π0]'*x0)
end

function forward_dual(stages; debug=0, normalize=false)
  ϵ = 1e-1

  gamma0 = 1.0
  state0 = Float64[]
  for (i,stage) in enumerate(stages)
    if i > 1
      set_initial_state!(stage, state0, gamma0)
    end
    JuMP.optimize!(stage)

    # Regularize probabilities, so that there's a chance to sample every scenario
    gammas = JuMP.value.(stage.ext[:vars][2])
    gammas_r = gammas .+ ϵ*gamma0
    if debug > 0
      println("  ", gammas)
    end
    j = choose(gammas_r, norm=sum(gammas_r))
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    gamma0 = gammas[j]
    if normalize && (gamma0 != 0)
      state0 ./= gamma0
      gamma0   = 1
    end
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0, prob $gamma0")
      if debug > 1
        println("                         local decisions: ", JuMP.value.(stage[:λ]))
        println("                                        : ", JuMP.value.(stage[:ζ]))
        println("                                        : ", JuMP.value.(stage[:ξ]))
      end
    end
  end
end

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
end

function add_cut_dual!(stage, next)
  ref_π = JuMP.FixRef.(next.ext[:vars][3])
  ref_γ = JuMP.FixRef.(next.ext[:vars][4])
  π0 = JuMP.value.(ref_π)
  γ0 = JuMP.value(ref_γ)
  mul_π = JuMP.dual.(ref_π)
  mul_γ = JuMP.dual.(ref_γ)
  cst = JuMP.dual_objective_value(next)

  z = stage.ext[:vars][8]
  π = stage.ext[:vars][1]
  γ = stage.ext[:vars][2]

  n_scen = length(z)
  for j = 1:n_scen
    JuMP.@constraint(stage, z[j] >= cst + mul_π'*(π[:,j] .- π0) + mul_γ*(γ[j] - γ0))
  end
end

function backward(stages)
  for i in 1:(length(stages)-1)
    add_cut!(stages[i], stages[i+1])
  end
end

function backward_dual(stages)
  for i in 1:(length(stages)-1)
    add_cut_dual!(stages[i], stages[i+1])
  end
end

function primalsolve(M, nstages, risk, solver, state0, niters;
                     verbose=false, ub=false)
  pb = mk_primal_decomp(M, nstages, risk)
  for m in pb
    JuMP.set_optimizer(m, solver)
  end

  println("********")
  println(" PRIMAL ")
  println("********")
  trajs = []
  for i = 1:niters
    push!(trajs, forward(pb, state0; return_traj=true))
    backward(pb)
    verbose && println("Iteration $i: LB = ", JuMP.objective_value(pb[1]))
  end
  if verbose
    println()
  else
    println("Lower bound: ", JuMP.objective_value(pb[1]))
  end

  if ub
    println("********************")
    println(" PRIMAL Upper Bounds")
    println("********************")
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
    return pb, trajs, stages, Ubs
  else
    return pb, trajs
  end
end

function dualsolve(M, nstages, risk, solver, state0, niters; verbose=false)
  pb = mk_dual_decomp(M, nstages, risk)
  for m in pb
    JuMP.set_optimizer(m, solver)
  end
  println("********")
  println("  DUAL  ")
  println("********")
  for i = 1:niters
    forward_dual(pb)
    backward_dual(pb)
    println("Iteration $i: UB = ", -JuMP.objective_value(pb[1]))
  end
  if verbose
    println()
  else
    println("Upper bound: ", -JuMP.objective_value(pb[1]))
  end

  return pb
end
