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
function forward(pb, state0; debug=0,return_traj =false)
  traj = [state0]
  branches = Int64[]

  for (i,stage) in enumerate(pb)
    set_initial_state!(stage, state0)
    opt_recover(stage, "primal_fw", "Primal, forward: Failed to solve for $(i)-th stage.\nInitial state $(state0)")

    j = choose(stage.ext[:prob])
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    push!(traj, state0)
    push!(branches, j)
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0")
      if debug > 1
        println("                         local decision: ", JuMP.value.(stage[:y]))
      end
    end
  end

  push!(pb.ext[:trajs], traj)
  push!(pb.ext[:branches], branches)
  return_traj && return traj
end


"Modify objective function of first dual stage to account for initial state decision"
function init_dual(pb, x0)
  m1 = pb[1]
  o1 = JuMP.objective_function(m1)
  JuMP.set_objective_function(m1, o1 - m1[:π0]'*x0)
end

function normalize_state(γ0, π0, solvertol, stage; debug=0)
  if γ0 < solvertol
    if debug > 0
      println("Setting γ0 = $(γ0) to zero, below solver tolerance $(solvertol).")
      println("    Control dual var:", JuMP.value.(stage[:ξ]))
      println("  Scenario variables:", JuMP.value.(stage[:λ][:,j]))
      println("                    :", JuMP.value.(stage[:π][:,j]))
    end
    return (0.0, π0)
  else
    if debug > 1
      println("Normalizing γ0 = $(γ0) to one")
    end
    return (1.0, π0 ./ γ0)
  end
end

""" Make a foward pass in the dual

`ϵ` is a regularization parameter that guarantees selection of all nodes
with positive probability
if `normalize`, we set `γ=1` (and `π` accordingly) if `γ > solvertol`.
"""
function forward_dual(pb; debug=0, normalize=false, solvertol=1e-5, ϵ = 1e-2)
  traj = DualState[]
  branches = Int64[]

  gamma0 = 1.0
  state0 = Float64[]
  for (i,stage) in enumerate(pb)
    if i > 1
      set_initial_state!(stage, state0, gamma0)
    end
    opt_recover(stage, "dual_fw", "Dual, forward: Failed to solve for $(i)-th stage.\nInitial state $(state0), $(gamma0)")
    if i == 1
        traj = [(1, JuMP.value.(stage.ext[:vars][3]))]
    end

    # Study Lipschitz constant constraints
    L = stage.ext[:lip]
    nx, nscen = size(stage[:π])
    iteration = stage.ext[:iteration]
    opt_pis = JuMP.value.(stage[:π])
    opt_gammas = JuMP.value.(stage[:γ])
    for j in 1:nscen
      for ix in 1:nx
        π_ij = opt_pis[ix,j]
        if isapprox(π_ij,  L*opt_gammas[j]) || isapprox(π_ij, -L*opt_gammas[j])
          stage.ext[:info][(iteration, ix, j)] = π_ij
        end
      end
    end
    # Regularize probabilities, so that there's a chance to sample every scenario
    gammas = JuMP.value.(stage.ext[:vars][2])
    gammas_r = gammas .+ ϵ*gamma0
    j = choose(gammas_r, norm=sum(gammas_r))
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    gamma0 = gammas[j]
    if normalize
      gamma0, state0 = normalize_state(gamma0, state0, solvertol, stage; debug)
    end
    push!(traj, (gamma0, state0))
    push!(branches, j)
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
    println("Finished at γ = 0")
  end
  push!(pb.ext[:trajs], traj)
  push!(pb.ext[:branches], branches)
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

""" Update primal value functions in forward fashion """
function update_vf(pb)
  for i in 1:(length(pb)-1)
    add_cut!(pb[i], pb[i+1])
  end
end

""" Perform a primal backward pass, re-solving problems with new cuts """
function backward(pb)
  for i in (length(pb)-1):-1:1
    opt_recover(pb[i+1], "primal_bw", "Primal, backward: Failed to solve for $(i+1)-th stage.")
    add_cut!(pb[i], pb[i+1])
  end
end

""" Update dual value functions in forward fashion """
function update_vf_dual(pb)
  for i in 1:(length(pb)-1)
    add_cut_dual!(pb[i], pb[i+1])
  end
end

""" Perform a dual backward pass, re-solving problems with new cuts """
function backward_dual(pb)
  for i in (length(pb)-1):-1:1
    opt_recover(pb[i+1], "dual_bw", "Dual, backward: Failed to solve for $(i+1)-th stage.")
    add_cut_dual!(pb[i], pb[i+1])
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
                     verbose=false, nprint=10, backward_solve=true)
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
    if backward_solve
      dt += @elapsed backward(pb)
    else
      dt += @elapsed update_vf(pb)
    end
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
  pb = mk_primal_decomp(M, nstages, risk)
  for m in pb
    JuMP.set_optimizer(m, solver)
  end
  Ubs = convex_ub(pb, trajs)
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
    pb = mk_primal_decomp(M, nstages, risk)
    for m in pb
      JuMP.set_optimizer(m, solver)
    end
    Ubs = convex_ub(pb, trajs[1:n])
    end
    verbose && println("Iteration $n: (Philpott) UB = $(Ubs[1,1])")
    push!(ub,(n,Ubs[1,1]))
    push!(times, dt)
  end
  verbose && println()
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
function dualsolve(M::MSLBO, nstages, risk, solver, state0, niters;
                   verbose=false, nprint=10, backward_solve=true, ϵ = 1e-2)
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
    for m in pb
      m.ext[:iteration] = i
    end
    dt = @elapsed forward_dual(pb; normalize=true, ϵ)
    ub = -JuMP.objective_value(pb[1])
    push!(ubs, ub)
    if verbose && (i % nprint == 0)
      println("Iteration $i: D-UB = ", ub)
    end
    if backward_solve
      dt += @elapsed backward_dual(pb)
    else
      dt += @elapsed update_vf_dual(pb)
    end
    push!(times, dt)
  end
  if verbose
    println()
  else
    println("Upper bound: ", ubs[end])
  end

  return pb, ubs, times
end
