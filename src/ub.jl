# Used to show previous states before solver failure
function debug_print(xs, x0)
  function inner()
    println("Previous states:")
    for xi in xs
      if xi == x0
        break
      end
      println(xi)
    end
  end
end

function bellman_convex_ub(stage,xs, xs_next,zs)
  # stage should be the optimization problem without any cut
  # xs is a collection of points in the current step at which the new upper bound is computed
  # xs_next is a collection of states at next step
  # zs is the associated upper value at xs

  z = stage.ext[:vars][5] # value next stage
  x = stage.ext[:vars][1]
  state_dim = size(x,1)
  L = stage.ext[:lip]
  n_scen = length(z)
  n_traj = length(zs)
  @variable(stage,s[k=1:n_traj,j=1:n_scen] >= 0) # convex combination coeff
  reg = @variable(stage,[1:state_dim,1:n_scen]) # regularization for each scenario

  regabs = @variable(stage,[1:state_dim,1:n_scen]) # modeling absolute value
  @constraint(stage, regabs .>=  reg)
  @constraint(stage, regabs .>= -reg)
  for j = 1:n_scen
    @constraint(stage, z[j] >= sum([s[k,j]*zs[k] for k in 1:n_traj]) + L*sum(regabs[:,j]))
    @constraint(stage, x[:,j] .== reg[:,j] + sum([s[k,j]*xs_next[k] for k in 1:n_traj]))
    @constraint(stage, sum([s[k,j] for k in 1:n_traj]) == 1)
  end

  ubs = []
  for x0 in xs
    set_initial_state!(stage, x0)
    opt_recover(stage, "ub_inner", "Failed to solve for initial state $(x0)",
                extra_print=debug_print(xs, x0))
    push!(ubs, JuMP.objective_value(stage))
  end
  return (ubs)
end

function convex_ub(stages,traj)
  #TODO : assuming final cost to be 0
    n_traj = length(traj)
    T = length(stages)

    Ubs = zeros((n_traj,T+1))

    for i in T:-1:1
      #print("  Evaluating at stage $(i): ")
      xs = [traj[l][i] for l in 1:n_traj]
      xs_next = [traj[l][i+1] for l in 1:n_traj]
      stage = stages[i]
      Ubs[:,i] = bellman_convex_ub(stage,xs, xs_next,Ubs[:,i+1])
      #println()
    end
    return Ubs
end
