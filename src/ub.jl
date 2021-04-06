import JuMP

function bellman_convex_ub(stage,xs, xs_next,zs)
  # stage should be the optimization problem without any cut
  # x0 is a collection of points in the current step at which the new upper bound is computed
  # xs is a collection of states at next step
  # zs is the associated upper value at xs

  z = stage.ext[:vars][5] # value next stage
  x = stage.ext[:vars][1]
  n_scen = length(z)
  n_traj = length(zs)
  JuMP.@variable(stage,s[k=1:n_traj,j=1:n_scen] >= 0) # convex combination coeff

  for j = 1:n_scen #TODO add regularization
    JuMP.@constraint(stage, z[j] >= sum([s[k,j]*zs[k] for k in 1:n_traj]))
    JuMP.@constraint(stage, x[:,j] .== sum([s[k,j]*xs_next[k] for k in 1:n_traj]))
    JuMP.@constraint(stage, sum([s[k,j] for k in 1:n_traj]) == 1)
  end

  ubs = []
  for x0 in xs
    set_initial_state!(stage, x0)
    JuMP.optimize!(stage)
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
      xs = [traj[l][i] for l in 1:n_traj]
      xs_next = [traj[l][i+1] for l in 1:n_traj]
      stage = stages[i]
      Ubs[:,i] = bellman_convex_ub(stage,xs, xs_next,Ubs[:,i+1])
    end
    return Ubs
end
