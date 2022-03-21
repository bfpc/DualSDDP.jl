# The DP model for a multi-stage stochastic linear program with SWI noise
#
# V_t(x_{t-1}) = min. \rho [ c_t^j y_t^j + V_{t+1}(x_t^j) ]
#                s.t. A_t^j x_t^j + B_t^j x_{t-1} + T_t^j y_t^j = d_t^j
#                     0 <= x_t^j <= Ux_t
#                     0 <= y_t^j <= Uy_t
#
# where \rho is the Risk Measure
#

struct Cut
    intercept :: Float64
    slope :: Vector{Float64}
end

struct Vertex
    value :: Float64
    point :: Vector{Float64}
end

struct IO_stage 
    outer :: Vector{Model}
    inner :: Vector{Model}
    prob :: Vector{Float64}
    n_branches :: Int  
    risk :: Function
    cuts :: Vector{Cut}
    inner_vertices :: Vector{Vertex}
    ub_model
end

function mk_primal_io(M::MSLBO, T::Int, risk)
    stages = IO_stage[]
    for t in 1:T
      outer = mk_primal_outer(M,t,T)
      inner = mk_primal_inner(M,t,T)
      prob = M.prob(t)
      n_branches = length(prob)
      risk = risk

      stage = IO_stage(outer,inner,prob,n_branches,risk,[],[])
      push!(stages,stage)
    end
    return stages 
end

function mk_primal_outer(M::MSLBO,t::Int, T::Int)
      prob = M.prob(t)
      n = length(prob)

      nx = size(M.A(t,1),2)
      nxprev = size(M.B(t,1),2)
      ny = length(M.c(t,1))
      
      scen_probs = Model[]
      for j in 1:n #problem for one child
        m = Model()
  
        x = @variable(m, x[i=1:nx])
        y = @variable(m, y[i=1:ny])
        x0 = @variable(m, x0[i=1:nxprev])
  
        set_lower_bound.(x,0)
        set_upper_bound.(x,M.Ux(t))
        set_lower_bound.(y,0)
        set_upper_bound.(y,M.Uy(t))
  
        @constraint(m, M.A(t,j)*x + M.B(t,j)*x0 + M.T(t,j)*y .== M.d(t,j))
        θ = @variable(m, θ)
        set_lower_bound(θ, M.lb(t,T))
    
        @objective(m, Min, M.c(t,j)'*y + θ)
        
        m.ext[:vars] = (x, x0, y, θ)
        m.ext[:lip]  = M.Lip(t, T)
        push!(scen_probs, m)
      end
    return scen_probs
end

function mk_primal_inner(M::MSLBO, t::Int, T::Int)
    prob = M.prob(t)
    n = length(prob)

    nx = size(M.A(t,1),2)
    nxprev = size(M.B(t,1),2)
    ny = length(M.c(t,1))

    scen_probs = Model[]
    for j in 1:n #problem for one child
      m = Model()

      x = @variable(m, x[i=1:nx])
      y = @variable(m, y[i=1:ny])
      x0 = @variable(m, x0[i=1:nxprev])

      set_lower_bound.(x,0)
      set_upper_bound.(x,M.Ux(t))
      set_lower_bound.(y,0)
      set_upper_bound.(y,M.Uy(t))

      @constraint(m, M.A(t,j)*x + M.B(t,j)*x0 + M.T(t,j)*y .== M.d(t,j))
      _z = @variable(m, _z)
      
      δ = @variable(m, δ[1:nx])
      δ_abs = @variable(m,δ_abs[1:nx] >= 0)
      @constraint(m,δ_abs.>= δ)
      @constraint(m,δ_abs.>= -δ)

      @constraint(m,z_lb,_z >= M.Lip(t,T) * sum(δ_abs))
      @constraint(m,x_cc, x .== δ)
      #@variable(m,0 <= σ[1:nb_iter] <=1)
      #@constraint(m,σ_cc, 1 == sum(σ))

      @objective(m, Min, M.c(t,j)'*y + _z)
      
      m.ext[:vars] = (x, x0, y, _z)
      m.ext[:lip]  = M.Lip(t, T)
      push!(scen_probs, m)
    end
  return scen_probs
end

function gap(stage::IO_stage,state)
    lb = -∞
    for (intercept,slope) in stage.cuts
        lb = max(lb,intercept + slope'*state)
    end

    m = stage.ub_model 
    JuMP.fix.(m[:x],state)
    opt_recover(m, "problem_child_gap", "Primal: Failed to compute upper bound at $(state0)")
    ub = JuMP.objective_value(m)
    return ub - lb
end

function compute_next_states(stage::IO_stage,curr_state)
    next_states=Vector{Float64}[]

    n = stage.n_branches
    for j in 1:n
        m = stage.outer[j]
        set_initial_state!(m,curr_state)
        opt_recover(m, "problem_child_fw", "Primal, forward: Failed to solve for initial state $(curr_state)")
        state_j = JuMP.value.(m.ext[:vars][1])
        push!(next_states,state_j)
    end

    return next_states
end

function choose_problem_child(stage::IO_stage,curr_state)
    next_states = compute_next_states(stage,curr_state)
    worst_gap = -1
    next_state = curr_state
    for state in next_states
        g = gap(stage,state)
        if g > worst_gap
            worst_gap = g
            next_state = state
        end
    end 
    return next_state
end 
  
function forward(stages::Vector{IO_stage}, state0; debug=0)
    curr_state = state0
    for (t,stage) in enumerate(stages)        
        curr_state = choose_problem_child(stage,curr_state)
        
      if debug > 0
        println("Going out from stage $t, state $curr_state")
      end
    end
end

function update_approximations(stages::Vector{IO_stage})
    curr_state = state0
    for (t,stage) in enumerate(stages)        
        curr_state = choose_problem_child(stage,curr_state)
        
      if debug > 0
        println("Going out from stage $t, state $curr_state")
      end
    end
end
  
  
  function compute_cut(stage::IO_stage, next_stage::IO_stage)
    slopes = []
    intercepts = []
    for next in next_stage.outer
        ref = JuMP.FixRef.(next.ext[:vars][2])
        x0 = JuMP.value.(ref)
        multipliers = JuMP.dual.(ref)
        cst = JuMP.dual_objective_value(next)
        intercept = cst - x0'*multipliers
        push!(slopes,multipliers)
        push!(intercepts,intercept)
    end

    #TODO copier change of measure de Oscar


    # Save cut coefficients
     push!(stage.ext[:cuts], Cut(ra_intercept,ra_slope))
  end

  
  function add_cut!(stage, next)
    # TODO a terminer
    for j = 1:n_scen
        z = stage.ext[:vars][5]
        x = stage.ext[:vars][1]
        JuMP.@constraint(stage, z[j] >= cst + multipliers'*(x[:,j] .- x0))
      end
  end
  
  
  function backward(stages)
    for i in 1:(length(stages)-1)
      add_cut!(stages[i], stages[i+1])
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
    lbs = Float64[]
    for i = 1:niters
      push!(trajs, forward(pb, state0; return_traj=true))
      lb = JuMP.objective_value(pb[1])
      push!(lbs, lb)
      verbose && println("Iteration $i: LB = ", lb)
      backward(pb)
    end
    if verbose
      println()
    else
      println("Lower bound: ", lbs[end])
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
      return pb, trajs, lbs, stages, Ubs
    else
      return pb, trajs, lbs
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
        print("  Evaluating at stage $(i): ")
        xs = [traj[l][i] for l in 1:n_traj]
        xs_next = [traj[l][i+1] for l in 1:n_traj]
        stage = stages[i]
        Ubs[:,i] = bellman_convex_ub(stage,xs, xs_next,Ubs[:,i+1])
        println()
      end
      return Ubs
  end
    