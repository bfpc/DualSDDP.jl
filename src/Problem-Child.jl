# The DP model for a multi-stage stochastic linear program with SWI noise
#
# V_t(x_{t-1}) = min. \rho [ c_t^j y_t^j + V_{t+1}(x_t^j) ]
#                s.t. A_t^j x_t^j + B_t^j x_{t-1} + T_t^j y_t^j = d_t^j
#                     0 <= x_t^j <= Ux_t
#                     0 <= y_t^j <= Uy_t
#
# where \rho is the Risk Measure
#

""" Primal cuts"""
struct Cut
    intercept :: Float64
    slope :: Vector{Float64}
end

""" Inner approximation vertices"""
struct Vertex
    value :: Float64
    point :: Vector{Float64}
end

""" Object representing both inner and outer problem 

outer is a collection of JuMP models representing the outer problem for each branch
inner is a collection of JuMP models representing the inner problem for each branch
prob is the reference probability
n_branches is the number of branch per stage
risk is a function building a JuMP model
cuts is a collection of primal cuts
inner_vertices is a collection of vertices inside the epigraph
ub_model is a JuMP model to compute the upper bound
"""
struct IO_stage
    outer :: Vector{Model}
    inner :: Vector{Model}
    prob :: Vector{Float64}
    n_branches :: Int
    risk :: Function
    cuts :: Vector{Cut}
    inner_vertices :: Vector{Vertex}
    ub_model :: JuMP.Model
end
# TODO: constructor

import Base.show
function show(io::IO, ::MIME"text/plain", s::DualSDDP.IO_stage)
  println(io, "A $(s.n_branches)-branch stage,")
  println(io, "  with probabilities $(s.prob)")
  print(  io, "  and state dimension $(length(s.outer[1][:x]))")
end

""" Initialize the IO_stages from the model M with risk function risk"""
function mk_primal_io(M::MSLBO, T::Int, risk)
    stages = IO_stage[]
    for t in 1:T
        outer = mk_primal_outer(M,t,T)
        inner = mk_primal_inner(M,t,T)
        prob = M.prob(t)
        n_branches = length(prob)
        risk = risk
        ub_model = mk_primal_ub_model(M, t, T)
        if t == T
            for m in inner
                JuMP.set_objective_coefficient(m,m[:_z],0)
            end
            m = ub_model
            JuMP.set_objective_coefficient(m,m[:_z],0)
        end
        stage = IO_stage(outer,inner,prob,n_branches,risk,Cut[],Vertex[],ub_model)
        push!(stages,stage)
    end
    return stages
end

""" Construct the initial primal outer problem"""
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

""" Construct the initial upper bound problem"""
function mk_primal_ub_model(M::MSLBO, t::Int, T::Int)
    nx = size(M.A(t,1),2)
    m = Model()
    @variable(m, x[i=1:nx])
    @variable(m,0 <= σ0 <= 1)
    @variable(m,_z)

    δ = @variable(m, δ[1:nx])
    δ_abs = @variable(m,δ_abs[1:nx] >= 0)

    @constraint(m,δ_abs.>= δ)
    @constraint(m,δ_abs.>= -δ)
    @constraint(m,z_lb, σ0*M.ub(t,T) +  M.Lip(t,T) * sum(δ_abs) <= _z)
    @constraint(m,x_cc, δ .== x)
    @constraint(m,σ_cc, σ0 == 1)

    @objective(m, Min, _z)
    return m
end

""" Construct the initial inner problem"""
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
        @variable(m, _z)

        δ = @variable(m, δ[1:nx])
        δ_abs = @variable(m,δ_abs[1:nx] >= 0)
        @constraint(m,δ_abs.>= δ)
        @constraint(m,δ_abs.>= -δ)

        @variable(m,0 <= σ0 <=1)
        @constraint(m,z_lb, M.Lip(t,T) * sum(δ_abs) + σ0*M.ub(t,T) <= _z)
        @constraint(m,x_cc, δ .== x)
        @constraint(m,σ_cc, σ0 == 1)

        @objective(m, Min, M.c(t,j)'*y + _z)

        m.ext[:vars] = (x, x0, y, _z)
        m.ext[:lip]  = M.Lip(t, T)
        push!(scen_probs, m)
    end
    return scen_probs
end

""" Compute, for a given state value, the gap between the inner and outer approximation
required by problem-child selection procedure
"""
function gap(stage::IO_stage,state)
    lb = -Inf
    for c in stage.cuts
        intercept = c.intercept
        slope = c.slope
        lb = max(lb,intercept + slope'*state)
    end

    m = stage.ub_model
    JuMP.fix.(m[:x],state)
    opt_recover(m, "problem_child_gap", "Primal: Failed to compute upper bound at $(state)")
    ub = JuMP.objective_value(m)
    return ub - lb
end

""" Return a vector of all possible next next states computed from curr_states
using the outer approximation """
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

""" Compute the problem-child next state 

Given a curr state and current inner and outer cost-to-go update_approximations
compute the next states according to all branches using outer cost-to-go 
return the next state with the largest difference between current upper and lower approx
"""
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

"""
Make a forward phase according to the problem child strategy
return trajectory 
"""
function forward(stages::Vector{IO_stage}, state0; debug=0)
    curr_state = state0
    traj = [curr_state]
    for (t,stage) in enumerate(stages)
        curr_state = choose_problem_child(stage,curr_state)
        push!(traj,curr_state)
        if debug > 0
            println("Going out from stage $t, state $curr_state")
        end
    end
    return traj
end

"""
Update inner and outer approximation along a given trajectories 
(e.g. computed by forward) 
"""
function update_approximations(stages::Vector{IO_stage},traj,solver; backward_solve=true)
    T = length(stages)
    if backward_solve
        stage_list = T-1:-1:1
    else
        stage_list = 1:T-1
    end
    for t in stage_list
        next_state = traj[t+1]
        c = compute_cut(stages[t],stages[t+1],next_state,solver)
        add_cut!(stages[t], c)
        v = compute_vertex(stages[t],stages[t+1],next_state,solver)
        add_vertex!(stages[t], v)
    end
end

""" Compute a new inner vertex"""
function compute_vertex(stage::IO_stage, next_stage::IO_stage,next_state, solver)
    values = []
    for next in next_stage.inner
        set_initial_state!(next, next_state)
        opt_recover(next, "primal_vertex_pc", "Primal, vertex: Failed to solve for initial state $(next_state)")
        value = JuMP.objective_value(next)
        push!(values,value)
    end

    # change of measure
    risk = next_stage.risk
    m = Model(solver)
    risk(m,values,next_stage.prob)
    JuMP.optimize!(m)
    ra_value = JuMP.objective_value(m)

    # Save vertex coefficients
    v = Vertex(ra_value,next_state)
    push!(stage.inner_vertices, v)
    return v
end

""" Calculate  ρ(Z) """
function compute_risk(values::Array, probabilities::Array, risk, solver)
    m = Model(solver)
    risk(m, values, probabilities)
    JuMP.optimize!(m)
    return JuMP.objective_value(m)
end

""" Add the same vertices in every inner model (one per branch) """
function add_vertex!(stage::IO_stage, vertex::Vertex)
    for m in stage.inner
        add_vertex!(m,vertex)
    end
    add_vertex!(stage.ub_model,vertex)
end

""" Add a vertex to an inner model"""
function add_vertex!(m::JuMP.Model, vertex::Vertex)
    σk = @variable(m)
    set_lower_bound.(σk,0)
    set_upper_bound.(σk,1)
    set_normalized_coefficient(m[:σ_cc], σk, 1)

    xk = vertex.point
    for i in 1:length(xk)
        set_normalized_coefficient(m[:x_cc][i], σk, xk[i])
    end

    vk = vertex.value
    set_normalized_coefficient(m[:z_lb], σk, vk)
end

""" Compute a primal cut at next_state"""
function compute_cut(stage::IO_stage, next_stage::IO_stage,next_state, solver)
    slopes = []
    values = []
    for next in next_stage.outer
        set_initial_state!(next, next_state)
        opt_recover(next, "primal_cut_pc", "Primal, cut: Failed to solve for initial state $(next_state)")
        ref = JuMP.FixRef.(next.ext[:vars][2])
        multipliers = JuMP.dual.(ref)
        value = JuMP.dual_objective_value(next)
        push!(slopes,multipliers)
        push!(values,value)
    end

    n = next_stage.n_branches
    # change of measure
    risk = next_stage.risk
    m = Model(solver)
    @variable(m, t[1:n])
    @constraint(m, delta, t .>= values)
    risk(m,t,next_stage.prob)
    JuMP.optimize!(m)
    ra_value = JuMP.objective_value(m)
    γs = JuMP.dual.(delta)
    ra_slope = sum(γs[i]*slopes[i] for i in 1:n)
    ra_intercept = ra_value - ra_slope'*next_state

    # Save cut coefficients
    c = Cut(ra_intercept,ra_slope)
    push!(stage.cuts, c)
    return c
end

function add_cut!(stage::IO_stage, cut::Cut)
    for m in stage.outer
        #TODO : finding variables by name
        z = m.ext[:vars][4]
        x = m.ext[:vars][1]
        JuMP.@constraint(m, z >= cut.intercept + cut.slope'*x)
    end
end

"""Solve the problem M using the problem-child approach

M is an MSLBO 
nstages is the number of stages considered
risk is a function building the risk function
solver is a linear solver 
state0 is the initaial state 
niters is the number of iterations ran before stopping
"""
function problem_child_solve(M, nstages, risk, solver, state0, niters;verbose=false, nprint=10, backward_solve=true)
    pb = mk_primal_io(M, nstages, risk)
    for stage in pb
        for m in stage.inner
            JuMP.set_optimizer(m, solver)
        end
        for m in stage.outer
            JuMP.set_optimizer(m, solver)
        end
        JuMP.set_optimizer(stage.ub_model, solver)
    end

    println("********************")
    println(" PROBLEM CHILD ")
    println("********************")
    trajs = []
    times = Float64[]
    lbs = Float64[]
    p_ubs = Float64[]
    for i = 1:niters
        dt = @elapsed traj = forward(pb, state0)
        push!(trajs, traj)
        # Calculate LowerBound from first-stage outer optimal values
        first_stage = pb[1]
        values = [JuMP.objective_value(sp) for sp in first_stage.outer]
        lb = compute_risk(values, first_stage.prob, first_stage.risk, solver)
        push!(lbs, lb)
        if verbose && (i % nprint == 0)
            print("Iteration $i: LB = ", lb)
        end

        # Add cuts/vertices to inner/outer approximations at the states
        # in the forward trajectory
        dt += @elapsed update_approximations(pb,trajs[end],solver; backward_solve)

        # Calculate UpperBound from first-stage inner optimal values
        values = Float64[]
        for m in first_stage.inner
            set_initial_state!(m, state0)
            dt += @elapsed opt_recover(m, "problem_child_1_st", "Primal, upper bound: Failed to solve initial problem")
            push!(values, JuMP.objective_value(m))
        end
        p_ub = compute_risk(values, first_stage.prob, first_stage.risk, solver)
        push!(p_ubs,p_ub)
        if verbose && (i % nprint == 0)
            println(" P-UB = ", p_ub)
        end
        push!(times,dt)
    end
    if verbose
        println()
    else
        println("Lower bound: ", lbs[end]," Upper bounds: ", p_ubs[end])
    end
    return pb, lbs, p_ubs, times
end

"""
    write_policy_to_file(model::Vector{IO_stage}, filename::String)

Write cuts and verticies that form the policy in `model` to `filename`,
in JSON format, similar to `SDDP.write_cuts_to_file`.
"""
function write_policy_to_file(model::Vector{IO_stage}, filename::String)
    policy = Dict{String,Any}[]
    for (i, stage) in enumerate(model)
        node_info = Dict(
            "node" => string(i),
            "single_cuts" => Dict{String,Any}[],
            "vertices" => Dict{String,Any}[],
        )

        for cut in stage.cuts
          push!(node_info["single_cuts"], to_dict(cut))
        end
        for vertex in stage.inner_vertices
          push!(node_info["vertices"], to_dict(vertex))
        end

        push!(policy, node_info)
    end
    open(filename, "w") do io
        return write(io, JSON.json(policy))
    end
    return
end

function to_dict(cut::Cut)
    Dict(
         "intercept" => cut.intercept,
         "coefficients" => Dict(["x$(n)" => slope_n for (n, slope_n) in enumerate(cut.slope)]),
        )
  end

  function to_dict(v::Vertex)
    Dict(
         "value" => v.value,
         "state" => Dict(["x$(n)" => x_n for (n, x_n) in enumerate(v.point)]),
        )
  end
