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
    ub_model :: JuMP.Model
end

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

function update_approximations(stages::Vector{IO_stage},traj,solver)
    T = length(stages)
    for t in 1:T-1
        next_state = traj[t+1]
        c = compute_cut(stages[t],stages[t+1],next_state,solver)
        add_cut!(stages[t], c)
        v = compute_vertex(stages[t],stages[t+1],next_state,solver)
        add_vertex!(stages[t], v)
    end
end

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

function add_vertex!(stage::IO_stage, vertex::Vertex)
    for m in stage.inner
        add_vertex!(m,vertex)
    end
    add_vertex!(stage.ub_model,vertex)
end

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


function problem_child_solve(M, nstages, risk, solver, state0, niters;
                     verbose=false)
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
    lbs = Float64[]
    p_ubs = Float64[]
    for i = 1:niters
        push!(trajs, forward(pb, state0))
         #TODO missing first stage risk measure
        lb = JuMP.objective_value(pb[1].outer[1])
        push!(lbs, lb)
        verbose && print("Iteration $i: LB = ", lb)
        update_approximations(pb,trajs[end],solver)
        
       

        m = pb[1].inner[1]
        set_initial_state!(m,state0)
        opt_recover(m, "problem_child_1_st", "Primal, Upperbound: Failed to solve initial problem")
        p_ub = JuMP.objective_value(m)
        verbose && println(" P-UB = ", p_ub)
        push!(p_ubs,p_ub)
    end
    if verbose
        println()
    else
        println("Lower bound: ", lbs[end]," Upper bounds: ", p_ubs[end])
    end
    return pb

end
