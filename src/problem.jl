# The DP model for a multi-stage stochastic linear program with SWI noise
#
# V_t(x_{t-1}) = min. \rho [ c_t^j y_t^j + V_{t+1}(x_t^j) ]
#                s.t. A_t^j x_t^j + B_t^j x_{t-1} + T_t^j y_t^j = d_t^j
#                     0 <= x_t^j <= Ux_t
#                     0 <= y_t^j <= Uy_t
#
# where \rho is the Risk Measure
#
#
# D_t(π_{t-1}, γ_{t-1}) = min. ζ_{t-1} Ux_{t-1} + Σ [ λ_t^j d_t^j + ξ_t^j Uy_t + ]
#                                                   [    + D_{t+1}(π_t^j, γ_t^j) ]
#                 s.t. (γ^j) ∈ γ_{t-1} Q
#                      Σ [ λ_t^j B_t^j ] + ζ_{t-1}     ⩾ π_{t-1}
#                        - λ_t^j A_t^j                 = π_t^j
#                      γ^j c_t^j + λ_t^j T_t^j + ξ_t^j ⩾ d_t^j
#                      0 <= ζ_{t-1}
#                      0 <= ξ_t^j
#
# where Q is the dual polyhedron for the Risk Measure \rho.

""" Build the initial primal JuMP model for each stage"""
function mk_primal_decomp(M::MSLBO, T::Int, risk)
  stages = Model[]
  for t in 1:T
    prob = M.prob(t)

    m = Model()
    n = length(prob)
    nx = size(M.A(t,1),2)
    nxprev = size(M.B(t,1),2)
    ny = length(M.c(t,1))

    x = @variable(m, x[i=1:nx,j=1:n])
    y = @variable(m, y[i=1:ny,j=1:n])
    x0 = @variable(m, x0[i=1:nxprev])

    set_lower_bound.(x,0)
    set_upper_bound.(x,M.Ux(t))
    set_lower_bound.(y,0)
    set_upper_bound.(y,M.Uy(t))

    for j = 1:n
      @constraint(m, M.A(t,j)*x[:,j] + M.B(t,j)*x0 + M.T(t,j)*y[:,j] .== M.d(t,j))
    end

    _t = @variable(m, _t[j=1:n])
    _z = @variable(m, _z[j=1:n])
    set_lower_bound.(_z, M.lb(t,T))

    @constraint(m, [j=1:n], _t[j] >= M.c(t,j)'*y[:,j] + _z[j])
    risk(m, _t, prob)

    m.ext[:vars] = (x, x0, y, _t, _z)
    m.ext[:prob] = prob
    m.ext[:lip]  = M.Lip(t, T)
    m.ext[:cuts] = []
    push!(stages, m)
  end
  return stages
end

""" set the initial state of a given stage"""
function set_initial_state!(stage::Model, x0::Vector{T}) where T <: Real
  fix.(stage.ext[:vars][2], x0)
  return
end

""" Build the initial dual JuMP model for each stage"""
function mk_dual_decomp(M::MSLBO, T::Int, dualrisk)
  stages = Model[]
  for t in 1:T
    prob = M.prob(t)
    m = Model()
    n = length(prob)
    nx = length(M.Ux(t))
    ny = length(M.Uy(t))
    nd = length(M.d(t,1))

    λ = @variable(m, λ[i=1:nd,j=1:n])
    π = @variable(m, π[i=1:nx,j=1:n])
    γ = @variable(m, γ[j=1:n])
    ζ = @variable(m, ζ[i=1:nx])
    ξ = @variable(m, ξ[i=1:ny,j=1:n])
    π0 = @variable(m, π0[i=1:nx])
    γ0 = @variable(m, γ0)

    set_lower_bound.(γ,0)
    set_lower_bound.(ζ,0)
    set_lower_bound.(ξ,0)

    # Box constraints for π, using (primal) Lipschitz constant
    L = M.Lip(t, T)
    @constraint(m, [j=1:n], π[:,j] .>= -L*γ[j])
    @constraint(m, [j=1:n], π[:,j] .<=  L*γ[j])

    for j = 1:n
      @constraint(m, γ[j]*M.c(t,j) + M.T(t,j)'*λ[:,j] + ξ[:,j] .>= 0)
      @constraint(m, - M.A(t,j)'*λ[:,j] .== π[:,j])
    end

    z = @variable(m, z[j=1:n])
    ub = M.ub(t,T)
    @constraint(m, z .>= -ub*γ)

    dualrisk(m, γ, prob, γ0)

    if t == 1
      @objective(m, Min, sum(λ[:,j]'*M.d(t,j) + ξ[:,j]'*M.Uy(t) + z[j] for j=1:n))
      @constraint(m, sum(M.B(t,j)'*λ[:,j] for j=1:n) .== π0)
      fix(γ0, 1)
    else
      @constraint(m, sum(M.B(t,j)'*λ[:,j] for j=1:n) + ζ .>= π0)
      @objective(m, Min, ζ'*M.Ux(t-1) + sum(λ[:,j]'*M.d(t,j) + ξ[:,j]'*M.Uy(t) + z[j] for j=1:n))
    end

    # Value function for the last stage is known
    if t == T
      @variable(m, ζ_end[i=1:nx,j=1:n] >= 0)
      @constraint(m, ζ_end .>= π)
      @constraint(m, [j=1:n], z[j] == ζ_end[:,j]'*M.Ux(t))
    end

    m.ext[:vars] = (π, γ, π0, γ0, λ, ζ, ξ, z)
    m.ext[:prob] = prob
    m.ext[:cuts] = []

    push!(stages, m)
  end
  return stages
end

""" Set initial dual state """
function set_initial_state!(stage::Model, π0::Vector{T}, γ0::Real) where T <: Real
  fix.(stage.ext[:vars][3], π0)
  fix(stage.ext[:vars][4], γ0)
  return
end


