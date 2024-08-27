module Hydro1d

include("hydro_scen.jl")

# inivol
inivol = [83.222]

using DualSDDP: MSLBO

"""
Build the MSLBO struct for the hydro problem

Parameters:
- `A::Vector{Matrix{Float64}}`: A matrix at each stage
- `B::Vector{Matrix{Float64}}`: B matrix at each stage
- `T::Vector{Matrix{Float64}}`: T matrix at each stage
- `c::Vector{Vector{Float64}}`: c vector at each stage
- `d::Vector{Vector{Float64}}`: d vector at each stage
- `Ux::Float64`: Upper bound on the positive state x
- `Uy::Vector{Float64}`: Upper bound on the positive control y
- `lb::Float64`: Lower bound at each stage
- `ub::Float64`: Upper bound on the value of the problem
- `Lip::Float64`: Upper bound on the Lipschitz constant
- `prob::Vector{Float64}`: reference probability over branches
- `n_stages::Int=4`: number of stages for testing before running the algorithm
"""
function build(A::Vector{Matrix{Float64}},
                B::Vector{Matrix{Float64}},
                T::Vector{Matrix{Float64}},
                c::Vector{Vector{Float64}},
                d::Vector{Vector{Float64}},
                Ux::Float64,
                Uy::Vector{Float64},
                lb::Float64,
                ub::Float64,
                Lip::Float64,
                prob::Vector{Float64},
                n_stages::Int=4)
  """
  Returns the A matrix at stage t,
  Corresponds to the matrix A in the equation Ax_t + Bx_{t-1} + Ty = d

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `A[t]::Matrix{Float64}`: A matrix at stage t
  """
  function A_func(t::Int, i::Int)
    return A[t]
  end

  """
  Returns the B matrix at stage t,
  Corresponds to the matrix B in the equation Ax_t + Bx_{t-1} + Ty = d

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `B[t]::Matrix{Float64}`: B matrix at stage t
  """
  function B_func(t::Int, i::Int)
    return B[t]
  end

  """
  Returns the T matrix at stage t,
  Corresponds to the matrix T in the equation Ax_t + Bx_{t-1} + Ty = d

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `T[t]::Matrix{Float64}`: T matrix at stage t
  """
  function T_func(t::Int, i::Int)
    return T[t]
  end
  
  """
  Returns the c vector at stage t,
  Corresponds to the marginal cost of y_t

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `c[t]::Vector{Float64}`: c vector at stage t
  """
  function c_func(t::Int, i::Int)
    return c[t]
  end

  """
  Returns the d vector at stage t,
  Corresponds to the vector of demmand (d) in the equation Ax_t + Bx_{t-1} + Ty = d

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `d::Vector{Float64}`: d vector at stage t
  """
  function d_func(t::Int, i::Int)
    return d[t]
  end

  """
  Returns the Upper bound on the positive state x
  Observations:
  It convert the scalar Ux to a vector of size 1 for access in the MSLBO
  
  Parameters:
  - `t::Int`: stage

  Returns:
  - `Ux::Vector{Float64}`: Upper bound on the positive state x
  """
  function Ux_func(t::Int)
    return [Ux]
  end

  """
  Returns the Upper bound on the positive control y

  Parameters:
  - `t::Int`: stage

  Returns:
  - `Uy::Vector{Float64}`: Upper bound on the positive control y
  """
  function Uy_func(t::Int)
    return Uy
  end

  """
  Returns the Lower bound at each stage

  Parameters:
  - `t::Int`: stage

  Returns:
  - `lb::Float64`: Lower bound at each stage
  """
  function lb_builder(t::Int, nstages::Int)
    return lb
  end

  """
  Returns the Upper bound on the value of the problem

  Parameters:
  - `t::Int`: stage

  Returns:
  - `ub::Float64`: Upper bound on the value of the problem
  """
  function ub_builder(t::Int, nstages::Int)
    return ub
  end

  """
  Returns the Upper bound on the Lipschitz constant

  Parameters:
  - `t::Int`: stage

  Returns:
  - `Lip::Float64`: Upper bound on the Lipschitz constant
  """
  function Lip_builder(t::Int, nstages::Int)
    return Lip
  end

  """
  Returns the reference probability over branches

  Parameters:
  - `t::Int`: stage

  Returns:
  - `prob::Vector{Float64}`: reference probability over branches
  """
  function prob_builder(t::Int)
    if t == 1
      return prob
    else
      return ones(Main.nscen)/Main.nscen
    end
  end

  if size(A, 1) != n_stages
    error("The number of stages must be equal to the number of A matrices")
  elseif size(B, 1) != n_stages
    error("The number of stages must be equal to the number of B matrices")
  elseif size(T, 1) != n_stages
    error("The number of stages must be equal to the number of T matrices")
  elseif size(c, 1) != n_stages
    error("The number of stages must be equal to the number of c vectors")
  elseif size(d, 1) != n_stages
    error("The number of stages must be equal to the number of d vectors")
  end

  return MSLBO(A_func, B_func, T_func, c_func, d_func,
                Ux_func, Uy_func, lb_builder, ub_builder,
                Lip_builder, prob_builder)
end

A = [reshape([1.0, 0], 2, 1), reshape([1.0, 0], 2, 1),
      reshape([1.0, 0], 2, 1), reshape([1.0, 0], 2, 1)]

B = [reshape([-1.0, 0], 2, 1), reshape([-1.0, 0], 2, 1),
      reshape([-1.0, 0], 2, 1), reshape([-1.0, 0], 2, 1)]

T = [reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5)]

c = [[0.1, 1, 5, 10, 50], [0, 1, 5, 10, 50],
      [0, 1, 5, 10, 50], [0, 1, 5, 10, 50]]

d = [[0, 75.0], [0, 75.0], [0, 75.0], [0, 75.0]]

Ux = 100.0
Uy = [60.0, 200, 15, 15, 75]

lb = 0.0
ub = 75*50.0
Lip = 50.0
prob = [1.0]

M = build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, 4)

end