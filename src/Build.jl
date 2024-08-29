module Build

using DualSDDP: MSLBO

"""
Build the MSLBO struct for the hydro problem by passing vetor of matrices

Parameters:
- `A::Vector{Matrix{Float64}}`: A matrix at each stage
- `B::Vector{Matrix{Float64}}`: B matrix at each stage
- `T::Vector{Matrix{Float64}}`: T matrix at each stage
- `c::Vector{Vector{Float64}}`: marginal cost of y_t at each stage
- `d::Vector{Float64}`: d vector
- `Ux::Vector{Vector{Float64}}`: Upper bound on the positive state x at each stage
- `Uy::Vector{Vector{Float64}}`: Upper bound on the positive control y at each stage
- `lb::Float64`: Lower bound on the value of the problem
- `ub::Float64`: Upper bound on the value of the problem
- `Lip::Float64`: Upper bound on the Lipschitz constant
- `prob::Vector{Vector{Float64}}`: reference probability over branches
- `n_stages::Int=4`: number of stages for testing before running the algorithm

Returns:
- `MSLBO`: MSLBO struct for the hydro problem, with the functions defined by
 the parameters on the builder
"""
function build(A::Vector{Matrix{Float64}},
                B::Vector{Matrix{Float64}},
                T::Vector{Matrix{Float64}},
                c::Vector{Vector{Float64}},
                d::Vector{Float64},
                Ux::Vector{Vector{Float64}},
                Uy::Vector{Vector{Float64}},
                lb::Float64,
                ub::Float64,
                Lip::Float64,
                prob::Vector{Vector{Float64}},
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
  Returns the d vector,
  Corresponds to the vector of demmand (d) in the equation Ax_t + Bx_{t-1} + Ty = d

  Parameters:
  - `t::Int`: stage
  - `i::Int`: branch

  Returns:
  - `d::Vector{Float64}`: d vector
  """
  function d_func(t::Int, i::Int)
    return d
  end

  """
  Returns the Upper bound on the positive state x at stage t
  
  Parameters:
  - `t::Int`: stage

  Returns:
  - `Ux::Vector{Float64}`: Upper bound on the positive state x at stage t
  """
  function Ux_func(t::Int)
    return Ux[t]
  end

  """
  Returns the Upper bound on the positive control y at stage t

  Parameters:
  - `t::Int`: stage

  Returns:
  - `Uy::Vector{Float64}`: Upper bound on the positive control y at stage t
  """
  function Uy_func(t::Int)
    return Uy[t]
  end

  """
  Returns the Lower bound at each stage

  Parameters:
  - `t::Int`: stage

  Returns:
  - `lb::Float64`: Lower bound at each stage
  """
  function lb_func(t::Int, nstages::Int)
    return lb
  end

  """
  Returns the Upper bound on the value of the problem

  Parameters:
  - `t::Int`: stage

  Returns:
  - `ub::Float64`: Upper bound on the value of the problem
  """
  function ub_func(t::Int, nstages::Int)
    return ub
  end

  """
  Returns the Upper bound on the Lipschitz constant

  Parameters:
  - `t::Int`: stage

  Returns:
  - `Lip::Float64`: Upper bound on the Lipschitz constant
  """
  function Lip_func(t::Int, nstages::Int)
    return Lip
  end

  """
  Returns the reference probability over branches

  Parameters:
  - `t::Int`: stage

  Returns:
  - `prob::Vector{Float64}`: reference probability over branches
  """
  function prob_func(t::Int)
    prob[t]
  end

  for (matrix, name) in zip([A, B, T, c], ["A", "B", "T", "c"])
    if size(matrix, 1) != n_stages
      error("The number of stages ($n_stages) must be equal to the number of $name, $(size(matrix, 1)) were given")
    end
  end

  return MSLBO(A_func, B_func, T_func, c_func, d_func,
                Ux_func, Uy_func, lb_func, ub_func,
                Lip_func, prob_func)
end

"""
Struct to hold one stage of the MSLBO problem

Fields:
- `A::Matrix{Float64}`: A matrix at the stage
- `B::Matrix{Float64}`: B matrix at the stage
- `T::Matrix{Float64}`: T matrix at the stage
- `c::Vector{Float64}`: marginal cost of y_t at the stage
- `Ux::Float64`: Upper bound on the positive state x at the stage
- `Uy::Float64`: Upper bound on the positive control y at the stage
- `prob::Vector{Float64}`: reference probability over branches
"""
struct StageMLSBO
    A::Matrix{Float64}
    B::Matrix{Float64}
    T::Matrix{Float64}
    c::Vector{Float64}
    Ux::Vector{Float64}
    Uy::Vector{Float64}
    prob::Vector{Float64}
end

"""
Build the MSLBO struct for the hydro problem by passing each stage as a struct

Parameters:
- `stages::Vector{StageMLSBO}`: Vector of stages for the MSLBO problem
- `d::Vector{Float64}`: d vector
- `lb::Float64`: Lower bound on the value of the problem
- `ub::Float64`: Upper bound on the value of the problem
- `Lip::Float64`: Upper bound on the Lipschitz constant

Returns:
- `MSLBO`: MSLBO struct for the hydro problem, with the functions defined by
 the parameters on the builder
"""
function build(stages::Vector{StageMLSBO}, d::Vector{Float64}, lb::Float64,
                               ub::Float64, Lip::Float64)
  A = [stage.A for stage in stages]
  B = [stage.B for stage in stages]
  T = [stage.T for stage in stages]
  c = [stage.c for stage in stages]
  Ux = [stage.Ux for stage in stages]
  Uy = [stage.Uy for stage in stages]
  prob = [stage.prob for stage in stages]
  return build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, length(stages))
end
  
  end