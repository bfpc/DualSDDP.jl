module Build

using DualSDDP: MSLBO

"""
Build the MSLBO struct for the hydro problem by passing vetor of matrices

Parameters:
- `A::Vector{Matrix{Float64}}`: A matrix at each stage
- `B::Vector{Matrix{Float64}}`: B matrix at each stage
- `T::Vector{Matrix{Float64}}`: T matrix at each stage
- `c::Vector{Vector{Float64}}`: marginal cost of y_t at each stage
- `d::Vector{Vector{Float64}}`: demmand (d) at each stage
- `Ux::Vector{Float64}`: Upper bound on the positive state x (same for all stages)
- `Uy::Vector{Float64}`: Upper bound on the positive control y (same for all stages)
- `lb::Float64`: Lower bound
- `ub::Float64`: Upper bound on the value of the problem
- `Lip::Float64`: Upper bound on the Lipschitz constant
- `prob::Vector{Float64}`: reference probability over branches (Fix!!!)
- `n_stages::Int=4`: number of stages for testing before running the algorithm

Returns:
- `MSLBO`: MSLBO struct for the hydro problem, with the functions defined by
 the parameters on the builder
"""
function build(A::Vector{Matrix{Float64}},
                B::Vector{Matrix{Float64}},
                T::Vector{Matrix{Float64}},
                c::Vector{Vector{Float64}},
                d::Vector{Vector{Float64}},
                Ux::Vector{Float64},
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
    return Ux
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
    if t == 1
      return prob
    else
      return ones(Main.nscen)/Main.nscen
    end
  end

  for (matrix, name) in zip([A, B, T,c ,d], ["A", "B", "T", "c", "d"])
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
- `d::Vector{Float64}`: demmand (d) at the stage
- `Ux::Float64`: Upper bound on the positive state x
- `Uy::Float64`: Upper bound on the positive control y
- `lb::Float64`: Lower bound
- `ub::Float64`: Upper bound on the value of the problem
- `Lip::Float64`: Upper bound on the Lipschitz constant
- `prob::Vector{Float64}`: reference probability over branches
"""
struct StageMLSBO
    A::Matrix{Float64}
    B::Matrix{Float64}
    T::Matrix{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    Ux::Float64
    Uy::Float64
    lb::Float64
    ub::Float64
    Lip::Float64
    prob::Vector{Float64}
end

"""
Build the MSLBO struct for the hydro problem by passing each stage as a struct

Parameters:
- `stages::Vector{StageMLSBO}`: Vector of stages for the MSLBO problem

Returns:
- `MSLBO`: MSLBO struct for the hydro problem, with the functions defined by
 the parameters on the builder
"""
function build(stages::Vector{StageMLSBO})
  A = [stage.A for stage in stages]
  B = [stage.B for stage in stages]
  T = [stage.T for stage in stages]
  c = [stage.c for stage in stages]
  d = [stage.d for stage in stages]
  Ux = stages[1].Ux
  Uy = stages[1].Uy
  lb = stages[1].lb
  ub = stages[1].ub
  Lip = stages[1].Lip
  prob = stages[1].prob
  return build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, length(stages))
end
  
  end