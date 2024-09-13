"""
Build the MSLBO struct by passing vectors of matrices.

  The (risk-neutral) MSLBO corresponds to the following problem:
  min. sum_{t=1}^{T} E [ < c_t | y_t > ]
  s.t.
  A_t x_t + B_t x_{t-1} + T_t y_t = d_t^j
  0 <= x_t <= Ux
  0 <= y_t <= Uy

  This interface is less general than the Functional form in the struct MSLBO,
  but can be easier to use in simple problems.
  In particular, it only deals with uncertainties in the RHS of the constraints,
  and the bounds on the state and control are constant over time.

Parameters:
- `A::Vector{Matrix{Float64}}`: A matrix at each stage
- `B::Vector{Matrix{Float64}}`: B matrix at each stage
- `T::Vector{Matrix{Float64}}`: T matrix at each stage
- `c::Vector{Vector{Float64}}`: marginal cost of y_t at each stage
- `d::Vector{Vector{Vector{Float64}}}`: d vector at each stage and branch j
- `Ux::Vector{Vector{Float64}}`: Upper bound on the positive state x at each stage
- `Uy::Vector{Vector{Float64}}`: Upper bound on the positive control y at each stage
- `lb::Vector{Float64}`: Lower bound on the value of the problem at each stage
- `ub::Vector{Float64}`: Upper bound on the value of the problem at each stage
- `Lip::Vector{Float64}`: Upper bound on the Lipschitz constant at each stage
- `prob::Vector{Vector{Float64}}`: reference probability at each stage over branches
- `n_stages::Int`: number of stages of the problem
"""
function build(A::Vector{Matrix{Float64}},
                B::Vector{Matrix{Float64}},
                T::Vector{Matrix{Float64}},
                c::Vector{Vector{Float64}},
                d::Vector{Vector{Vector{Float64}}},
                Ux::Vector{Vector{Float64}},
                Uy::Vector{Vector{Float64}},
                lb::Vector{Float64},
                ub::Vector{Float64},
                Lip::Vector{Float64},
                prob::Vector{Vector{Float64}},
                n_stages::Int)
  function A_func(t::Int, i::Int)
    return A[t]
  end

  function B_func(t::Int, i::Int)
    return B[t]
  end

  function T_func(t::Int, i::Int)
    return T[t]
  end

  function c_func(t::Int, i::Int)
    return c[t]
  end

  function d_func(t::Int, i::Int)
    return d[t][i]
  end

  function Ux_func(t::Int)
    return Ux[t]
  end

  function Uy_func(t::Int)
    return Uy[t]
  end

  function lb_func(t::Int, nstages::Int)
    return lb[t]
  end

  function ub_func(t::Int, nstages::Int)
    return ub[t]
  end

  function Lip_func(t::Int, nstages::Int)
    return Lip[t]
  end

  function prob_func(t::Int)
    prob[t]
  end

  vectors = [A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob]
  names = ["A", "B", "T", "c", "d", "Ux", "Uy", "lb", "ub", "Lip", "prob"]
  for (vector, name) in zip(vectors, names)
    if length(vector) != n_stages
      error("The number of stages ($n_stages) must be equal to the number of entries in $name; $(length(vector)) were given.")
    end
  end

  return MSLBO(A_func, B_func, T_func, c_func, d_func,
                Ux_func, Uy_func, lb_func, ub_func,
                Lip_func, prob_func)
end

"""
Struct to hold one Linear Bellman Operator, with only RHS uncertainties.

  The LBO builds the recourse problems Q(x, j) defined by:
  Q(x, j) =
    min.  < c | y > + FutureCost(z)
    s.t.  A z + B x + T y = d^j
          0 <= x <= Ux
          0 <= y <= Uy

  The value function  Q(x)  will only be determined after choosing a risk measure.

Fields:
- `A::Matrix{Float64}`: A matrix at the stage
- `B::Matrix{Float64}`: B matrix at the stage
- `T::Matrix{Float64}`: T matrix at the stage
- `c::Vector{Float64}`: marginal cost of y at the stage
- `d::Vector{Vector{Float64}}`: d vector at the stage (demand) over branches
- `Ux::Vector{Float64}`: Upper bound on the positive state x
- `Uy::Vector{Float64}`: Upper bound on the positive control y
- `lb::Float64`: Lower bound on the value of the problem at the stage
- `ub::Float64`: Upper bound on the value of the problem at the stage
- `Lip::Float64`: Upper bound on the Lipschitz constant
- `prob::Vector{Float64}`: reference probability over branches
"""
struct SimpleLBO
    A::Matrix{Float64}
    B::Matrix{Float64}
    T::Matrix{Float64}
    c::Vector{Float64}
    d::Vector{Vector{Float64}}
    Ux::Vector{Float64}
    Uy::Vector{Float64}
    lb::Float64
    ub::Float64
    Lip::Float64
    prob::Vector{Float64}
end

"""
Build the MSLBO struct from the list of SimpleLBO stages.

Parameters:
- `stages::Vector{SimpleLBO}`: Vector of stages for the MS problem

Returns:
- `MSLBO`: MSLBO in functional form
"""
function build(stages::Vector{SimpleLBO})
  A = [stage.A for stage in stages]
  B = [stage.B for stage in stages]
  T = [stage.T for stage in stages]
  c = [stage.c for stage in stages]
  d = [stage.d for stage in stages]
  Ux = [stage.Ux for stage in stages]
  Uy = [stage.Uy for stage in stages]
  lb = [stage.lb for stage in stages]
  ub = [stage.ub for stage in stages]
  Lip = [stage.Lip for stage in stages]
  prob = [stage.prob for stage in stages]
  return build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, length(stages))
end
