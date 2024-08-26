module Hydro1d

include("hydro_scen.jl")

# inivol
inivol = [83.222]

using DualSDDP: MSLBO

function build(A::Vector{Matrix{Float64}},
                B::Vector{Matrix{Float64}},
                T::Vector{Matrix{Float64}},
                c::Vector{Matrix{Float64}},
                d::Vector{Matrix{Float64}},
                Ux::Float64,
                Uy::Float64,
                lb::Float64,
                ub::Float64,
                Lip::Float64,
                prob::Float64)
  function A_builder(t::Int, i::Int)
    return A[i]
  end

  function B_builder(t::Int, i::Int)
    return B[i]
  end

  function T_builder(t::Int, i::Int)
    return T[i]
  end
  
  function c_builder(t::Int, i::Int)
    return c[i]
  end

  function d_builder(t::Int, i::Int)
    return d[i]
  end

  function Ux_builder(t::Int)
    return Ux
  end

  function Uy_builder(t::Int)
    return Uy
  end

  function lb_builder(t::Int, nstages::Int)
    return lb
  end

  function ub_builder(t::Int, nstages::Int)
    return ub
  end

  function Lip_builder(t::Int, nstages::Int)
    return Lip
  end

  function prob_builder(t::Int)
    return prob
  end

  return MSLBO(A_builder, B_builder, T_builder, c_builder, d_builder,
                Ux_builder, Uy_builder, lb_builder, ub_builder,
                Lip_builder, prob_builder)
end

A = [reshape([1.0, 0], 2, 1)]
B = [reshape([-1.0, 0], 2, 1)]
T = [reshape([1.0, 1, 0, 0, 0, 1, 0, 1, 1, 1], 2, 5)]
c = [reshape([0.0, 1.0, 5.0, 10.0, 50.0], 1, 5)]
d = [reshape([0.0, 75.0], 1, 2)]
Ux = 100.0
Uy = 60.0
lb = 0.0
ub = 75*50.0
Lip = 50.0
prob = 1.0

M = build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob)

end
