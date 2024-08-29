# Helper types
State = Vector{Float64}
DualState = Tuple{Float64, Vector{Float64}}

# Multi-Stage Linear Bellman Operator
#   A,B,T,c,d are functions of (time, branch)
#   Ux, Uy are functions of time
#   lb, ub, Lip are functions of (time, Horizon)
#   prob is a function of time
struct MSLBO
  A :: Function # Ax_t + Bx_{t-1} + Ty = d
  B :: Function
  T :: Function
  c :: Function # marginal cost of y_t
  d :: Function
  Ux :: Function # upper bound on the positive state x
  Uy :: Function # upper bound on the positive control y
  lb :: Function # lower bound at each stage
  ub :: Function # upper bound on the value of the problem
  Lip :: Function # upper bound on the Lipschitz constant
  prob :: Function # reference probability over branches
end

struct PrimalCut
  value :: Float64
  slope :: Vector{Float64}
  x0    :: Vector{Float64}
end

struct DualCut
  value   :: Float64
  slope_π :: Vector{Float64}
  slope_γ :: Float64
  π0      :: Vector{Float64}
  γ0      :: Float64
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
"""
struct StageMLSBO
    A::Matrix{Float64}
    B::Matrix{Float64}
    T::Matrix{Float64}
    c::Vector{Float64}
    Ux::Vector{Float64}
    Uy::Vector{Float64}
end

# Multi-stage stochastic problems, with some memory in .ext
struct MSSP <: AbstractArray{Model, 1}
  stages::Vector{Model}
  ext::Dict
end

MSSP(stages::Vector{Model}) = MSSP(stages, Dict())

Base.size(pb::MSSP) = Base.size(pb.stages)
Base.IndexStyle(::Type{<:MSSP}) = IndexLinear()
Base.getindex(pb::MSSP, i::Int) = pb.stages[i]
