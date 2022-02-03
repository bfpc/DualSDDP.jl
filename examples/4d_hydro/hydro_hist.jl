include(src * "structs.jl")

module Hydro_Hist

using LinearAlgebra: diagm

include("hydro_data.jl")
include("hydro_demand.jl")
include("hist_scen.jl")

# array with vector for maximum inflow
Maxinflows = reshape(maximum(eafs, dims=2), (4,12))


Id = diagm(ones(n))
I5 = diagm(n_markets, n, ones(n))
Z  = zeros(n,n)
Z5 = zeros(n_markets,n)

function A(t::Int, i::Int)
  return [Id; Z5]
end

function B(t::Int, i::Int)
  return [-Id; Z5]
end

function T(t::Int, i::Int)
  return [Id Id zeros(n,n_thermal) zeros(n,n_def) zeros(n,n_xch);
          I5 Z5 M_thermal          M_deficit      M_xch]
end

function c(t::Int, i::Int)
  return [zeros(n); spill_cost*ones(n);
          thermal_c; deficit_c; xch_cost]
end

function d(t::Int, i::Int)
  p = (t-1)%12 + 1
  return [eafs[1:n,i,p]; demand[t]]
end

function Ux(t::Int)
  return Maxvol
end
function Uy(t::Int)
  p = (t-1)%12 + 1
  max_deficit_t = [l*d for d in demand[t] for l in deficit_levels]
  return [Maxturb; Maxinflows[1:n,p]; GTmax .- GTmin;
          max_deficit_t; max_xch]
end

function lb(t::Int, nstages::Int)
  return 0
end
function ub(t::Int, nstages::Int)
  tot = 0.
  for t in nstages:-1:t
    tot += sum(c(t,1) .* Uy(t))
  end
  return tot
end

function Lip(t::Int, nstages::Int)
  tot = 0.
  for t in nstages:-1:t
    tot += maximum(c(t,1))
  end
  return tot
end

function prob(t::Int)
  if t == 1
    return [1]
  else
    n = size(eafs,2)
    return ones(n)/n
  end
end

M = Main.MSLBO(A,B,T,c,d,Ux,Uy,lb,ub,Lip,prob)

end
