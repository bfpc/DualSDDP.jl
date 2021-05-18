include(src * "structs.jl")

module HydroAR1

using LinearAlgebra: diagm

include("hydro_data.jl")
include("hydro_demand.jl")
include("ar1_model.jl")
include("ar1_scen.jl")

# array with vector for maximum inflow
include("max_ar1.jl")
Maxinflows = max_us(Phi0, Phi1, noise, u0, nstages)


Id = diagm(ones(n))
I5 = diagm(n_markets, n, ones(n))
Z  = zeros(n,n)
Z5 = zeros(n_markets,n)

function A(t::Int, i::Int)
  return [Id -Id; Z Id; Z5 Z5]
end

function B(t::Int, i::Int)
  return [-Id Z; Z -Phi1[t]*noise[t][i]; Z5 Z5]
end

function T(t::Int, i::Int)
  return [Id Id zeros(n,n_thermal) zeros(n,n_def) zeros(n,n_xch);
          Z  Z  zeros(n,n_thermal) zeros(n,n_def) zeros(n,n_xch);
          I5 Z5 M_thermal          M_deficit      M_xch]
end

function c(t::Int, i::Int)
  return [zeros(n); spill_cost*ones(n);
          thermal_c; deficit_c; xch_cost]
end

function d(t::Int, i::Int)
  return [zeros(n); noise[t][i]*Phi0[t]; demand[t]]
end

function Ux(t::Int)
  return [Maxvol; Maxinflows[t]]
end
function Uy(t::Int)
  max_deficit_t = [l*d for d in demand[t] for l in deficit_levels]
  return [Maxturb; Maxinflows[t]; GTmax .- GTmin;
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

function prob(t::Int)
  return scen_prob[t]
end

M = Main.MSLBO(A,B,T,c,d,Ux,Uy,lb,ub,prob)

end
