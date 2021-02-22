include("structs.jl")

module Hydro1d

nscen = 10

function A(t::Int, i::Int)
  return reshape([1.0, 0], 2, 1)
end

function B(t::Int, i::Int)
  return reshape([-1.0, 0], 2, 1)
end

function T(t::Int, i::Int)
  return [1.0 1 0 0 0; 1 0 1 1 1]
end

function c(t::Int, i::Int)
  return [0, 1, 5, 10, 50]
end

import Random
Random.seed!(2)
af = 40 .+ 20*randn(nscen)
function d(t::Int, i::Int)
  return [af[i], 75]
end

function Ux(t::Int)
  return [100]
end
function Uy(t::Int)
  return [60, 200, 15, 15, 75]
end

function lb(t::Int)
  return 0
end

function prob(t::Int)
  return ones(nscen)/nscen
end

M = Main.MSLBO(A,B,T,c,d,Ux,Uy,lb,prob)

end
