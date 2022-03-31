module Evil10d

n = Main.dim
import Random
Random.seed!(42)
Cost_matrix = rand(n,100)

using LinearAlgebra
using DualSDDP: MSLBO

function A(t::Int, i::Int)
  return diagm(ones(n))
end

function B(t::Int, i::Int)
  return -diagm(ones(n))
end

function T(t::Int, i::Int)
  return diagm(ones(n))
end

function c(t::Int, i::Int)
  return Cost_matrix[:,t]
end

function d(t::Int, i::Int)
  if t == 1
    return zeros(n)
  else
    return Main.ξ[:,i]
  end
end

function Ux(t::Int)
  return 5*ones(n)
end
function Uy(t::Int)
  return 2*ones(n)
end

function lb(t::Int, nstages::Int)
  return 0
end
function ub(t::Int, nstages::Int)
  return (nstages - t + 1)*2*n
end

function Lip(t::Int, nstages::Int)
  return (nstages - t + 1)*Main.lip_factor
end

function prob(t::Int)
  if t == 1
    return [1.0]
  else
    p = Main.p
    v = [p^i for i in 1:Main.nscen]
    return v./sum(v)
  end
end

M = MSLBO(A,B,T,c,d,Ux,Uy,lb,ub,Lip,prob)

end
