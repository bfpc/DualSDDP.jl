include(src * "structs.jl")

module Hydro1d

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

function d(t::Int, i::Int)
  if t == 1
    return [0, 75.0]
  else
    return [Main.inflows[i], 75.0]
  end
end

function Ux(t::Int)
  return [100]
end
function Uy(t::Int)
  return [60, 200, 15, 15, 75]
end

function lb(t::Int, nstages::Int)
  return 0
end
function ub(t::Int, nstages::Int)
  return (nstages - t + 1)*75*50
end

function prob(t::Int)
  if t == 1
    return [1.0]
  else
    return ones(Main.nscen)/Main.nscen
  end
end

M = Main.MSLBO(A,B,T,c,d,Ux,Uy,lb,ub,prob)

end
