module Evil1d

using DualSDDP: MSLBO

function A(t::Int, i::Int)
  return reshape([1.0],1, 1)
end

function B(t::Int, i::Int)
  return reshape([-1.0], 1, 1)
end

function T(t::Int, i::Int)
  return reshape([1.0],1, 1)
end

function c(t::Int, i::Int)
  return [1.0]
end

function d(t::Int, i::Int)
  if t == 1
    return [0.0]
  else
    return [Main.ξ[i]]
  end
end

function Ux(t::Int)
  return [10]
end
function Uy(t::Int)
  return [20]
end

function lb(t::Int, nstages::Int)
  return 0
end
function ub(t::Int, nstages::Int)
  return (nstages - t + 1)*20
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
