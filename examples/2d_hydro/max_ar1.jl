function max_us(Phi0, Phi1, noises, u0::Vector, nper::Int)
  per = length(Phi0)
  @assert per == length(Phi1)

  list_max = []

  for i = 1:nper
    t = ((i-1)%per) + 1
    u = zero(u0)
    forecast = Phi0[t] .+ Phi1[t]*u0
    for noise in noises[i]
      u = max.(u, noise*forecast)
    end
    u0 = u
    push!(list_max, u0)
  end

  return list_max
end
