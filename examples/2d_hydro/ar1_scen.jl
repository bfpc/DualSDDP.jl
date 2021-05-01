using LinearAlgebra: diagm

include("../readblock.jl")

file = open("scenarios_2REE_small.txt", "r")
idxs, data = readblock(file)

# scen_prob[t] is a vector with scenario probabilities
# noise[t][i] is a (diagonal) (n x n) matrix
function organize_arrays(idxs, data)
scen_prob = Vector{Float64}[]
noise     = []

ps_t     = Float64[]
noises_t = Matrix{Float64}[]
t        = 1
for (i, d) in zip(idxs, data)
  if i > t
    push!(scen_prob, copy(ps_t))
    push!(noise,     copy(noises_t))
    empty!(ps_t)
    empty!(noises_t)
    t = t+1
  end
  push!(ps_t, d[1])
  push!(noises_t, diagm(d[2:end]))
end
push!(scen_prob, copy(ps_t))
push!(noise,     copy(noises_t))
empty!(ps_t)
empty!(noises_t)

return scen_prob, noise
end

scen_prob, noise = organize_arrays(idxs, data)
for ps in scen_prob
  @assert sum(ps) ≈ 1.
end
