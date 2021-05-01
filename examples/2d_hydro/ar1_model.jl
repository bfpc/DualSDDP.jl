include("../readblock.jl")

# AR1 data
# Phi = AR coefficients
# Phi0 is an array with 12 n-vectors
# Phi1 is an array with 12 (n x n) matrices
#
# u0  = Last inflows, to bootstrap AR1
Phi0 = Array{Vector{Float64}}(undef, 12)
Phi1 = Array{Matrix{Float64}}(undef, 12)
file = open("linear_models_2REE.txt")
for i = 1:12
  idx, data = readblock(file)
  Phi0[i] = data[1]
  Phi1[i] = hcat(data[2:end]...)
end
idx, data = readblock(file)
u0 = data[1]
close(file)

# TODO: Verify dimension of Phi0, Phi1 correspond to state dimensions
