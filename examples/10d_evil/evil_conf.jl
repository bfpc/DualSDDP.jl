# Parameters
nstages = 20
alpha   = 0.4
niters  = 100
ub_step = 10

p = 0.3

using LinearAlgebra: diagm
nscen = 10

dim = 20 
ξ = diagm(ones(dim))

# inivol
inivol = 5.0*ones(dim)
