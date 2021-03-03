# Parameters
nstages = 4
beta    = 0.4
niters  = 40

# Uncertainty
import Random

nscen = 10
Random.seed!(2)
inflows = 40 .+ 20*randn(nscen)

# inivol
inivol = 83.222
