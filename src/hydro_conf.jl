# Parameters
nstages = 2
beta    = 1.0

# Uncertainty
import Random

nscen = 3
Random.seed!(2)
inflows = 40 .+ 20*randn(nscen)

# inivol
inivol = 83.222
