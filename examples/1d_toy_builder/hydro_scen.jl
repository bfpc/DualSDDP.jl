# Uncertainty
import Random

Random.seed!(2)
inflows = 40 .+ 20*randn(Main.nscen)
inflows = max.(inflows, 0)
