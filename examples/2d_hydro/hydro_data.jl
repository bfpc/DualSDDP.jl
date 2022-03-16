# Raw julia vectors, to be reformatted
include("data_2REE.jl")

# Manually set parameters
nstages = 24
ini_shift = 0

# Problem dimensions
# n: state variable (stored energy, past inflows)
# n_markets
n = length(HYDRO_PLANTS)
n_markets = length(SYSTEMS_DATA)

# Maximum stored energy
Maxvol = [x[6] for x in HYDRO_PLANTS]
# Maximum outflow
Maxturb = [x[5] for x in HYDRO_PLANTS]
# initial stored energy
inivol = [x[7] for x in HYDRO_PLANTS]

# spill cost
spill_cost = 0.001

# Thermal units: number, incidence, cost, bounds
n_thermal = length(THERMAL_PLANTS)

M_thermal = zeros(n_markets, n_thermal)
for i in 1:n_markets
  for (j,x) in enumerate(THERMAL_PLANTS)
    if x[2] == i
      M_thermal[i,j] = 1
    end
  end
end

thermal_c = [x[4] for x in THERMAL_PLANTS]

GTmax = [x[6] for x in THERMAL_PLANTS]
GTmin = [x[5] for x in THERMAL_PLANTS]

# Deficit (non-provided energy): number, incidence, cost, bounds
n_def_pat = length(ENERGY_DEFICITS_DATA)
n_def = n_def_pat * n_markets

M_deficit = kron(diagm(ones(n_markets)), ones(1,n_def_pat))

deficit_c = [x[3] for x in ENERGY_DEFICITS_DATA]
deficit_c = kron(ones(n_markets), deficit_c)

deficit_levels = [x[2] for x in ENERGY_DEFICITS_DATA]

# Exchange: number, incidence, cost, bounds
idx_xch = Tuple{Int,Int}[]
for i in 1:n_markets
  for j in 1:n_markets
    if i != j && ENERGY_EXCHANGES_CAPACITY_BETWEEN_SYSTEMS[i][j+1] > 0
      push!(idx_xch, (i,j))
    end
  end
end

n_xch     = length(idx_xch)

M_xch     = zeros(n_markets, n_xch)
for (k, (i,j)) in enumerate(idx_xch)
  M_xch[i,k] = -1
  M_xch[j,k] =  1
end

xch_cost  = [ENERGY_EXCHANGES_PENALTIES[i][j+1] for (i,j) in idx_xch]

max_xch   = [ENERGY_EXCHANGES_CAPACITY_BETWEEN_SYSTEMS[i][j+1] for (i,j) in idx_xch]

