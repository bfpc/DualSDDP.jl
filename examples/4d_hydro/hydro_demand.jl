include("demand.jl")

demand = Vector{Float64}[]
for i in 1:size(ENERGY_DEMAND,1)
  if i != ENERGY_DEMAND[i,1]
    println("Warning, demand for stage $(i) not defined.")
  end
  push!(demand, ENERGY_DEMAND[i,2:end])
end
