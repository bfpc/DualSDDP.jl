import JuMP, SDDP
using JuMP: @variable, @constraint
using SDDP: State, @stageobjective

include("hydro_conf.jl")

function build_hydro(sp, t)
  @variable(sp, 0 <= vol <= 100, State, initial_value=inivol)
  @variable(sp, 0 <= gh <= 60)
  @variable(sp, 0 <= spill <= 200)
  @variable(sp, 0 <= gt[i=1:2] <= 15)
  @variable(sp, 0 <= def <= 75)
  @variable(sp, inflow)

  @constraint(sp, sum(gt) + gh + def == 75)
  @constraint(sp, vol.in + inflow - gh - spill == vol.out)

  if t == 1
      JuMP.fix(inflow, 0.0)
  else
    SDDP.parameterize(sp, inflows) do observed_inflow
      JuMP.fix(inflow, observed_inflow)
    end
  end

  @stageobjective(sp, spill + 5*gt[1] + 10*gt[2] + 50*def)
end

import GLPK
solver = GLPK.Optimizer;

pb = SDDP.LinearPolicyGraph(build_hydro;
                            stages = nstages,
                            sense  = :Min,
                            optimizer = solver,
                            lower_bound = 0.0);

SDDP.train(pb, iteration_limit=niters, risk_measure=SDDP.AVaR(beta))

using Statistics: mean, std

simulations = SDDP.simulate(pb, 500,
                           [:vol, :gh, :spill, :gt, :def])

objective_values = [
    sum(stage[:stage_objective] for stage in sim) for sim in simulations
]

μ = round(mean(objective_values), digits = 2)

ci = round(1.96 * std(objective_values) / sqrt(500), digits = 2)

println("Confidence interval: ", μ, " ± ", ci)
println("        Lower bound: ", round(SDDP.calculate_bound(pb), digits = 2))
