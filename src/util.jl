function opt_recover(model::JuMP.Model, prefix::String, errmsg::String;
                     extra_print::Function=()->nothing)
    JuMP.optimize!(model)
    status = JuMP.primal_status(model)
    if status != JuMP.FEASIBLE_POINT
      # trying to recover
      JuMP.MOI.Utilities.reset_optimizer(model)
      JuMP.optimize!(model)
      status = JuMP.primal_status(model)

      if status != JuMP.FEASIBLE_POINT
        JuMP.write_to_file(model, "$(prefix).mps")
        JuMP.write_to_file(model, "$(prefix).lp")
        extra_print()
        error(errmsg * "\nSolver primal status = $(status)")
      end
    end

end

import JSON

"""
    write_cuts_to_file(model::MSSP, filename::String)

Write the cuts that form the policy in `model` to `filename` in JSON format,
similar to `SDDP.write_cuts_to_file`.
"""
function write_cuts_to_file(model::MSSP, filename::String)
    cuts = Dict{String,Any}[]
    for (i, stage) in enumerate(model)
        node_cuts = Dict(
            "node" => string(i),
            "single_cuts" => Dict{String,Any}[],
        )
        for cut in stage.ext[:cuts]
          push!(node_cuts["single_cuts"], to_dict(cut))
        end
        push!(cuts, node_cuts)
    end
    open(filename, "w") do io
        return write(io, JSON.json(cuts))
    end
    return
end

function to_dict(cut::PrimalCut)
  Dict(
       "intercept" => cut.value - cut.slope' * cut.x0,
       "coefficients" => Dict(["x$(n)" => slope_n for (n, slope_n) in enumerate(cut.slope)]),
       "state" => Dict(["x$(n)" => state_n for (n, state_n) in enumerate(cut.x0)]),
      )
end

function to_dict(cut::DualCut; coperspective=false)
  intercept = (coperspective ? 0 : cut.slope_γ)
  d = Dict(
       "intercept" => intercept,
       "coefficients" => Dict(["x$(n)" => slope_n for (n, slope_n) in enumerate(cut.slope_π)]),
       "state" => Dict(["x$(n)" => state_n for (n, state_n) in enumerate(cut.π0)]),
      )
  if coperspective
    d["coefficients"]["γ"] = cut.slope_γ
    d["state"]["γ"] = cut.γ0
  end
  return d
end
