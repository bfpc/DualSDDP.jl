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
