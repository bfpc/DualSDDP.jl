import JuMP
using JuMP: @variable, @constraint

function mk_primal_avar(beta)
  function primal_avar(m, t, ps)
    n = length(t)
    # RU representation: extra variables, constraints and objective function
    z = @variable(m, _z_ru)
    u = @variable(m, _u_ru[i=1:n] >= 0)
    # JuMP.set_name(z, "_z_ru")
    # for (i,ui) in enumerate(u)
    #   JuMP.set_name(ui, "_u_ru[$i]")
    # end
    # JuMP.set_lower_bound.(u, 0.0)
    @constraint(m, gamma[i=1:n], z + u[i] >= t[i])
    JuMP.@objective(m, Min, z + 1/beta*sum(ps' * u))
  end
end

function mk_copersp_avar(beta)
  function copersp_avar(m, gamma, ps, gamma_in)
    @constraint(m, _z_ru, sum(gamma) == gamma_in)
    @constraint(m, gamma .<= (gamma_in/beta) .* ps)
  end
end

# Do we need sum() over ps' * VAR ?
