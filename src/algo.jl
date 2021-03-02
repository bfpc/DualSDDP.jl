import JuMP

function choose(prob; norm=1.0)
  v = rand()*norm
  acc = 0
  for (j,p) in enumerate(prob)
    acc += p
    if acc > v
      return j
    end
  end
end

function forward(stages, state0; debug=0)
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0)
    JuMP.optimize!(stage)

    j = choose(stage.ext[:prob])
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0")
      if debug > 1
        println("                         local decision: ", JuMP.value.(stage[:y]))
      end
    end
  end
end

function forward_dual(stages, state0; debug=0)
  ϵ = 1e-1

  gamma0 = 1.0
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0, gamma0)
    JuMP.optimize!(stage)

    # Regularize probabilities, so that there's a chance to sample every scenario
    gammas = JuMP.value.(stage.ext[:vars][2]) .+ ϵ*gamma0
    gammas_r = gammas .+ ϵ*gamma0
    if debug > 0
      println("  ", gammas)
    end
    j = choose(gammas_r, norm=sum(gammas_r))
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    gamma0 = gammas[j]
    if debug > 0
      println("Going out from stage $i, branch $j, state $state0, prob $gamma0")
    end
  end
end
