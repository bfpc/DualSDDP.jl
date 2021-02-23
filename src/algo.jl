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

function forward(stages, state0)
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0)
    JuMP.optimize!(stage)

    j = choose(stage.ext[:prob])
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    println("Going out from stage $i, branch $j, state $state0")
  end
end

function forward_dual(stages, state0)
  gamma0 = 1.0
  for (i,stage) in enumerate(stages)
    set_initial_state!(stage, state0, gamma0)
    JuMP.optimize!(stage)

    gammas = JuMP.value.(stage.ext[:vars][2]) # + epsilon
    println("  ", gammas)
    j = choose(gammas, norm=sum(gammas))
    state0 = JuMP.value.(stage.ext[:vars][1][:,j])
    gamma0 = gammas[j]
    println("Going out from stage $i, branch $j, state $state0, prob $gamma0")
  end
end
