import JuMP
using JuMP: Model, @variable, @constraint, @objective
using LinearAlgebra: dot

""" Return a JuMP model with the primal value function
based on stored cuts"""
function primal_value_function(stage; maxcuts=nothing)
  nx = size(stage[:x],1)
  ub_x = JuMP.upper_bound.(stage[:x][:,1])
  lb = JuMP.lower_bound(stage[:_z][1])

  vf = Model()
  @variable(vf, x[1:nx])
  @variable(vf, z >= lb)

  @constraint(vf, x .>= 0)
  @constraint(vf, x .<= ub_x)

  cuts = stage.ext[:cuts]
  if maxcuts != nothing
    cuts = cuts[1:maxcuts]
  end
  for cut in cuts
    cst, multipliers, x0 = cut.value, cut.slope, cut.x0
    @constraint(vf, z >= cst + multipliers'*(x .- x0))
  end

  @objective(vf, Min, z)

  return vf
end

""" For plotting purpose in 2D """
function eval_vf(vf, x1s, x2s)
  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2)
  for i = 1:n1
    x = x1s[i]
    for j = 1:n2
      JuMP.fix.(vf[:x], [x, x2s[j]])
      JuMP.optimize!(vf)
      z[i,j] = JuMP.objective_value(vf)
    end
  end

  return z
end


function pvf_info(stage; maxcuts=nothing)
  lb = JuMP.lower_bound(stage[:_z][1])
  cuts = stage.ext[:cuts]
  if maxcuts != nothing
    cuts = cuts[1:maxcuts]
  end
  coefs = Tuple{Float64, Vector{Float64}}[]
  for cut in cuts
    cst, multipliers, x0 = cut.value, cut.slope, cut.x0
    push!(coefs, (cst - dot(multipliers, x0), multipliers))
  end
  return coefs, lb
end

function eval_pvf(vf_info, x1s, x2s)
  coefs, lb = vf_info

  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2)
  for i = 1:n1
    x1 = x1s[i]
    for j = 1:n2
      x2 = x2s[j]
      curmax = lb
      for k in 1:length(coefs)
        cst, multipliers = coefs[k]
        v = cst + multipliers[1]*x1 + multipliers[2]*x2
        curmax = max(curmax, v)
      end
      z[i,j] = curmax
    end
  end

  return z
end

function norm1_evolution(vf_info, x1s, x2s)
  coefs, lb = vf_info
  niters = length(coefs)

  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2, niters+1)
  for i = 1:n1
    x1 = x1s[i]
    for j = 1:n2
      x2 = x2s[j]
      curmax = lb
      for k in 1:length(coefs)
        z[i,j,k] = curmax
        cst, multipliers = coefs[k]
        v = cst + multipliers[1]*x1 + multipliers[2]*x2
        curmax = max(curmax, v)
      end
      z[i,j,end] = curmax
    end
  end

  return sum(z, dims=(1,2))[1,1,:]
end

function norm1_evolution(pb::Vector{Model}, x1s, x2s)
  zs = Vector{Float64}[]
  for stage in pb[1:end-1]
    vf_info = pvf_info(stage)
    push!(zs, norm1_evolution(vf_info, x1s, x2s))
  end
  return zs
end



function dual_value_function(stage, ub, Lip; tol=1e-6, maxcuts=nothing)
  nx = size(stage[:π],1)

  vf = Model()
  @variable(vf, π[1:nx])
  # Fix γ to one, by homogeneity, for a start
  # @variable(vf, γ >= 0)
  γ = 1.0
  @variable(vf, z >= -ub*γ)

  @constraint(vf, π .<=  Lip)
  @constraint(vf, π .>= -Lip)

  cuts = stage.ext[:cuts]
  if maxcuts != nothing
    cuts = cuts[1:maxcuts]
  end
  for cut in cuts
    cst, mul_π, mul_γ, π0, γ0 = cut.value, cut.slope_π, cut.slope_γ, cut.π0, cut.γ0
    b = cst - mul_π'*π0 - mul_γ*γ0
    if b > 0
      println("Cut with positive intercept $(b), truncating")
      b = 0.0
    end
    if b < -tol
      println("Cut with negative intercept $(b)")
    end
    @constraint(vf, z >= b + mul_π'*π + mul_γ * γ)
  end

  @objective(vf, Min, z)

  return vf
end

function eval_df(df, x1s, x2s)
  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2)
  for i = 1:n1
    x = x1s[i]
    for j = 1:n2
      JuMP.fix.(df[:π], [x, x2s[j]])
      JuMP.optimize!(df)
      z[i,j] = JuMP.objective_value(df)
    end
  end

  return z
end

# function dual_value_function(stage, ub, Lip; tol=1e-6, maxcuts=nothing)
function dvf_info(stage; maxcuts=nothing, tol=1e-6)
  cuts = stage.ext[:cuts]
  if maxcuts != nothing
    cuts = cuts[1:maxcuts]
  end
  # Constant, π multiplier, γ multiplier
  coefs = Tuple{Float64, Vector{Float64}, Float64}[]
  for cut in cuts
    cst, mul_π, mul_γ, π0, γ0 = cut.value, cut.slope_π, cut.slope_γ, cut.π0, cut.γ0
    b = cst - mul_π'*π0 - mul_γ*γ0
    if b > 0
      println("Cut with positive intercept $(b), truncating")
      b = 0.0
    end
    if b < -tol
      println("Cut with negative intercept $(b)")
    end
    push!(coefs, (b, mul_π, mul_γ))
  end
  return coefs
end

function eval_dvf(vf_info, x1s, x2s, ub; γ::Float64=1.0)
  coefs = vf_info

  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2)
  for i = 1:n1
    x1 = x1s[i]
    for j = 1:n2
      x2 = x2s[j]
      curmax = -ub
      for k in 1:length(coefs)
        b, mul_π, mul_γ = coefs[k]
        v = b + mul_π[1]*x1 + mul_π[2]*x2 + mul_γ * γ
        curmax = max(curmax, v)
      end
      z[i,j] = curmax
    end
  end

  return z
end

function norm1_dual_evolution(vf_info, x1s, x2s, ub; γ::Float64=1.0)
  coefs = vf_info

  n1 = length(x1s)
  n2 = length(x2s)
  z = zeros(n1, n2, length(coefs)+1)
  for i = 1:n1
    x1 = x1s[i]
    for j = 1:n2
      x2 = x2s[j]
      curmax = -ub
      for k in 1:length(coefs)
        z[i,j,k] = curmax
        b, mul_π, mul_γ = coefs[k]
        v = b + mul_π[1]*x1 + mul_π[2]*x2 + mul_γ * γ
        curmax = max(curmax, v)
      end
      z[i,j,end] = curmax
    end
  end

  return sum(z, dims=(1,2))[1,1,:]
end

function norm1_dual_evolution(pb::Vector{Model}, M, nsteps)
  zs = Vector{Float64}[]
  nstages = length(pb)
  for k = 1:nstages-1
    stage = pb[k]
    vf_info = dvf_info(stage)
    Lip = M.Lip(k,nstages)
    x1s = -Lip:Lip/nsteps:Lip
    x2s = -Lip:Lip/nsteps:Lip
    push!(zs, norm1_dual_evolution(vf_info, x1s, x2s, M.ub(k,nstages)))
  end
  return zs
end



function compare_vf(stage, x1s, x2s, n1, n2, solver; primal=true, ub=nothing)
  if primal
    vf_1 = pvf_info(stage, maxcuts=n1)
    vf_2 = pvf_info(stage, maxcuts=n2)
    z1 = eval_pvf(vf_1, x1s, x2s)
    z2 = eval_pvf(vf_2, x1s, x2s)
  else
    vf_1 = dvf_info(stage, maxcuts=n1)
    vf_2 = dvf_info(stage, maxcuts=n2)
    z1 = eval_dvf(vf_1, x1s, x2s, ub)
    z2 = eval_dvf(vf_2, x1s, x2s, ub)
  end

  fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10,5))
  box = [minimum(x1s), maximum(x1s), minimum(x2s), maximum(x2s)]
  i1 = ax1.imshow(z1', origin="lower", extent=box)
  plt.colorbar(i1, ax=ax1)
  ax1.set_title("Using $(n1) cuts")

  i2 = ax2.imshow(z2', origin="lower", extent=box)
  plt.colorbar(i2, ax=ax2)
  ax2.set_title("Using $(n2) cuts")

  i3 = ax3.imshow((z2.-z1)', origin="lower", extent=box)
  plt.colorbar(i3, ax=ax3)
  ax3.set_title("Difference")
end

function norm1(stage::Model, x1s, x2s; primal=true, ub=nothing)
  if primal
    vf = pvf_info(stage, maxcuts=n1)
  end
end
