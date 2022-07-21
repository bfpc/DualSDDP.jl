using Config: parse!, ConfigManager, total_combinations, load
using Printf: @printf

# Loader functions
function read_ab_all(cfg::ConfigManager)
  bounds_dict = Dict()
  N = total_combinations(cfg)

  for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]

    bounds_dict[(alpha,beta)] = (data["primal lb"], data["dual ub"], data["inner recursive bound"], data["inner recursive iters"], data["io lb"], data["io ub"])
  end
  return bounds_dict
end

function read_ablip_all(cfg::ConfigManager)
  bounds_dict = Dict()
  N = total_combinations(cfg)

  for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]
    lip   = cfg["parameters"]["Lip"]

    bounds_dict[(alpha,beta,lip)] = (data["primal lb"], data["dual ub"], data["inner recursive bound"], data["inner recursive iters"], data["io lb"], data["io ub"])
  end
  return bounds_dict
end

function read_abe_dualub(cfg::ConfigManager)
  bounds_dict = Dict()
  N = total_combinations(cfg)

  for idx=1:N
    parse!(cfg, idx)
    data = load(cfg)

    alpha = cfg["risk-aversion"]["alpha"]
    beta  = cfg["risk-aversion"]["beta"]
    epsilon = cfg["parameters"]["epsilon"]

    bounds_dict[(alpha, beta, epsilon)] = (data["dual ub"])
  end
  return bounds_dict
end

# Auxiliary functions
function sets_in_dict(bounds_dict)
  nparams = length(first(bounds_dict)[1])
  sets = [Set() for i=1:nparams]

  for k in keys(bounds_dict)
    for i in 1:nparams
      push!(sets[i], k[i])
    end
  end

  lists = [sort([s...]) for s in sets]
  return lists
end

function best_lb(bounds_dict, alpha, beta, params)
  lb = -Inf
  for p in params
    pd = bounds_dict[(alpha, beta, p)]
    lb = max(lb, pd[1][end], pd[5][end])
  end
  return lb
end

# Printing functions
function print_bounds(bounds_dict, param_name)
  alphas, betas, params = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      @printf("Experiment α = %5.3f β = %5.3f:\n", alpha, beta)
      for p in params
        @printf("    %s=%4d yields final bound at %.17f\n", param_name, p, bounds_dict[(alpha, beta, p)][2][end])
      end
    end
  end
end

function print_relgaps(bounds_dict, param_name)
  alphas, betas, params = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      lb = best_lb(bounds_dict, alpha, beta, params)
      @printf("Experiment α = %5.3f β = %5.3f:\n", alpha, beta)
      for p in params
        gap = bounds_dict[(alpha, beta, p)][2][end]/lb - 1
        @printf("    %s=%4d relative gap %5.2f%%\n", param_name, p, 100*gap)
      end
    end
  end
end

function table_relgaps(bounds_dict, param_name, ts)
  alphas, betas, params = sets_in_dict(bounds_dict)

  nalphas = length(alphas)
  nbetas = length(betas)
  @printf(" %11s │ %4s │    Relative gap at iteration\n", "", "")
  @printf(" %11s │ %4s │", "(α,β)", param_name)
  for t in ts
          @printf(" %8d │", t)
  end
  println()
  sepline = repeat("─", 13) * "┼" * repeat("─", 6) * repeat("┼" * repeat("─", 10), length(ts)) * "┤"
  println(sepline)
  for i in 1:nalphas
    alpha = alphas[i]
    for j in 1:nbetas
      beta = betas[j]
      lb = best_lb(bounds_dict, alpha, beta, params)
      for p in params
        if p != params[2]
          print("             │")
        else
          @printf("(%4.2f, %4.2f) │", alpha, beta)
        end
        @printf(" %4d │", p)
        for t in ts
          gap = bounds_dict[(alpha, beta, p)][2][t]/lb - 1
          @printf(" %8.2f │", 100*gap)
        end
        println()
      end
      println(sepline)
    end
  end
end

