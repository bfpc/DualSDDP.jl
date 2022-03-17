module DualSDDP

import JuMP
using JuMP: Model, @variable, @constraint, @objective,
            set_lower_bound, set_upper_bound, fix

include("structs.jl")
include("problem.jl")
include("risk_models.jl")
include("algo.jl")
include("ub.jl")
include("util.jl")

include("analysis.jl")

export MSLBO
export mk_primal_avar, mk_copersp_avar
export primalsolve, dualsolve

export pvf_info, dvf_info

end
