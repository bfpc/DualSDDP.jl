module DualSDDP

import JuMP
using JuMP: Model, @variable, @constraint, @objective,
            set_lower_bound, set_upper_bound, fix, set_normalized_coefficient

include("structs.jl")
include("problem.jl")
include("risk_models.jl")
include("algo.jl")
include("ub.jl")
include("util.jl")
include("Problem-Child.jl")

include("analysis.jl")

export MSLBO
export mk_primal_avar, mk_copersp_avar, mk_primal_io
export primalsolve, dualsolve, problem_child_solve, primalub

export pvf_info, dvf_info

end
