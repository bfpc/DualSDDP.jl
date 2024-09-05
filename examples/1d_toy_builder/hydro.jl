module Hydro1d

include("hydro_scen.jl")
include("hydro_conf.jl")

# inivol
inivol = [83.222]

using DualSDDP: MSLBO
using ..Build: build, StageMLSBO

A = [reshape([1.0, 0], 2, 1), reshape([1.0, 0], 2, 1),
      reshape([1.0, 0], 2, 1), reshape([1.0, 0], 2, 1)]

B = [reshape([-1.0, 0], 2, 1), reshape([-1.0, 0], 2, 1),
      reshape([-1.0, 0], 2, 1), reshape([-1.0, 0], 2, 1)]

T = [reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5),
      reshape([1.0 1 0 0 0; 1 0 1 1 1], 2, 5)]

c = [[0.1, 1, 5, 10, 50], [0, 1, 5, 10, 50],
      [0, 1, 5, 10, 50], [0, 1, 5, 10, 50]]

d = fill(fill([0.0, 75.0], 10), 4)

Ux = [100.0]
Uy = [60.0, 200, 15, 15, 75]

lb = fill(0.0, 4)
ub = fill(75*50, 4)
Lip = fill(50.0, 4)
prob = [[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]

#M = build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, 4)

stages = Vector{StageMLSBO}(undef, 4)
for t in 1:4
  stages[t] = StageMLSBO(A[t], B[t], T[t], c[t], d[t], lb[t], ub[t], Lip[t], prob[t])
end

M = build(stages, Ux, Uy)

end