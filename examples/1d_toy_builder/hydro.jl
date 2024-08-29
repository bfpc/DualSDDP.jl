module Hydro1d

include("hydro_scen.jl")

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

d = [0, 75.0]

Ux = [[100.0], [100.0], [100.0], [100.0]]
Uy = [[60.0, 200, 15, 15, 75], [60.0, 200, 15, 15, 75],
      [60.0, 200, 15, 15, 75], [60.0, 200, 15, 15, 75]]

lb = 0.0
ub = 75*50.0
Lip = 50.0
prob = [[0.5, 0.5], [0.5, 0.5], [0.5, 0.5], [0.5, 0.5]]

#M = build(A, B, T, c, d, Ux, Uy, lb, ub, Lip, prob, 4)

stages = Vector{StageMLSBO}(undef, 4)
for t in 1:4
  stages[t] = StageMLSBO(A[t], B[t], T[t], c[t], Ux[t], Uy[t], prob[t])
end

M = build(stages, d, lb, ub, Lip)

end