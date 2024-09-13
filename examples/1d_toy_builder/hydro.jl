module Hydro1d

include("hydro_scen.jl")

# inivol
inivol = [83.222]
nscen = Main.nscen
nstages = Main.nstages

using DualSDDP: MSLBO
using DualSDDP: build, SimpleLBO

A = reshape([1.0, 0], 2, 1)
B = reshape([-1.0, 0], 2, 1)
T = [1.0 1 0 0 0; 1 0 1 1 1]

c = [0., 1, 5, 10, 50]

d = [[[0, 75.0]], [[i, 75.0] for i in inflows], [[i, 75.0] for i in inflows], [[i, 75.0] for i in inflows]]

Ux = [100.0]
Uy = [60.0, 200, 15, 15, 75]

lb = fill(0.0, 4)
ub = 75*50*[4.0, 3, 2, 1]
Lip = 50*[4.0, 3, 2, 1]
prob = [[1.0], ones(nscen)/nscen, ones(nscen)/nscen, ones(nscen)/nscen]

stages = Vector{SimpleLBO}(undef, 4)
for t in 1:4
  stages[t] = SimpleLBO(A, B, T, c, d[t], Ux, Uy, lb[t], ub[t], Lip[t], prob[t])
end

M = build(stages)

end
