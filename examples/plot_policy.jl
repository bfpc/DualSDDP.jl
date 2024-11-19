using PyPlot: plot

"""
    plot_tangents(stage, xrange)

Plots PrimalCuts from `stage` problem as lines, over `xrange`.

# Example
```julia
plot_tangents(stage, [0, 100.0])
```
"""
function plot_tangents(stage, xrange)
  cuts = stage.ext[:cuts]
  for cut in cuts
    a = cut.slope
    b = cut.value - cut.slope'*cut.x0
    plot(xrange, a .* xrange .+ b, "C1")
  end
end

"""
    plot_dual_as_vertices(stage)

Plots DualCuts from `stage` problem as vertices for the upper bound, over `xrange`.

# Example
```julia
plot_dual_as_vertices(stage)
```
"""
function plot_dual_as_vertices(stage)
  dualcuts = stage.ext[:cuts]
  vertices = []
  for cut in dualcuts
    x = cut.slope_π
    y = - cut.slope_γ
    push!(vertices, [x,y])
  end
  plot(first.(vertices), last.(vertices), "o")
end

"""
    plot_policy(stage)

Plots both inner and outer policy from `stage` problem.
"""
function plot_policy(stage::DualSDDP.IO_stage)
  vertices = []
  for vertex in stage.inner_vertices
    push!(vertices, [vertex.point[1], vertex.value])
  end
  plot(first.(vertices), last.(vertices), "o")
  
  xmin = minimum(first.(vertices))
  xmax = maximum(first.(vertices))
  @show xmin, xmax
  for cut in stage.cuts
    a = cut.slope
    b = cut.intercept
    plot([xmin, xmax], a .* [xmin, xmax] .+ b, "C1")
  end
end
