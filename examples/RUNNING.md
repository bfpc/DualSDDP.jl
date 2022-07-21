Examples
========
* 1d model

  This is a toy model, with 1 state variable, 4 stages and 10 branches/stage
  Test cases:
  - config_1d_tiny.json
    Varies \alpha and \beta
    Useful for sanity check of installation.

  - config_1d.json
    Varies several \alpha and \beta
    Also used as base case for 1d_epsilon.

  - config_1d_epsilon.json
    Varies \alpha, \beta and the probability regularization \epsilon
    to explore the effects of different regularization in the forward sampling.
    Only for dual algorithm.

  - config_1d_nscen.json
    Varies \alpha, \beta and the branching number at each stage,
    to analyse its impact on the convergence and running times.

  - config_1d_lip.json
    Varies \alpha, \beta and the Lipschitz factor
    to analyse its impact on the convergence and running times.

* 2d model
  This is a small-scale problem, with 2 state variables
  and a modest number of local decision variables.
  This complements the tests varying \alpha, \beta
  and algorithm parameters (Lipschitz, epsilon)

  - config_2d.json
    Base setting: 82 branches per stage, 12 stges.

  - config_2d_epsilon.json

  - config_2d_lip.json

* 4d model
  A medium-scale problem: 4 state variables
  and hundreds of local decision variables

  - config_4d.json
    Basic setting: 82 branches per stage, 12 stages
    Varies \alpha, \beta and the Lipschitz constant


Running
=======

For the first run, one must install the DualSDDP package
in the current Project, issuing

    $ julia
    julia> import Pkg
    julia> Pkg.activate(".")
    julia> Pkg.add("..")

Then, to run the examples in this directory, one can usually issue

    $ julia batch.jl <config.json>

for a given JSON config file.


Analysis
========

Several analysis can be made from the results we have.
Most are accessed through the graph_* scripts, so for example

    julia graph_lip.jl config_2d_lip.json

will output a summary table and plot the evolution of upper bounds for
the dual algorithm, varying the Lipschitz constant.


* graph_bounds.jl will plot bounds for all algorithms (primal, dual, inner-outer)
* graph_times.jl for timing per iteration
* graph_epsilon.jl for analysis of bounds with different regularizaion \epsilon
