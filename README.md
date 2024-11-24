# DualSDDP

DualSDDP implements the Dual SDDP algorithm
for multistage stochastic linear problems with stagewise independent uncertainty,
allowing for polyhedral risk measures such as Expectation-AV@R.
The theory is developped in the paper
[Dual SDDP for risk-averse multistage stochastic programs](https://www.sciencedirect.com/science/article/abs/pii/S0167637723000603)
[[arXiv version](https://arxiv.org/abs/2107.10930)].

It also implements, for comparison:
- A primal SDDP algorithm
  * with outer approximations of the value function, and
  * stochastic sampling in the scenario tree;
- a deterministic DDP algorithm
  * with inner and outer approximations for the value function, and
  * worst-gap selection in the scenario tree;
- a Dynamic Programming upper bound
  * built from a set of states in each stage, in the spirit of [Philpott et al.](https://www.jstor.org/stable/23481808).

# Usage

- Create a module `M`, of type `MSLBO`, containing the problem data;
  * A simple example is `examples/1d_toy/hydro.jl`.
  * An interface for problems with only right-hand side uncertainty is provided through `src/build.jl`, and the analog of the above example can be found in `examples/1d_toy_builder/hydro.jl`.
- Determine the number of stages of the problem;
- Choose the E-AV@R parameters and create builder functions (primal or dual)
  using `mk_primal_avar` and `mk_copersp_avar` (respectively);
- Set the LP solver for subproblems;
- Set the initial value;
- Set the number of iterations of the algorithm;
- call `primalsolve`, `dualsolve` or `problem_child_solve`.

# Examples

The `examples` directory contains both "one shot" examples (in subdirectories),
as well as `Config.jl` batches and configs.

The first time you use the examples directory, you'll have to instantiate it.
The fastest way of doing so is using `Pkg` commands (in the `examples` directory):


```julia
pkg> activate .
pkg> dev ..
```

This will use the `Project.toml` we provide,
then the current checkout for the `DualSDDP` package,
and finally resolve all packages.

In so doing, it will also require that the `Gurobi` package is installed,
which was used for some of the examples.
If that is not desired, one can remove it from the `Packages.toml` file,
or simply issue `pkg> rm Gurobi` before `dev ..`.

## Running from the command line

The one-shot examples can be invoked, for example, by
`examples/1d_toy $ julia example.jl`.

The `Config.jl` examples were designed to run several cases in sequence,
likely on a remote server.
A sample command line invocation is
`examples $ julia batch.jl config_1d.json`

## Running from the REPL

For the one-shot examples, simply `cd` into their directory and include
the corresponding `example.jl` file:

```julia
shell> cd examples/1d_toy/
julia> include("example.jl")
```


For the `Config.jl` examples, set `ARGS[1]` to the desired configuration
file, and then run `include("batch.jl")`.
This will produce output files in `data/output/$(save_path)`, which can be
further processed using the `graph*.jl` scripts.
