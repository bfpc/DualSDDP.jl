# DualSDDP

DualSDDP implements the Dual SDDP algorithm
for multistage stochastic linear problems with stagewise independent uncertainty,
allowing for polyhedral risk measures such as Expectation-AV@R.
The theory is developped in https://arxiv.org/abs/2107.10930.

It also implements, for comparison:
- A primal SDDP algorithm
  * with outer approximations of the value function, and
  * stochastic sampling in the scenario tree.
- a deterministic DDP algorithm
  * with inner and outer approximations for the value function, and
  * worst-gap selection in the scenario tree.

# Usage

- Create a module `M`, of type `MSLBO`, containing the problem data;
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
or simply issue `pkg> rm Gurobi`.
