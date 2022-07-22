# Algorithm Description

We give further details on the dual dynamic programming algorithm
we use to solve the (regularized) dual problem in the recursion (13)
of https://arxiv.org/abs/2107.10930.

- The problem is **setup** with an approximation to the value functions
  which are guaranteed upper bounds for that function.

- The **forward pass** starts from $t = 1$ and solves equation (13)
  for increasing values of $t$.
  The first stage yields a starting value for $\pi$, the dual state,
  and fixes $\gamma$ (node probability) to one.
  For the following stages, it starts from the current pair $(\pi,\gamma)$
  and solves equation (13).
  It yields all solutions $(\pi_j,\gamma_j)$ for outgoing states at each scenario.

- The selection of the **next state** is weighted by the probabilities $\gamma_j$,
  which correspond to the (current) change-of-measure,
  and biases the selection towards states that most contribute to the
  value function.
  To ensure that all branches are possible,
  we actually use $\gamma_j + \epsilon$ as weights.

- If desired (and as default for the main loop),
  the probabilities $\gamma$ are **normalized** at each stage to one, or zero.
  This usually improves numerical stability.

- Then the algorithm **adds cuts** to the value functions for all stage,
  from the first to the last one.
  Since the forward step must calculate solutions to all scenarios,
  we chose not to add cuts in a ``backwards'' fashion,
  since the solution of the forward pass already allows to construct
  a valid cut from its value and optimal mutiplier.

- The algorithm performs the ``double sweep'' of forward solves
  and cut addition for a given number of iterations.
