# Models

We consider a hydro-thermal generation planning model.

## Problem description

Each hydro reservoir has a certain stored volume $v_\text{ini}$,
and will get a random inflow $I(\omega)$ along the planning period.
The hydro plants have a maximum power generation $\bar{g^h}$,
and can further spill water (at unlimited capacity),
so that the reservoir stored volume at the end of the planning period is
$$v_{\text{end}} = v_\text{ini} + I(\omega) - g^h - \text{spill}$$
for each reservoir.

Given this dynamic, let's now expose the static problem
for a given time period.
Each thermal unit $k$ has a minimum and maximum power generation,
$\underline{g^t_k}$ and $\bar{gt_k}$,
and a unit cost $c_k$.
Combining hydro and thermal generation at each system $s$,
and allowing for power exchanges between systems,
we must supply for an energy demand of $D_s$.

Unattended demand, also known as load shedding,
incurs a unit cost of $M$.
A more realistic load shedding model includes several _levels_
of load shedding, with increasing unit costs per level.

## Optimization model

We have



$$\begin{align*}
\min \qquad  & \mathbb{E}\left[ \sum_{\tau=1}^T c \cdot g^t_\tau + M \text{def}_\tau \right] \\
\text{s.t.} \qquad & \underline{g^t} \leq g^t_\tau \leq \bar{g^t} \\
            & 0 \leq g^h_\tau \leq \bar{g^h} \\
            & 0 \leq \text{def}_\tau \\
            & v_{\text{end},\tau} = v_{\text{ini},\tau} + I_\tau(\omega) - g^h_\tau - \text{spill}_\tau \\
            & M_H g_{\tau}^h + M_T g^t_{\tau} + \text{def}_\tau + X^{ch}_{\tau} = D_\tau
\end{align*}$$

Here, most variables are vectors:
$c$ for unit costs for thermal generation,
$g^t_\tau$ for planned thermal generation at each plant at stage $\tau$,
$g^h$ for planned hydro generation, and so on.
The matrices $M_H$ and $M_T$ are incidence matrices
that indicate to which system a given plant is attached.

## The test cases

We have three models.
- A 1D toy model: one hydro reservoir, two thermal plants;
  as slack variables, it has spillage and load shedding.

- A 2D model: two interconnected systems, each one with a reservoir,
  the first has two thermal plants, the second has only one;
  each system has four levels of load shedding.

- A 4D model: 4 interconnected systems, with a transshipment node;
  one consolidated hydro plant per system, (39,15,39,33) thermal plants,
  four levels of load shedding.

The parameters are given in each of the dedicated directories,
`1d_toy`, `2d_hydro`, `4d_hydro`.

