# FEProblemSolve

The `FEProblemSolve` class has two main roles:

- handle a variety of parameters for the [linear and nonlinear solves](systems/NonlinearSystem.md)
- encapsulate the solve function call with a [geometric grid sequencing calls](syntax/Executioner/index.md#Grid Sequencing) for nonlinear problems
- encapsulate each nonlinear system within a multi-system fixed point loop. This loop is nested within the geometric grid sequencing.

The `FEProblemSolve` is a solve executioner nested inside most executioners,
such as [Steady](executioners/Steady.md) and [Transient](executioners/Transient.md) but notably *not* in the [Eigenvalue](executioners/Eigenvalue.md) executioner.

## Multi-system solve capabilities id=multi-system_solve

If using multiple nonlinear systems, the default behavior of the `FEProblemSolve` will be to solve them one by one,
in the order that they were specified, without iterating between systems.

If the [!param](/Executioner/Steady/multi_system_fixed_point) parameter is set to true, this solve will be iterated.
The user must pass a convergence object to the [!param](/Executioner/Steady/multi_system_fixed_point_convergence)
to let the `FEProblemSolve` know when to terminate the fixed point loop.

The fixed point iterations can be accelerated using the
[!param](/Executioner/Steady/multi_system_fixed_point_algorithm) parameter:

- `relaxation` (default) applies simple under/over-relaxation to each system solution, blending the
  new solution with the previous iterate.
- `secant` applies a projected (Aitken) secant acceleration to the coupled system solutions as a
  whole, once per fixed point sweep. A single step size is obtained from the inner products of the
  coupled fixed point residual change, `mu = (dx.dg) / (dg.dg)`, which reduces to the scalar secant
  method in 1D. Applying the acceleration to the whole coupled state (rather than to each system
  independently) is what captures the coupling between systems; a per-system secant would model
  each system as a self-map and is not robust. The first sweep, having no history to build the
  secant slope, takes a Picard step.

In both cases the [!param](/Executioner/Steady/multi_system_fixed_point_relaxation_factor) is used to
under/over-relax the update. A value of 1 disables relaxation, so the `secant` algorithm with a
factor of 1 performs an un-relaxed secant update.

!alert note
Options are currently limited for setting a multi-system fixed point convergence. We do not recommend using the
nonlinear residual with a [VariableResidual.md] postprocessor or a [DefaultNonlinearConvergence.md] as these
are not re-computed the end of a multi-system fixed point iteration.
