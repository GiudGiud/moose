# =====================================================================
# Gradient-free optimization search of
# one rate constant (k1 of A -> B) from the observed final [A], [B] and [C]
# in the network A -> B -> C (forward_optimization.i).
#
# The rate constant k2 of B -> C is held at its reference value.
#
# The searched parameters are log10(k), so the rate constants stay positive and
# the search spans decades. The objective is half the sum of squared misfits
# against synthetic observations computed with k1 = 2 /s and k2 = 0.5 /s, i.e.
# log10(k1) = 0.301, log10(k2) = -0.301. It is minimized with the TAO Nelder-Mead
# simplex method, which needs no gradient and therefore no adjoint problem.
# =====================================================================

[Optimization]
[]

[OptimizationReporter]
  type = GeneralOptimization
  objective_name = objective_value
  parameter_names = 'logk1'
  num_values = '1'
  initial_condition = '0'
[]

[MultiApps]
  [forward]
    type = FullSolveMultiApp
    input_files = forward_optimization.i
    cli_args = "objective=0.5*((A-0.135339805736)^2+(B-0.628266304216)^2+(C-0.236393890048)^2);Outputs/console=false"
    execute_on = FORWARD
  []
[]

[Transfers]
  [to_forward]
    type = MultiAppReporterTransfer
    to_multi_app = forward
    from_reporters = 'OptimizationReporter/logk1'
    to_reporters = 'params/logk1'
  []
  [from_forward]
    type = MultiAppReporterTransfer
    from_multi_app = forward
    from_reporters = 'objective/value'
    to_reporters = 'OptimizationReporter/objective_value'
  []
[]

[Executioner]
  type = Optimize
  tao_solver = taonm
  # For Nelder-Mead the "gradient" norm is a measure of the simplex size; the
  # default tolerance (1e-8) stops ~1e-3 away from the minimum in log10(k)
  petsc_options_iname = '-tao_gatol'
  petsc_options_value = '1e-16'
  verbose = true
[]

[Outputs]
  csv = true
[]
