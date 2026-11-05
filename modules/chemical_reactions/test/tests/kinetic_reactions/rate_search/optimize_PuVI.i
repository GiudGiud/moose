# =====================================================================
# Gradient-free optimization search of one rate constant of the plutonium
# radiolysis network (forward_PuVI.i), that of
#   PuO22p + H -> PuO2p + Hp
# from the observed final log10([PuO2p]).
#
# The searched parameter is log10(k), starting one decade below the reference
# value k = 2e8 dm^3 mol^-1 s^-1 (log10(k) = 8.301) used to compute the
# synthetic observation. The objective, half the squared misfit of
# log10([PuO2p]), is minimized with the TAO Nelder-Mead simplex method, which
# needs no gradient and therefore no adjoint problem.
# =====================================================================

[Optimization]
[]

[OptimizationReporter]
  type = GeneralOptimization
  objective_name = objective_value
  parameter_names = 'logk'
  num_values = '1'
  initial_condition = '7.3'
[]

[MultiApps]
  [forward]
    type = FullSolveMultiApp
    input_files = forward_PuVI.i
    execute_on = FORWARD
  []
[]

[Transfers]
  [to_forward]
    type = MultiAppReporterTransfer
    to_multi_app = forward
    from_reporters = 'OptimizationReporter/logk'
    to_reporters = 'params/logk'
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
  # For Nelder-Mead the "gradient" norm is a measure of the simplex size. The
  # search stops ~1e-5 from the reference log10(k), where the objective is still
  # ~1e-11: the forward solve tolerances limit the attainable accuracy
  petsc_options_iname = '-tao_gatol'
  petsc_options_value = '1e-12'
  verbose = true
[]

[Outputs]
  csv = true
  print_nonlinear_converged_reason = false
  print_linear_converged_reason = false
[]
