# =====================================================================
# Forward model for the optimization-module rate searches: the consecutive
# first-order network A -> B -> C in a single 0-D element, run to t = 1 s, in
# two experiments solved side by side:
#   experiment 0 (species A, B, C)    starts from [A] = 1 M
#   experiment 1 (species A1, B1, C1) starts from [B] = 1 M
#
# The optimizer sets log10 of the searched rate constants in 'params/logk1'
# and/or 'params/logk2'. Rate constants that are not searched keep their
# reference values k1 = 2 /s and k2 = 0.5 /s, also used to compute the
# observations. The objective, set by the main application with cli_args, is
# half the sum of squared misfits of the observed final concentrations.
# =====================================================================

objective = '0.5*(A-0.135339805736)^2'

[Mesh]
  # Single element (0-D homogeneous kinetics)
  [single]
    type = ElementGenerator
    nodal_positions = '0 0 0
                       1 0 0'
    element_connectivity = '0 1'
    elem_type = EDGE2
  []
[]

[Functions]
  [k1]
    type = ParsedOptimizationFunction
    expression = '10^logk1'
    param_symbol_names = 'logk1'
    param_vector_name = 'params/logk1'
  []
  [k2]
    type = ParsedOptimizationFunction
    expression = '10^logk2'
    param_symbol_names = 'logk2'
    param_vector_name = 'params/logk2'
  []
[]

[Physics]
  [KineticReactions]
    [experiment0]
      species = 'A B C'
      diffusivities = '0 0 0'
      initial_conditions = '1 0 0'
      reactions = 'A -> B [k=k1]
                   B -> C [k=k2]'
    []
    [experiment1]
      species = 'A1 B1 C1'
      diffusivities = '0 0 0'
      initial_conditions = '0 1 0'
      reactions = 'A1 -> B1 [k=k1]
                   B1 -> C1 [k=k2]'
    []
  []
[]

[Reporters]
  [params]
    type = ConstantReporter
    real_vector_names = 'logk1 logk2'
    # Reference values log10(2) and log10(0.5)
    real_vector_values = '0.301029995663981; -0.301029995663981'
  []
[]

[Postprocessors]
  [A]
    type = ElementAverageValue
    variable = A
  []
  [B]
    type = ElementAverageValue
    variable = B
  []
  [C]
    type = ElementAverageValue
    variable = C
  []
  [C1]
    type = ElementAverageValue
    variable = C1
  []
  [objective]
    type = ParsedPostprocessor
    expression = '${objective}'
    pp_names = 'A B C C1'
  []
[]

[Executioner]
  type = Transient
  scheme = bdf2
  solve_type = NEWTON
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  dt = 0.01
  end_time = 1
[]

[Outputs]
  # console = false
[]
