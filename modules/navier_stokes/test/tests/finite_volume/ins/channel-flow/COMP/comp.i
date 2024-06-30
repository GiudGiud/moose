[Mesh]
  [gen]
    type = FileMeshGenerator
    file = '10s.e'
    use_for_exodus_restart = true
  []
[]

[AuxVariables]
  [v1x]
    type = MooseVariableFVReal
    initial_from_file_timestep = LATEST
    initial_from_file_var = 'vel_x'
  []
  [v1y]
    type = MooseVariableFVReal
    initial_from_file_timestep = LATEST
    initial_from_file_var = 'vel_y'
  []
  [p1]
    type = MooseVariableFVReal
    initial_from_file_timestep = LATEST
    initial_from_file_var = 'pressure'
  []
  [v2x]
    type = MooseVariableFVReal
  []
  [v2y]
    type = MooseVariableFVReal
  []
  [p2]
    type = MooseVariableFVReal
  []
  [compx]
    type = MooseVariableFVReal
  []
  [compy]
    type = MooseVariableFVReal
  []
  [compp]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [v2x]
    type = SolutionAux
    solution = uo
    variable = v2x
    from_variable = 'vel_x'
  []
  [v2y]
    type = SolutionAux
    solution = uo
    variable = v2y
    from_variable = 'vel_y'
  []
  [p2]
    type = SolutionAux
    solution = uo
    variable = p2
    from_variable = 'pressure'
  []
  [compx]
    type = ParsedAux
    expression = '(v2x - v1x)'
    coupled_variables = 'v2x v1x'
    variable = 'compx'
  []
  [compy]
    type = ParsedAux
    expression = '(v2y - v1y)'
    coupled_variables = 'v2y v1y'
    variable = 'compy'
  []
  [compp]
    type = ParsedAux
    expression = '(p2 - p1)'
    coupled_variables = 'p2 p1'
    variable = 'compp'
  []
[]

[UserObjects]
  [uo]
    type = SolutionUserObject
    mesh = '01s.e'
    timestep = LATEST
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Postprocessors]
  [min_x]
    type = ElementExtremeValue
    value_type = min
    variable = compx
  []
  [min_y]
    type = ElementExtremeValue
    value_type = min
    variable = compy
  []
  [min_p]
    type = ElementExtremeValue
    value_type = min
    variable = compp
  []
  [max_x]
    type = ElementExtremeValue
    variable = compx
  []
  [max_y]
    type = ElementExtremeValue
    variable = compy
  []
  [max_p]
    type = ElementExtremeValue
    variable = compp
  []
  [l2xdif]
    type = ElementL2Norm
    variable = compx
    execute_on = 'initial timestep_end'
  []
  [l2ydif]
    type = ElementL2Norm
    variable = compy
    execute_on = 'initial timestep_end'
  []
  [l2pdif]
    type = ElementL2Norm
    variable = compp
    execute_on = 'initial timestep_end'
  []
[]

[Outputs]
  exodus = true
[]
