[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 10.0
    nx = 16
    ny = 16
  []
[]

[Variables]
  [velocity]
    family = LAGRANGE_VEC
  []
  [p]
  []
  [temperature]
  []
[]

[AuxVariables]
  [vel_x]
    family = LAGRANGE
  []
  [vel_y]
    family = LAGRANGE
  []
  [energy]
  []
[]

[ICs]
  [velocity]
    type = VectorConstantIC
    x_value = 1e-15
    y_value = 1e-15
    variable = velocity
  []
[]

[Kernels]
  # Mass equation
  [mass]
    type = INSADMass
    variable = p
  []
  [mass_pspg]
    type = INSADMassPSPG
    variable = p
  []

  # Momentum equation
  [momentum_convection]
    type = INSADMomentumAdvection
    variable = velocity
  []
  [momentum_viscous]
    type = INSADMomentumViscous
    variable = velocity
  []
  [momentum_pressure]
    type = INSADMomentumPressure
    variable = velocity
    p = p
    integrate_p_by_parts = true
  []
  [momentum_supg]
    type = INSADMomentumSUPG
    variable = velocity
    velocity = velocity
  []

  # Energy equation
  [temperature_advection]
    type = INSADEnergyAdvection
    variable = temperature
  []
  [temperature_conduction]
    type = ADHeatConduction
    variable = temperature
    thermal_conductivity = 'k'
  []
  [temperature_source]
    type = INSADEnergySource
    variable = temperature
    source_function = 100
  []
  [temperature_supg]
    type = INSADEnergySUPG
    variable = temperature
    velocity = velocity
  []
[]

[AuxKernels]
  [extract_x]
    type = VectorVariableComponentAux
    vector_variable = velocity
    variable = vel_x
    component = x
  []
  [extract_y]
    type = VectorVariableComponentAux
    vector_variable = velocity
    variable = vel_y
    component = y
  []
  # [compute_energy]
  #   type =
  #   variable = energy
  # []
[]

[BCs]
  # [no_slip]
  #   type = VectorFunctionDirichletBC
  #   variable = velocity
  #   boundary = 'right left'
  # []

  [inlet]
    type = VectorFunctionDirichletBC
    variable = velocity
    boundary = 'bottom'
    function_y = 1
  []

  [pressure_BC]
    type = DirichletBC
    variable = p
    boundary = 'top'
    value = 0
  []

  [temperature_bot]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 100
  []
  # [temperature_top]
  #   type = NeumannBC
  #   variable = temperature
  #   boundary = 'top'
  #   # value = 100
  # []
[]

[Materials]
  [const]
    type = ADGenericConstantMaterial
    prop_names = 'rho mu cp k'
    prop_values = '1  1  1  .1'
  []
  [ins_mat]
    type = INSADStabilized3Eqn
    velocity = velocity
    pressure = p
    temperature = temperature
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -sub_pc_factor_levels -ksp_gmres_restart'
  petsc_options_value = 'asm      6                     200'
  line_search = 'none'
  nl_rel_tol = 1e-12
  nl_max_its = 6
[]

[Postprocessors]
  [inlet_mass]
    type = MassFlowRate
    boundary = bottom
    vel_x = vel_x
    vel_y = vel_y
    rho = 1
  []
  [outlet_mass]
    type = MassFlowRate
    boundary = top
    vel_x = vel_x
    vel_y = vel_y
    rho = 1
  []
  [inlet_advected_energy]
    type = EnergyFlowRate
    boundary = bottom
    vel_x = vel_x
    vel_y = vel_y
    energy = temperature
  []
  [outlet_advected_energy]
    type = EnergyFlowRate
    boundary = top
    vel_x = vel_x
    vel_y = vel_y
    energy = temperature
  []
  [inlet_diffused_energy]
    type = ADSideDiffusiveFluxIntegral
    boundary = bottom
    diffusivity = 'k'
    variable = temperature
  []
  [outlet_diffused_energy]
    type = ADSideDiffusiveFluxIntegral
    boundary = top
    diffusivity = 'k'
    variable = temperature
  []
  # [left_diffused_energy]
  #   type = ADSideDiffusiveFluxIntegral
  #   boundary = right
  #   diffusivity = 'k'
  #   variable = temperature
  # []
[]

[Outputs]
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
[]
