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

[GlobalParams]
  enthalpy = enthalpy
  fluid_properties = ideal_gas
  rho = rho
  rhoE = rhoE
  rhou = rhou
  rhov = rhov
  rho_u = rhou
  rho_v = rhov
  temperature = temperature
  vel_x = vel_x
  vel_y = vel_y
[]

[Variables]
  [vel_x]
    family = LAGRANGE
    order = FIRST
  []
  [vel_y]
    family = LAGRANGE
    order = FIRST
  []
  [p]
  []
  [temperature]
  []
[]

[AuxVariables]
  [enthalpy]
  []
  [rho]
  []
  [rhoE]
  []
  [rhou]
  []
  [rhov]
  []
[]

[Kernels]
  # Mass equation
  [mass]
    type = NSSUPGMass
    variable = p
  []

  # Momentum equation
  [momentum_convection]
    type = MomentumConvectiveFlux
    variable = vel_x
  []
  [momentum_viscous_x]
    type = NSMomentumViscousFlux
    variable = vel_x
    component = 0
  []
  [momentum_viscous_y]
    type = NSMomentumViscousFlux
    variable = vel_y
    component = 1
  []
  # [momentum_pressure]
  #   type = NSADMomentumPressure
  #   variable = velocity
  #   p = p
  #   integrate_p_by_parts = true
  # []
  [momentum_supgx]
    type = NSSUPGMomentum
    variable = vel_x
    component = 0
  []
  [momentum_supgy]
    type = NSSUPGMomentum
    variable = vel_y
    component = 1
  []

  # Energy equation
  [temperature_advection]
    type = TotalEnergyConvectiveFlux
    variable = temperature
  []
  [temperature_conduction]
    type = ADHeatConduction
    variable = temperature
    thermal_conductivity = 'k'
  []
  # [temperature_source]
  #   type = INSADEnergySource
  #   variable = temperature
  #   source_function = 100
  # []
  [temperature_supg]
    type = NSSUPGEnergy
    variable = temperature
  []
[]

[BCs]
  # [no_slip]
  #   type = DirichletBC
  #   variable = vel_y
  #   boundary = 'right left'
    # value = 0
  # []

  [inlet]
    type = DirichletBC
    variable = vel_y
    boundary = 'bottom'
    value  = 1
  []

  [pressure_BC]
    type = DirichletBC
    variable = p
    boundary = 'top'
    value = 0
  []

  [temperature_inlet]
    type = DirichletBC
    variable = temperature
    boundary = 'bottom'
    value = 100
  []
[]

[Modules]
  [./FluidProperties]
    [ideal_gas]
      type = IdealGasFluidProperties
      gamma = 1.4
    []
  []
[]

[Materials]
  [fluid]
    type = Air
    rho = rho
    rhou = rhou
    rhov = rhov
    rhoE = rhoE
    vel_x = vel_x
    vel_y = vel_y
    temperature = temperature
    enthalpy = enthalpy
    # This value is not used in the Euler equations, but it *is* used
    # by the stabilization parameter computation, which it decreases
    # the amount of artificial viscosity added, so it's best to use a
    # realistic value.
    dynamic_viscosity = 0.0
    fluid_properties = ideal_gas
  []
  [thermal_conductivity]
    type = ADGenericConstantMaterial
    prop_names = 'k'
    prop_values = '1'
  []
[]

[Executioner]
  type = Transient

  dtmin = 1.e-5
  dtmax = 1.e-2
  start_time = 0.0
  num_steps = 1000
  nl_rel_tol = 1e-5

  nl_max_its = 5
  l_tol = 1e-4 # Relative linear tolerance for each Krylov solve
  l_max_its = 100 # Number of linear iterations for each Krylov solve

  # Specify the order as FIRST, otherwise you will get warnings in DEBUG mode...
  [Quadrature]
    type = TRAP
    order = FIRST
  []

  # Increase time step
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 5.e-5
    growth_factor = 2
  []
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
  # [wall_diffused_energy]
  #   type = ADSideDiffusiveFluxIntegral
  #   boundary = 'left right'
  #   diffusivity = 'k'
  #   variable = temperature
  # []
[]

[Outputs]
  exodus = true
[]

[Debug]
  show_var_residual_norms = true
  show_actions = true
[]
