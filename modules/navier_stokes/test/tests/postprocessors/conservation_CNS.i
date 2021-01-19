[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0
  xmax = 10.0
  ymax = 10
  nx = 10
  ny = 10
[]

[Modules]
  [FluidProperties]
    [ideal_gas]
      type = IdealGasFluidProperties
      gamma = 1.4
    []
  []

  [CompressibleNavierStokes]
    # steady-state or transient
    equation_type = transient

    # fluid
    fluid_properties = ideal_gas

    # boundary conditions
    no_penetration_boundary = 'top bottom'
    static_pressure_boundary = 'right'
    static_pressure = '100000' # Pa
    stagnation_boundary = 'left'
    stagnation_flow_direction = '1 0'
    stagnation_pressure = '100001' #Pa
    stagnation_temperature = 300

    # variable types, scalings and initial conditions
    family = LAGRANGE
    order = FIRST
    total_energy_scaling = 9.869232667160121e-6
    initial_pressure = 101325.
    initial_temperature = 300.
    initial_velocity = '1 0 0'
  []
[]


# [Kernels]
#   [temperature_conduction]
#     type = ADHeatConduction
#     variable = temperature
#     thermal_conductivity = '1'
#   []
# []


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
[]


[Preconditioning]
  [SMP_PJFNK]
    type = SMP
    full = true
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


[BCs]
  [no_slip]
    type = DirichletBC
    variable = rhov
    value = 0
    boundary = 'top bottom'
  []
[]


[Postprocessors]
  [inlet_mass_vol]
    type = VolumetricFlowRate
    boundary = left
    vel_x = vel_x
    vel_y = vel_y
    advected_quantity = 1
  []
  [inlet_mass]
    type = MassFlowRate
    boundary = left
    vel_x = vel_x
    vel_y = vel_y
    rho = 1
  []
  [outlet_mass]
    type = MassFlowRate
    boundary = right
    vel_x = vel_x
    vel_y = vel_y
    rho = 1
  []
  [inlet_momentum]
    type = MomentumFlowRate
    momentum_component = x
    boundary = left
    vel_x = vel_x
    vel_y = vel_y
    rho_u = rhou
    rho_v = rhov
    rho = 1
  []
  [outlet_momentum]
    type = MomentumFlowRate
    momentum_component = x
    boundary = right
    vel_x = vel_x
    vel_y = vel_y
    rho_u = rhou
    rho_v = rhov
    rho = 1
  []
  [inlet_advected_energy]
    type = EnergyFlowRate
    boundary = left
    vel_x = vel_x
    vel_y = vel_y
    energy = temperature
  []
  [outlet_advected_energy]
    type = EnergyFlowRate
    boundary = right
    vel_x = vel_x
    vel_y = vel_y
    energy = temperature
  []
[]


[Outputs]
  interval = 10
  exodus = true
[]


[Debug]
  show_actions = true
[]
