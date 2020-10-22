mu=1
rho=1
advected_interp_method='average'
velocity_interp_method='rc'
force_boundary_execution=false
mass_boundaries_to_force='bottom top'
momentum_boundaries_to_force=''

[Mesh]
  file = diverging.msh
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = true
  coord_type = 'RZ'
[]

[Variables]
  [u]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1e-15
  []
  [v]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1e-15
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[FVKernels]
  [mass]
    type = NSFVKernel
    variable = pressure
    advected_quantity = 1
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
    ghost_layers = 2
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${mass_boundaries_to_force}
  []

  [u_advection]
    type = NSFVKernel
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
    ghost_layers = 2
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
  [u_pressure]
    type = FVMomPressure
    variable = u
    momentum_component = 'x'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
  [u_pressure_rz]
    type = FVMomPressureRZ
    variable = u
    p = pressure
  []

  [v_advection]
    type = NSFVKernel
    variable = v
    advected_quantity = 'rhov'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
    ghost_layers = 2
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
  [v_pressure]
    type = FVMomPressure
    variable = v
    momentum_component = 'y'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    force_boundary_execution = ${force_boundary_execution}
    boundaries_to_force = ${momentum_boundaries_to_force}
  []
[]

[FVBCs]
  active = 'inlet-u inlet-v free-slip-wall-u free-slip-wall-v'
  [inlet-u]
    type = FVDirichletBC
    boundary = 'bottom'
    variable = u
    value = 0
  []
  [inlet-v]
    type = FVDirichletBC
    boundary = 'bottom'
    variable = v
    value = 1
  []
  [free-slip-wall-u]
    type = FVPenaltyFreeSlipBC
    boundary = 'right'
    variable = u
    momentum_component = x
    u = u
    v = v
    penalty = 1e3
  []
  [free-slip-wall-v]
    type = FVPenaltyFreeSlipBC
    boundary = 'right'
    variable = v
    momentum_component = y
    u = u
    v = v
    penalty = 1e3
  []
  [no-slip-wall-u]
    type = FVDirichletBC
    boundary = 'right'
    variable = u
    value = 0
  []
  [no-slip-wall-v]
    type = FVDirichletBC
    boundary = 'right'
    variable = v
    value = 0
  []
  [outlet-p]
    type = FVDirichletBC
    boundary = 'top'
    variable = pressure
    value = 0
  []
[]

[Materials]
  [rho]
    type = ADGenericConstantMaterial
    prop_names = 'rho'
    prop_values = ${rho}
  []
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      1000               lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  csv = true
  [dof]
    type = DOFMap
    execute_on = 'initial'
  []
[]
