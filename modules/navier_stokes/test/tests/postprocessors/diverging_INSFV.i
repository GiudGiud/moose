mu=1.1
rho=1
advected_interp_method='average'
velocity_interp_method='average'

[Mesh]
  inactive = 'mesh internal_boundary_bot internal_boundary_top'
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1'
    dy = '1 1 1'
    ix = '5'
    iy = '5 5 5'
    subdomain_id = '0
                    1
                    2'
  []
  [internal_boundary_bot]
    type = SideSetsBetweenSubdomainsGenerator
    input = mesh
    new_boundary = 'internal_bot'
    primary_block = 0
    paired_block = 1
  []
  [internal_boundary_top]
    type = SideSetsBetweenSubdomainsGenerator
    input = internal_boundary_bot
    new_boundary = 'internal_top'
    primary_block = 1
    paired_block = 2
  []
  [diverging_mesh]
    type = FileMeshGenerator
    file = 'expansion_quad.e'
  []
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = true
[]

[Variables]
  [u]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1
  []
  [v]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = 1
  []
  [pressure]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [temperature]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxVariables]
  [advected_density]
    order = CONSTANT
    family = MONOMIAL
    fv = true
    initial_condition = ${rho}
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
    flow_boundaries = 'bottom top'
    no_slip_wall_boundaries = 'left right'
  []

  [u_advection]
    type = INSFVMomentumAdvection
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
    flow_boundaries = 'bottom top'
    no_slip_wall_boundaries = 'left right'
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    p = pressure
  []

  [v_advection]
    type = INSFVMomentumAdvection
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
    flow_boundaries = 'bottom top'
    no_slip_wall_boundaries = 'left right'
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    p = pressure
  []

  [temp_advection]
    type = INSFVEnergyAdvection
    vel = 'velocity'
    variable = temperature
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = ${mu}
    rho = ${rho}
    flow_boundaries = 'bottom top'
    no_slip_wall_boundaries = 'left right'
  []
  [temp_source]
    type = FVBodyForce
    variable = temperature
    function = 10
    block = 1
  []
[]

[FVBCs]
  inactive = 'noslip'
  [inlet-v]
    type = FVDirichletBC
    boundary = 'bottom'
    variable = v
    value = 1
  []
  [noslip]  #bad noslip
    type = FVDirichletBC
    boundary = 'left right'
    variable = v
    value = 0
  []
  [inlet-and-walls-u]
    type = FVDirichletBC
    boundary = 'bottom left right'
    variable = u
    value = 0
  []
  [outlet_p]
    type = FVDirichletBC
    boundary = 'top'
    variable = pressure
    value = 0
  []
  [inlet_temp]
    type = FVDirichletBC
    boundary = 'bottom'
    variable = temperature
    value = 300
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
    temperature = 'temperature'
    rho = ${rho}
  []
  [advected_material_property]
    type = ADGenericConstantMaterial
    prop_names = 'advected_rho cp'
    prop_values ='${rho} 1'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-12
[]

[Postprocessors]
  [inlet_mass_variable]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_variable = advected_density
    fv = true
  []
  [inlet_mass_constant]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_variable = ${rho}
    fv = true
  []
  [inlet_mass_matprop]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_mat_prop = 'advected_rho'
    fv = true
  []
  [mid1_mass]
    type = VolumetricFlowRate
    boundary = internal_bot
    vel_x = u
    vel_y = v
    fv = true
  []
  [mid2_mass]
    type = VolumetricFlowRate
    boundary = internal_top
    vel_x = u
    vel_y = v
    fv = true
  []
  [outlet_mass]
    type = BoundaryVolumetricFlowRate
    boundary = top
    vel_x = u
    vel_y = v
    fv = true
  []

  [inlet_momentum_x]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_variable = u
    fv = true
  []
  [mid1_momentum_x]
    type = VolumetricFlowRate
    boundary = internal_bot
    vel_x = u
    vel_y = v
    advected_variable = u
    fv = true
  []
  [mid2_momentum_x]
    type = VolumetricFlowRate
    boundary = internal_top
    vel_x = u
    vel_y = v
    advected_variable = u
    fv = true
  []
  [outlet_momentum_x]
    type = BoundaryVolumetricFlowRate
    boundary = top
    vel_x = u
    vel_y = v
    advected_variable = u
    fv = true
  []

  [inlet_momentum_y]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_variable = v
    fv = true
  []
  [mid1_momentum_y]
    type = VolumetricFlowRate
    boundary = internal_bot
    vel_x = u
    vel_y = v
    advected_variable = v
    fv = true
  []
  [mid2_momentum_y]
    type = VolumetricFlowRate
    boundary = internal_top
    vel_x = u
    vel_y = v
    advected_variable = v
    fv = true
  []
  [outlet_momentum_y]
    type = BoundaryVolumetricFlowRate
    boundary = top
    vel_x = u
    vel_y = v
    advected_variable = v
    fv = true
  []

  [inlet_advected_energy]
    type = BoundaryVolumetricFlowRate
    boundary = bottom
    vel_x = u
    vel_y = v
    advected_mat_prop = 'rho_cp_temp'
    fv = true
  []
  [mid1_advected_energy]
    type = VolumetricFlowRate
    boundary = internal_bot
    vel_x = u
    vel_y = v
    advected_mat_prop = 'rho_cp_temp'
    fv = true
  []
  [mid2_advected_energy]
    type = VolumetricFlowRate
    boundary = internal_top
    vel_x = u
    vel_y = v
    advected_mat_prop = 'rho_cp_temp'
    fv = true
  []
  [outlet_advected_energy]
    type = BoundaryVolumetricFlowRate
    boundary = top
    vel_x = u
    vel_y = v
    advected_mat_prop = 'rho_cp_temp'
    fv = true
  []
[]

[Outputs]
  exodus = true  #change this
  csv = true
  # [console_mass]
  #   type = Console
  #   start_step = 1
  #   show = 'inlet_mass_variable inlet_mass_constant inlet_mass_matprop mid1_mass mid2_mass outlet_mass'
  # []
  # [console_momentum_x]
  #   type = Console
  #   start_step = 1
  #   show = 'inlet_momentum_x mid1_momentum_x mid2_momentum_x outlet_momentum_x'
  # []
  # [console_momentum_y]
  #   type = Console
  #   start_step = 1
  #   show = 'inlet_momentum_y mid1_momentum_y mid2_momentum_y outlet_momentum_y'
  # []
  [console_energy]
    type = Console
    start_step = 1
    show = 'inlet_advected_energy mid1_advected_energy mid2_advected_energy outlet_advected_energy'
  []
[]
