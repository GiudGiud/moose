mu = 1
rho = 1

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1 0.5 1'
    dy = '0.5 0.5'
    ix = '8 5 8'
    iy = '8 8'
    subdomain_id = '0 1 2
                    1 2 1'
  []

  [diodes_blocking]
    type = ParsedGenerateSideset
    input = cmg
    combinatorial_geometry = 'x>1.4999 & x<1.501 & y>0.0999'
    new_sideset_name = diode_block
  []

  [diodes_flowing]
    type = SideSetsBetweenSubdomainsGenerator
    input = diodes_blocking
    primary_block = 1
    paired_block = 0
    new_boundary = diode_flow
  []

  [top_outlet]
    type = ParsedGenerateSideset
    input = diodes_flowing
    combinatorial_geometry = 'x>2.499 & y>0.4999'
    new_sideset_name = top_right
  []

  [bottom_outlet]
    type = ParsedGenerateSideset
    input = top_outlet
    combinatorial_geometry = 'x>2.499 & y<0.50001'
    new_sideset_name = bottom_right
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = 'upwind'
  velocity_interp_method = 'rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 1e-6
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    rho = ${rho}
    boundaries_to_avoid = 'diode_flow diode_block'
  []

  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    rho = ${rho}
    momentum_component = 'x'
    boundaries_to_avoid = 'diode_flow diode_block'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
    boundaries_to_avoid = 'diode_block'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    rho = ${rho}
    momentum_component = 'y'
    boundaries_to_avoid = 'diode_flow diode_block'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
    boundaries_to_avoid = 'diode_block'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
  []
[]

[FVBCs]
  [walls_u]
    type = INSFVNoSlipWallBC
    variable = vel_x
    boundary = 'top bottom'
    function = 0
  []
  [walls_v]
    type = INSFVNoSlipWallBC
    variable = vel_y
    boundary = 'top bottom'
    function = 0
  []

  [inlet_u]
    type = INSFVInletVelocityBC
    variable = vel_x
    boundary = 'left'
    function = 1
  []
  [inlet_v]
    type = INSFVInletVelocityBC
    variable = vel_y
    boundary = 'left'
    function = 0
  []

  [outlet]
    type = INSFVOutletPressureBC
    variable = pressure
    boundary = 'right'
    function = 1
  []
[]

[FVInterfaceKernels]
  [mass_diode_block]
    type = INSFVMassWallFlowDiode
    variable1 = pressure
    subdomain1 = 1
    subdomain2 = 2
    boundary = diode_block
    rho = ${rho}
    invert_diode = true
  []
  [mass_diode_flow]
    type = INSFVMassWallFlowDiode
    variable1 = pressure
    subdomain1 = 0
    subdomain2 = 1
    boundary = diode_flow
    rho = ${rho}
  []

  [momx_diode_block]
    type = INSFVMomentumWallFlowDiode
    variable1 = vel_x
    subdomain1 = 1
    subdomain2 = 2
    boundary = diode_block
    rho = ${rho}
    momentum_component = 'x'
    invert_diode = true
  []
  [momx_diode_flow]
    type = INSFVMomentumWallFlowDiode
    variable1 = vel_x
    subdomain1 = 0
    subdomain2 = 1
    boundary = diode_flow
    rho = ${rho}
    momentum_component = 'x'
  []
  [momy_diode_block]
    type = INSFVMomentumWallFlowDiode
    variable1 = vel_y
    subdomain1 = 1
    subdomain2 = 2
    boundary = diode_block
    rho = ${rho}
    momentum_component = 'y'
    invert_diode = true
  []
  [momy_diode_flow]
    type = INSFVMomentumWallFlowDiode
    variable1 = vel_y
    subdomain1 = 0
    subdomain2 = 1
    boundary = diode_flow
    rho = ${rho}
    momentum_component = 'y'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'asm      lu           NONZERO                   200'
  line_search = 'none'

  nl_abs_tol = 1e-14
[]

[Postprocessors]
  [mdot_top]
    type = VolumetricFlowRate
    boundary = 'top_right'
    vel_x = vel_x
    vel_y = vel_y
  []
  [mdot_bottom]
    type = VolumetricFlowRate
    boundary = 'bottom_right'
    vel_x = vel_x
    vel_y = vel_y
  []
[]

[Outputs]
  exodus = true
[]
