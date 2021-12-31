mu=1.1
rho=1.1
advected_interp_method='average'
velocity_interp_method='rc'

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 1
    dx = '1 1'
    ix = '30 30'
    subdomain_id = '1 2'
  []
  [jump_sideset]
    type = SideSetsBetweenSubdomainsGenerator
    input = mesh
    primary_block = 1
    paired_block = 2
    new_boundary = 'jump'
  []
[]

[Variables]
  [u]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1
  []
  [pressure]
    type = INSFVPressureVariable
  []
[]

[AuxVariables]
  [porosity]
    type = MooseVariableFVReal
  []
[]

[ICs]
  inactive = 'porosity_continuous'
  [porosity_1]
    type = ConstantIC
    variable = porosity
    block = 1
    value = 1
  []
  [porosity_2]
    type = ConstantIC
    variable = porosity
    block = 2
    value = 0.5
  []
  [porosity_continuous]
    type = FunctionIC
    variable = porosity
    block = '1 2'
    function = smooth_jump
  []
[]

[Functions]
  [smooth_jump]
    type = ParsedFunction
    value = '1 - 0.5 * 1 / (1 + exp(-30*(x-1)))'
  []
[]

[FVKernels]
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    vel = 'velocity'
    pressure = pressure
    u = u
    mu = ${mu}
    rho = ${rho}
    porosity = porosity
  []

  [u_advection]
    type = PINSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = 'velocity'
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    pressure = pressure
    u = u
    mu = ${mu}
    rho = ${rho}
    porosity = porosity
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    porosity = porosity
  []
  [u_pressure]
    type = PINSFVMomentumPressureFlux
    variable = u
    pressure = pressure
    porosity = porosity
    momentum_component = 'x'
    boundary_to_skip = 'jump'
  []
[]

[FVInterfaceKernels]
  [jump]
    type = PINSFVPorosityJumpv2
    variable1 = pressure #### Change this??
    # variable2 = pressure
    boundary = 'jump'
    subdomain1 = '1'
    subdomain2 = '2'
    pressure = 'pressure'
    rho = ${rho}
    porosity = porosity
    superficial_vel_x = u
    momentum_component = 'x'
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u
    function = '1'
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = 1
  []
[]

[Materials]
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    pressure = 'pressure'
    rho = ${rho}
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      100                lu           NONZERO'
  line_search = 'none'
[]

[Postprocessors]
  [inlet_p]
    type = SideAverageValue
    variable = 'pressure'
    boundary = 'left'
  []
  [outlet-u]
    type = SideIntegralVariablePostprocessor
    variable = u
    boundary = 'right'
  []
[]

[Outputs]
  exodus = true
  csv = false
[]
