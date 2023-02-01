# This test solves two identical heat conduction problems, one created with THM
# components, and one with the constituent lower-level objects and FileMeshComponent.

rho = 8000
cp = 500
k = 15

initial_T = 1000
T_left = 1005
T_right = 300
htc_right = 1000

[Variables]
  # [T_moose]
  #   block = 'hs_external:block_a'
  #   initial_condition = ${initial_T}
  # []
  [./vel_x]
    # order = SECOND
    order = FIRST
    family = LAGRANGE
    block = 'hs_external:block_a'
  [../]
  [./vel_y]
    # order = SECOND
    order = FIRST
    family = LAGRANGE
    block = 'hs_external:block_a'
  [../]
  [./p]
    # order = FIRST
    order = CONSTANT
    # family = LAGRANGE
    family = MONOMIAL
    block = 'hs_external:block_a'
  [../]
[]

[Kernels]
  # [time_derivative]
  #   type = ADHeatConductionTimeDerivative
  #   variable = T_moose
  #   block = 'hs_external:block_a'
  #   density_name = density
  #   specific_heat = specific_heat
  # []
  # [heat_conduction]
  #   type = ADHeatConduction
  #   variable = T_moose
  #   block = 'hs_external:block_a'
  #   thermal_conductivity = thermal_conductivity
  # []
  [./mass]
    type = INSMass
    variable = p
    u = vel_x
    v = vel_y
    pressure = p
  [../]
  [./x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_x
    u = vel_x
    v = vel_y
    pressure = p
    component = 0
  [../]
  [./y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = vel_y
    u = vel_x
    v = vel_y
    pressure = p
    component = 1
  [../]
[]

[BCs]
  # [dirichlet_bc]
  #   type = ADFunctionDirichletBC
  #   variable = T_moose
  #   boundary = 'hs_external:left'
  #   function = ${T_left}
  # []
  # [convection_bc]
  #   type = ADConvectionHeatTransferBC
  #   variable = T_moose
  #   boundary = 'hs_external:right'
  #   T_ambient = ${T_right}
  #   htc_ambient = ${htc_right}
  # []
  [./x_no_slip]
    type = DirichletBC
    variable = vel_x
    boundary = 'hs_external:top hs_external:bottom'
    value = 0.0
  [../]
  [./y_no_slip]
    type = DirichletBC
    variable = vel_y
    boundary = 'hs_external:left hs_external:top hs_external:bottom'
    value = 0.0
  [../]
  [./x_inlet]
    type = FunctionDirichletBC
    variable = vel_x
    boundary = 'hs_external:left'
    function = 'inlet_func'
  [../]
[]

[Materials]
  [prop_mat]
    type = ADGenericConstantMaterial
    prop_names = 'density specific_heat thermal_conductivity'
    prop_values = '${rho} ${cp} ${k}'
  []
  [./const]
    type = GenericConstantMaterial
    # block = 0
    prop_names = 'rho mu'
    prop_values = '1  1'
  [../]
[]

[Components]
  [hs_external]
    type = FileMeshComponent
    file = 'mesh_in.e'
    position = '0 0 0'
    # second_order = true
  []
  [hs]
    type = HeatStructurePlate
    position = '0 0 0'
    orientation = '1 0 0'

    length = 5.0
    n_elems = 10

    names = 'blk'
    widths = '1.0'
    n_part_elems = '2'

    depth = 1.0

    initial_T = ${initial_T}
  []
  [start]
    type = HSBoundarySpecifiedTemperature
    hs = hs
    boundary = 'hs:start'
    T = ${T_left}
  []
  [end]
    type = HSBoundaryAmbientConvection
    hs = hs
    boundary = 'hs:end'
    T_ambient = ${T_right}
    htc_ambient = ${htc_right}
  []
[]

[Functions]
  [./inlet_func]
    type = ParsedFunction
    value = '-4 * (y - 0.5)^2 + 1'
  [../]
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  # scheme = 'bdf2'

  start_time = 0
  dt = 1.0
  num_steps = 5
  abort_on_solve_fail = true

  solve_type = 'NEWTON'

  # petsc_options_iname = '-pc_type -pc_factor_shift'
  # petsc_options_value = 'lu NONZERO'
  petsc_options = '-pc_svd_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift'
  petsc_options_value = 'svd NONZERO'
[]

[Outputs]
  [exodus]
    type = Exodus
    file_base = 'file_mesh_component'
  []
[]

