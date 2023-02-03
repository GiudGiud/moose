# This test solves two identical heat conduction problems, one created with THM
# components, and one with the constituent lower-level objects and FileMeshComponent.

rho = 1000
cp = 500
k = 15

initial_T = 300

D_h = 5

[GlobalParams]
  gravity_vector = '0 -9.81 0'

  initial_p = 1e6
  initial_T = 300
  initial_vel = 1

  closures = simple_closures

  pressure = pressure
[]

###########################
# Multi-D physics
###########################

[Variables]
  [v_x]
    order = SECOND
    # order = FIRST
    family = LAGRANGE
    block = 'hs_external:block_a'
  []
  [v_y]
    order = SECOND
    # order = FIRST
    family = LAGRANGE
    block = 'hs_external:block_a'
  []
  [pressure]
    order = FIRST
    # order = CONSTANT
    family = LAGRANGE
    # family = MONOMIAL
    block = 'hs_external:block_a'
  []
[]

[Kernels]
  [mass]
    type = INSMass
    variable = pressure
    u = v_x
    v = v_y
    pressure = pressure
    rho_name = rho_ns
    mu_name = mu_ns
  []
  [x_momentum_space]
    type = INSMomentumLaplaceForm
    variable = v_x
    u = v_x
    v = v_y
    pressure = pressure
    component = 0
    rho_name = rho_ns
    mu_name = mu_ns
  []
  [y_momentum_space]
    type = INSMomentumLaplaceForm
    variable = v_y
    u = v_x
    v = v_y
    pressure = pressure
    component = 1
    rho_name = rho_ns
    mu_name = mu_ns
  []
[]

[BCs]
  [x_no_slip]
    type = DirichletBC
    variable = v_x
    boundary = 'hs_external:top hs_external:bottom'
    value = 0.0
  []
  [y_no_slip]
    type = DirichletBC
    variable = v_y
    boundary = 'hs_external:left hs_external:top hs_external:bottom'
    value = 0.0
  []
[]

###########################
# System
###########################

[Components]
  [hs_external]
    type = FileMeshComponent
    file = 'mesh_in.e'
    position = '0 0 0'
  []


  # [left_wall]
  #   type = SolidWall1Phase
  #   input = pipe:in
  # []
  [inlet]
    type = InletDensityVelocity1Phase
    input = 'pipe:in'
    rho = ${rho}
    vel = 1
  []

  [pipe]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '1 0 0'
    length = 1
    n_elems = 1

    A = 1.0e-4
    D_h = dh_fn

    f = 0.0

    fp = eos
  []

  # [right_wall]
  #   type = SolidWall1Phase
  #   input = pipe:out
  # []
  [outlet]
    type = Outlet1Phase
    input = 'pipe:out'
    p = 1e5
  []
[]

###########################
# Coupling NS -> THM
###########################

[Postprocessors]
  [pressure_pp]
    type = SideAverageValue
    variable = pressure
    boundary = hs_external:left
    execute_on = LINEAR
  []
  [velocity_pp]
    type = SideAverageValue
    variable = v_x
    boundary = hs_external:right
    execute_on = LINEAR
  []
[]

[ControlLogic]
  # [outlet_pressure]
  #   type = ParsedFunctionControl
  #   function = 'pressure_pp'
  #   symbol_names = 'pressure_pp'
  #   symbol_values = 'pressure_pp'
  # []
  # [inlet_velocity]
  #   type = ParsedFunctionControl
  #   function = '2'
  #   symbol_names = 'v'
  #   symbol_values = 'velocity_pp'
  # []

  # [set_outlet_pressure]
  #   type = SetComponentRealValueControl
  #   component = outlet
  #   parameter = p
  #   value = outlet_pressure:value
  # []
  # [set_inlet_velocity]
  #   type = SetComponentRealValueControl
  #   component = inlet
  #   parameter = vel
  #   value = inlet_velocity:value
  # []
[]

[Controls]
  [set_outlet_pressure]
    type = RealFunctionControl
    parameter = 'Components/outlet/p'
    function = 'pressure_pp_fun'
    # execute_on = 'initial timestep_begin'
    execute_on = LINEAR
  []
  [set_inlet_velocity]
    type = RealFunctionControl
    parameter = 'Components/inlet/vel'
    function = 'velocity_pp_fun'
    # execute_on = 'initial timestep_begin'
    execute_on = LINEAR
  []
[]

[Functions]
  [pressure_pp_fun]
    type = ParsedFunction
    expression = 'pressure_pp'
    symbol_names = 'pressure_pp'
    symbol_values = 'pressure_pp'
  []
  [velocity_pp_fun]
    type = ParsedFunction
    expression = 'velocity_pp'
    symbol_names = 'velocity_pp'
    symbol_values = 'velocity_pp'
  []
[]

###########################
# Coupling THM -> NS
###########################

[BCs]
  [x_inlet]
    type = PostprocessorDirichletBC
    variable = v_x
    boundary = 'hs_external:left'
    postprocessor = pipe_out_v
  []
  [p_outlet]
    type = PostprocessorDirichletBC
    variable = pressure
    boundary = 'hs_external:right'
    postprocessor = pipe_in_p
  []
[]

[Postprocessors]
  [m_dot_out]
    type = ADFlowBoundaryFlux1Phase
    boundary = 'outlet'
    equation = mass
    outputs = none
    execute_on = LINEAR
  []
  [pipe_out_v]
    type = ParsedPostprocessor
    function = 'm_dot_out / rho'
    pp_names = 'm_dot_out'
    constant_names = 'rho'
    constant_expressions = '${rho}'
    execute_on = LINEAR
  []
  [pipe_in_p]
    type = SideAverageValue
    variable = p
    boundary = inlet
    execute_on = LINEAR
  []
[]

###########################
# Properties
###########################

[Materials]
  [prop_mat]
    type = ADGenericConstantMaterial
    prop_names = 'density specific_heat thermal_conductivity'
    block = 'pipe'
    prop_values = '${rho} ${cp} ${k}'
  []
  [const_material]
    type = GenericConstantMaterial
    block = 'hs_external:block_a'
    prop_names = 'rho_ns mu_ns'
    prop_values = '${rho} 1'
  []
[]

[FluidProperties]
  # [eos]
  #   type = StiffenedGasFluidProperties
  #   gamma = 2.35
  #   cv = 1816.0
  #   q = -1.167e6
  #   p_inf = 1.0e9
  #   q_prime = 0
  # []
  [eos]
    type = IdealGasFluidProperties
    # type = Water97FluidProperties
    # gamma = 2.35
    # cv = 1816.0
    # q = -1.167e6
    # p_inf = 1.0e9
    # q_prime = 0
    # rho_0 = ${rho}
    # e_0 = 1000
    # T_0 = 300
    # p_0 = 1e5
    # a2 = 0
    # beta = 0
  []
[]

[Closures]
  [simple_closures]
    type = Closures1PhaseSimple
  []
[]

[Functions]
  [inlet_func]
    type = ParsedFunction
    expression = '-4 * (y - 0.5)^2 + 1'
  []
  [dh_fn]
    type = ConstantFunction
    value = ${D_h}
  []
[]

###########################
# Numerical solve
###########################

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
  num_steps = 10
  # abort_on_solve_fail = true
  nl_abs_tol = 1e-12

  line_search = none

  solve_type = 'NEWTON'

  # petsc_options_iname = '-pc_type -pc_factor_shift -pc_factor_mat_solver_type'
  # petsc_options_value = 'lu NONZERO mumps'
  petsc_options = '-pc_svd_monitor'
  petsc_options_iname = '-pc_type -pc_factor_shift'
  petsc_options_value = 'svd NONZERO'
  # petsc_options = '-snes_test_jacobian'
  # petsc_options_iname = '-snes_test_error'
  # petsc_options_value = '1e-8'
[]

###########################
# Postprocessing
###########################

[Outputs]
  [exodus]
    type = Exodus
    file_base = 'file_mesh_component'
  []
[]

