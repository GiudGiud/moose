# This test solves two identical heat conduction problems, one created with THM
# components, and one with the constituent lower-level objects and FileMeshComponent.

rho = 8000
cp = 500
k = 15

initial_T = 1000
T_left = 1005
T_right = 300
htc_right = 1000

D_h = 5

[GlobalParams]
  gravity_vector = '0 -9.81 0'

  initial_p = 1e6
  initial_T = 453.1
  initial_vel = 0.0

  closures = simple_closures

  pressure = pressure
[]

###########################
# Multi-D physics
###########################

[Variables]
  [temp]
    block = 'hs_external:block_a'
  []
[]

[Kernels]
  [mass]
    type = Diffusion
    variable = temp
  []
[]

[BCs]
  [no_slip]
    type = DirichletBC
    variable = temp
    boundary = 'hs_external:top hs_external:bottom'
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
    # second_order = true
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
# Coupling THM -> NS
###########################

[BCs]
  [inlet]
    type = PostprocessorDirichletBC
    variable = temp
    boundary = 'hs_external:left'
    postprocessor = pipe_out_v
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
  [eos]
    type = StiffenedGasFluidProperties
    gamma = 2.35
    cv = 1816.0
    q = -1.167e6
    p_inf = 1.0e9
    q_prime = 0
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

