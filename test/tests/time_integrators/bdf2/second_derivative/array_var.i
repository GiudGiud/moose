[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 3
[]

[Variables]
  [u]
    components = 2
  []
[]

[Kernels]
  [test]
    type = PrintArrayTimeDerivatives
    variable = u
    use_displaced_mesh = false
  []
  [u_diffusion]
    type = ArrayDiffusion
    variable = u
    diffusion_coefficient = u_dc
    use_displaced_mesh = false
  []
  [u_time]
    type = ArrayTimeDerivative
    variable = u
    use_displaced_mesh = false
  []
[]

[ICs]
  [u]
    type = ArrayFunctionIC
    variable = u
    function = '.3*(-x+y-0.1*z) 0.19*(-x+y+z)'
  []
[]

[Materials]
  [u_dc]
    type = GenericConstantArray
    prop_name = u_dc
    prop_value = '1 1'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  dt = 0.1
  dtmin = 0.1
  num_steps = 3
  solve_type = 'PJFNK'
  petsc_options = '-snes_test_jacobian -snes_test_jacobian_view'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'

  # scheme = newmark-beta
  [TimeIntegrator]
    type = BDF2
    compute_second_derivative = true
  []
[]

[Outputs]
  exodus = true
[]
