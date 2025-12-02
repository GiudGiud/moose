# Test for coupledVectorValuesOld for constant monomials

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
  []
[]

[ICs]
  [ics]
    type = VectorFunctionIC
    variable = var
    function_x = 'x + 0.2 * y + 0.5 * z'
    function_y = '2*(x + 2*y + z)'
    function_z = '3*(x + y + 3*z)'
  []
[]

[AuxVariables]
  [var]
    order = FIRST
    family = MONOMIAL_VEC
  []

  [var_x_grad]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_y_grad]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_z_grad]
    order = FIRST
    family = MONOMIAL_VEC
  []
[]

[AuxKernels]
  [grad_x_comp]
    type = CoupledVectorAux
    variable = var_x_grad
    coupled = 'var'
    operator = 'gradient'
    component_1 = '0'
  []
  [grad_y_comp]
    type = CoupledVectorAux
    variable = var_y_grad
    coupled = 'var'
    operator = 'gradient'
    component_1 = '1'
  []
  [grad_z_comp]
    type = CoupledVectorAux
    variable = var_z_grad
    coupled = 'var'
    operator = 'gradient'
    component_1 = '2'
  []
[]

[Executioner]
  type = Transient
  num_steps=1
[]

[Problem]
  solve = false
[]

[Outputs]
  exodus = true
  hide = 'var'
  execute_on = 'FINAL'
[]
