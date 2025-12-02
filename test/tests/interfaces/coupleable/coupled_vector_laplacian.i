# Test for coupledVectorValuesOld for constant monomials

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 1
    ny = 1
    nz = 1
  []
  order = SECOND
[]

[ICs]
  [ics]
    type = VectorFunctionIC
    variable = var
    function_x = 'y * (x + 0.2 * y + 0.5 * z)'
    function_y = 'y * (2*(x + 2*y + z))'
    function_z = 'z * 3*(x + y + 3*z)'
  []
[]

[AuxVariables]
  [var]
    order = SECOND
    family = MONOMIAL_VEC
  []

  [var_xx_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_xy_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_xz_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_yx_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_yy_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_yz_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_zx_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_zy_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
  [var_zz_lapl]
    order = FIRST
    family = MONOMIAL_VEC
  []
[]

[AuxKernels]
  [lapl_x_comp]
    type = CoupledVectorAux
    variable = var_xx_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '0'
    component_2 = '0'
  []
  [lapl_xy_comp]
    type = CoupledVectorAux
    variable = var_xy_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '0'
    component_2 = '1'
  []
  [lapl_xz_comp]
    type = CoupledVectorAux
    variable = var_xz_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '0'
    component_2 = '2'
  []
  [lapl_yx_comp]
    type = CoupledVectorAux
    variable = var_yx_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '1'
    component_2 = '0'
  []
  [lapl_yy_comp]
    type = CoupledVectorAux
    variable = var_yy_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '1'
    component_2 = '1'
  []
  [lapl_yz_comp]
    type = CoupledVectorAux
    variable = var_yz_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '1'
    component_2 = '2'
  []
  [lapl_zx_comp]
    type = CoupledVectorAux
    variable = var_zx_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '2'
    component_2 = '0'
  []
  [lapl_zy_comp]
    type = CoupledVectorAux
    variable = var_zy_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '2'
    component_2 = '1'
  []
  [lapl_zz_comp]
    type = CoupledVectorAux
    variable = var_zz_lapl
    coupled = 'var'
    operator = 'laplacian'
    component_1 = '2'
    component_2 = '2'
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
