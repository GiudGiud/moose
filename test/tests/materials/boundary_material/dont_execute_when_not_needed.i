[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2

    xmin = 0
    xmax = 1
    nx = 4

    ymin = 0
    ymax = 1
    ny = 3
  []
  [interface]
    type = ParsedGenerateSideset
    input = 'gmg'
    combinatorial_geometry = 'x>0.49999 & x<0.50001'
    new_sideset_name = 'interface'
  []
[]

[Variables]
   [temp]
    initial_condition = 1
  []
[]

[Materials]
  [diff_coeff]
    type = ErrorMaterial
    sidesets_to_error_on = 'left interface'
    block = '0'
  []
[]

[Kernels]
  [heat]
    type = MatCoefDiffusion
    variable = temp
    conductivity = 'matp'
  []
[]

[BCs]
  [right_t]
    type = DirichletBC
    boundary =  right
    value    =  2
    variable =  temp
  []
  [left_t]
    type = NeumannBC
    boundary =  left
    value    =  2
    variable =  temp
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  num_steps = 1
  end_time = 1
[]

[Outputs]
  exodus = true
[]
