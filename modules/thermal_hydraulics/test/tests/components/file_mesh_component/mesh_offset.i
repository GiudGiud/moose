# This file generates the mesh for the FileMeshComponent test.

[Mesh]
  # second_order = true
  [block_1st]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 1
    ny = 1
    xmin = 0
    xmax = 5.0
    ymin = 2.0
    ymax = 3.0
    subdomain_ids = '3'
  []
  [block_2nd]
    type = FileMeshGenerator
    file = 'mesh_in.e'
  []

  [both]
    type = CombinerGenerator
    inputs = 'block_1st block_2nd'
  []
[]


[Variables]
  [v_x]
    order = SECOND
    # order = FIRST
    family = LAGRANGE
    block = 0
  []
[]

[Problem]
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
  [Quadrature]
    custom_blocks = '0 3'
    custom_orders = 'first second'
    element_order = third
  []
[]
