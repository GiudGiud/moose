# This file generates the mesh for the FileMeshComponent test.

[Mesh]
    second_order = true
    [gen_mesh_mg]
      type = FileMeshGenerator
      file = 'mesh_in.e'
    []
  []
  