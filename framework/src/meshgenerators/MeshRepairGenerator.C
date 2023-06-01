//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MeshRepairGenerator.h"

#include "CastUniquePointer.h"
#include "MooseMeshUtils.h"

#include "libmesh/mesh_modification.h"

registerMooseObject("MooseApp", MeshRepairGenerator);

InputParameters
MeshRepairGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addRequiredParam<MeshGeneratorName>("input", "The mesh we want to modify");

  params.addClassDescription("Fixes the mesh.");

  return params;
}

MeshRepairGenerator::MeshRepairGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters), _input(getMesh("input"))
{
}

std::unique_ptr<MeshBase>
MeshRepairGenerator::generate()
{
  std::unique_ptr<MeshBase> mesh = std::move(_input);
  std::cout << mesh->n_elem() << std::endl;

  // Flip negative elements
  MeshTools::Modification::orient_elements(*mesh);

  // Delete 0 volume elements
  for (auto & elem : mesh->active_element_ptr_range())
    if (elem->volume() < 1e-14)
      mesh->delete_elem(elem);
  mesh->contract();
  std::cout << mesh->n_elem() << std::endl;

  // Delete unsupported elements
  for (auto & elem : mesh->active_element_ptr_range())
    if (elem->dim() != 3)
    {
      // std::cout << "deleting " << elem->get_info() << std::endl;
      mesh->delete_elem(elem);
    }
  mesh->contract();

  std::cout << "Left : " << mesh->n_elem() << std::endl;

  mesh->set_isnt_prepared();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
