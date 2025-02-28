//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MeshGenerator.h"

/**
 * MeshGenerator for re-numbering or re-naming boundaries
 */
class HexagonalAssemblyMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  HexagonalAssemblyMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  Real _side_length;
  unsigned int _num_rows;
  unsigned int _num_cols;
  subdomain_id_type _assembly_block_id;
};
