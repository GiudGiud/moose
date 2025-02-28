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
class SolidConjugateMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  SolidConjugateMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  unsigned int _nx;
  unsigned int _ny;
  unsigned int _nz;
  Real _x_length;
  Real _y_length;
  Real _z_length;
  subdomain_id_type _block_id;
};
