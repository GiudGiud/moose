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
class GlobalSketchMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  GlobalSketchMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  Real _width;
  Real _height;
  Real _chimney_radius;
  Point _chimney_center;
  subdomain_id_type _solid_block_id;
  subdomain_id_type _chimney_block_id;
  unsigned int _nx;
  unsigned int _ny;
};
