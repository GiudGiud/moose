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
class QuadrilateralExpansionConeMeshGenerator : public MeshGenerator
{
public:
  static InputParameters validParams();

  QuadrilateralExpansionConeMeshGenerator(const InputParameters & parameters);

  std::unique_ptr<MeshBase> generate() override;

protected:
  const Point _apex;
  const Real _inner_radius;
  const Real _outer_radius;
  const Real _theta_min;
  const Real _theta_max;
  const unsigned int _num_radial_layers;
  const unsigned int _num_angular_divisions;
  const Real _radial_expansion;
  const subdomain_id_type _block_id;
};
