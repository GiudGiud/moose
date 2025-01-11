//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorizedKernel.h"

/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class VectorizedDiffusion : public VectorizedKernel
{
public:
  static InputParameters validParams();

  VectorizedDiffusion(const InputParameters & parameters);

protected:
  virtual void computeQpResiduals() override;

  virtual void computeQpJacobians() override;
};
