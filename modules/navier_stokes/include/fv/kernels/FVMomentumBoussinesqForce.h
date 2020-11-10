//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "NSFVBase.h"

/**
 * Computes a body force due to density changes and gravity
 */
class FVMomentumBoussinesqForce : public FVElementalKernel
{
public:
  FVMomentumBoussinesqForce(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  ADReal computeQpResidual();

  const Real _gravity_value;
  const Real _T_reference;
  const ADMaterialProperty<Real> & _alpha;
  const ADVariableValue & _T;
};
