//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVMomentumGravityForce.h"
#include "FEProblemBase.h"

registerMooseObject("NavierStokesApp", FVMomentumGravityForce);

InputParameters
FVMomentumGravityForce::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription("Computes a body force due to gravity.");
  params.addRequiredParam<Real>("gravity", "Gravity magnitude in variable direction");
  return params;
}

FVMomentumGravityForce::FVMomentumGravityForce(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _gravity_value(getParam<Real>("gravity")),
    _rho(getADMaterialProperty<Real>("rho"))
{
}

ADReal
FVMomentumGravityForce::computeQpResidual()
{
  return -_rho[_qp] * _gravity_value;
}
