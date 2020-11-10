//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVMomentumBoussinesqForce.h"
#include "FEProblemBase.h"

registerMooseObject("NavierStokesApp", FVMomentumBoussinesqForce);

InputParameters
FVMomentumBoussinesqForce::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription("Computes a body force due to gravity and density changes.");
  params.addRequiredParam<Real>("gravity", "Gravity magnitude in variable direction");
  params.addRequiredParam<Real>("T_0", "Reference temperature");
  params.addRequiredCoupledVar("T", "Fluid temperature");
  return params;
}

FVMomentumBoussinesqForce::FVMomentumBoussinesqForce(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _gravity_value(getParam<Real>("gravity")),
    _T_reference(getParam<Real>("T_0")),
    _alpha(getADMaterialProperty<Real>("alpha")),
    _T(adCoupledValue("T"))
{
}

ADReal
FVMomentumBoussinesqForce::computeQpResidual()
{
  return - _gravity_value * _alpha[_qp] * (_T[_qp] - _T_reference);
}
