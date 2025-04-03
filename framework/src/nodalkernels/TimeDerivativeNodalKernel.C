//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TimeDerivativeNodalKernel.h"

registerMooseObject("MooseApp", TimeDerivativeNodalKernel);

InputParameters
TimeDerivativeNodalKernel::validParams()
{
  InputParameters params = TimeNodalKernel::validParams();
  params.addClassDescription(
      "Forms the contribution to the residual and jacobian of the time derivative term from an ODE "
      "being solved at all nodes.");
  params.addCoupledVar("nodal_mass", "If specified, multiplies by nodal mass");
  return params;
}

TimeDerivativeNodalKernel::TimeDerivativeNodalKernel(const InputParameters & parameters)
  : TimeNodalKernel(parameters),
    _nodal_mass(isCoupled("nodal_mass") ? &coupledValue("nodal_mass") : nullptr)
{
}

Real
TimeDerivativeNodalKernel::computeQpResidual()
{
  if (_nodal_mass)
    return _u_dot[_qp] * (*_nodal_mass)[_qp];
  else
    return _u_dot[_qp];
}

Real
TimeDerivativeNodalKernel::computeQpJacobian()
{
  if (_nodal_mass)
    return _du_dot_du[_qp] * (*_nodal_mass)[_qp];
  else
    return _du_dot_du[_qp];
}
