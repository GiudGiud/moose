//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PrintArrayTimeDerivatives.h"

registerMooseObject("MooseTestApp", PrintArrayTimeDerivatives);

InputParameters
PrintArrayTimeDerivatives::validParams()
{

  InputParameters params = ArrayKernel::validParams();

  return params;
}

PrintArrayTimeDerivatives::PrintArrayTimeDerivatives(const InputParameters & parameters)
  : ArrayKernel(parameters), _u_dot(_var.uDot()), _u_dotdot(_var.uDotDot())
{
}

void
PrintArrayTimeDerivatives::computeQpResidual(RealEigenVector & residual)
{

  std::cerr << "_u_dot:\n" << _u_dot[_qp] << "\n";
  std::cerr << "_u_dotdot:\n" << _u_dotdot[_qp] << "\n";

  return;
}
