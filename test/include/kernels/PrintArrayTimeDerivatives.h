//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ArrayKernel.h"

class PrintArrayTimeDerivatives : public ArrayKernel
{

public:
  PrintArrayTimeDerivatives(const InputParameters & params);

  static InputParameters validParams();

protected:
  const ArrayVariableValue & _u_dot;

  const ArrayVariableValue & _u_dotdot;

  virtual void computeQpResidual(RealEigenVector & residual) override;
};
