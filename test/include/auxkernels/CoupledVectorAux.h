//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * Coupled auxiliary value of a vector variable
 */
class CoupledVectorAux : public VectorAuxKernel
{
public:
  static InputParameters validParams();

  CoupledVectorAux(const InputParameters & parameters);

  virtual ~CoupledVectorAux() {}

protected:
  virtual RealVectorValue computeValue() override;

  /// Operator being applied on this vector variable and coupled vector variable
  MooseEnum _operator;
  /// Value of the coupled vector variable at each Qp
  const VectorVariableValue & _coupled_val;
  /// Value of the gradient of the coupled vector variable at each Qp
  const VectorVariableGradient & _coupled_grad;
  /// Value of the laplacian of the coupled vector variable at each Qp
  const VectorVariableSecond & _coupled_lapl;
  /// Value of the time derivative of the coupled vector variable at each Qp
  const VectorVariableValue & _coupled_dot;
  /// Components when using a Gradient/second
  unsigned int _comp_1, _comp_2;
};
