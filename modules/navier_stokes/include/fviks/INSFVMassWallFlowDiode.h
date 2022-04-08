//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSFVAdvectionInterfaceKernel.h"

/**
 * An interface that only lets flow through in one direction. In that direction it allows an
 * advective mass flux. In the other, it becomes a wall boundary condition.
 */
class INSFVMassWallFlowDiode : public INSFVAdvectionInterfaceKernel
{
public:
  static InputParameters validParams();
  INSFVMassWallFlowDiode(const InputParameters & params);

protected:
  virtual ADReal computeQpResidual() override;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// Whether flow diode is inverted
  const int _inverted;
};
