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
#include "INSFVMomentumResidualObject.h"

/**
 * An interface that only lets flow through in one direction. In that direction it allows an
 * advective momentum flux. In the other, it does not allow anything through.
 */
class INSFVMomentumWallFlowDiode : public INSFVAdvectionInterfaceKernel,
                                   public INSFVMomentumResidualObject
{
public:
  static InputParameters validParams();
  INSFVMomentumWallFlowDiode(const InputParameters & params);
  void gatherRCData(const Elem &) override final {}
  void gatherRCData(const FaceInfo & fi) override final;

protected:
  virtual ADReal computeQpResidual() override;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// The a coefficient for the element
  ADReal _ae = 0;

  /// The a coefficient for the neighbor
  ADReal _an = 0;

  /// Whether flow diode is inverted
  const int _inverted;
};
