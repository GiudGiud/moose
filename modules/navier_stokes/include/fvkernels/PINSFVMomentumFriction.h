//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVElementalKernel.h"

/**
 * Imposes a friction force on the momentum equation in porous media.
 */
class PINSFVMomentumFriction : public FVElementalKernel
{
public:
  static InputParameters validParams();
  PINSFVMomentumFriction(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  /// Momentum equation component (x = 0, y = 1, z = 2)
  const unsigned int _component;
  /// The linear friction factor for laminar flow as a material property
  const ADMaterialProperty<Real> * const _linear_friction_matprop;
  /// The quadratic friction factor for turbulent flow as a material property
  const ADMaterialProperty<Real> * const _quadratic_friction_matprop;
  /// Darcy coefficient
  const ADMaterialProperty<RealVectorValue> * const _cL;
  /// Forchheimer coefficient
  const ADMaterialProperty<RealVectorValue> * const _cQ;
  /// Booleans to select the right model
  const bool _use_linear_friction_matprop;
  const bool _use_quadratic_friction_matprop;
  const bool _use_Darcy_friction_model;
  const bool _use_Forchheimer_friction_model;
  /// Porosity to compute the intersitial velocity from the superficial velocity
  const VariableValue & _eps;
  /// Momentum as a material property
  const ADMaterialProperty<Real> * const _momentum;
  /// Constant density, use only with incompressible flow
  const Real _rho;
};