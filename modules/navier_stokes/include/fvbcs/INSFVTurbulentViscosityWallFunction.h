

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FVDirichletBCBase.h"
#include "FVFluxBC.h"

/**
 * Applies a Dirichlet boundary condition with a value prescribed by a function
 */
class INSFVTurbulentViscosityWallFunction : public FVDirichletBCBase
{
public:
  INSFVTurbulentViscosityWallFunction(const InputParameters & parameters);

  static InputParameters validParams();

  ADReal boundaryValue(const FaceInfo & fi) const override;

private:
  /// the dimension of the simulation
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<ADReal> & _u_var;
  /// y-velocity
  const Moose::Functor<ADReal> * _v_var;
  /// z-velocity
  const Moose::Functor<ADReal> * _w_var;

  /// Density
  const Moose::Functor<ADReal> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<ADReal> & _mu;
  /// Trubulent Dynamic viscosity
  const Moose::Functor<ADReal> & _mu_t;

  /// Turbulent kinetic energy
  const Moose::Functor<ADReal> & _k;

  /// C_mu turbulent coefficient
  const Moose::Functor<ADReal> & _C_mu;

  /// Linearzied computation of y_plus ?
  const bool _linearized_yplus;

  /// Maximum mixing length allowed for the domain
  const Real _min_mixing_length;

  /// Non-equilibrium wall treatement ?
  const bool _non_equilibrium_treatement;
};
