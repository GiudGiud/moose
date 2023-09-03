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
 * Turbulent kinetic energy source and sink term
 */
class INSFVTKESourceSink : public FVElementalKernel
{
public:
  static InputParameters validParams();

  INSFVTKESourceSink(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

protected:
  /// The dimension of the simulation
  const unsigned int _dim;

  /// x-velocity
  const INSFVVelocityVariable * const _u_var;
  /// y-velocity
  const INSFVVelocityVariable * const _v_var;
  /// z-velocity
  const INSFVVelocityVariable * const _w_var;

  /// epsilon - dissipation rate of TKE
  const Moose::Functor<ADReal> & _epsilon;

  /// Density
  const Moose::Functor<ADReal> & _rho;

  /// Dynamic viscosity
  const Moose::Functor<ADReal> & _mu;

  /// Turbulent dynamic viscosity
  const Moose::Functor<ADReal> & _mu_t;

  /// Wall boundaries
  std::vector<BoundaryName> _wall_boundary_names;

  /// Maximum mixing length allowed for the domain
  const Real _max_mixing_length;

  /// Whether to use a linearized model
  const bool _linearized_model;

  /// Linearization coupled functor
  const Moose::Functor<ADReal> & _linear_variable;

  /// Whether to apply a realizable constraint
  const bool _realizable_constraint;

  /// Local relaxation factor
  const Real _rf;

  /// Whether to not apply an equilibrium treatment
  const bool _non_equilibrium_treatment;

  /// C_mu constant
  Real _C_mu;

  /// Maps for wall treatment
  std::map<const Elem *, bool> _wall_bounded;
  std::map<const Elem *, std::vector<Real>> _dist;
  std::map<const Elem *, std::vector<Point>> _normal;

  /// Storing current time
  Real _loc_dt;
  std::map<const Elem *, Real> _previous_production;
  std::map<const Elem *, Real> _previous_destruction;

  /// Constants of the method
  static constexpr Real _von_karman{0.4187};

  /// Storing a time point for XXXXXX
  Real _stored_time;

  /// Relaxation method for production and destruction
  const MooseEnum _relaxation_method;

  /// Number of iterations needed to activate the source in the k epsilon model
  unsigned int _iters_to_activate;

  /// Top bounds for turbulent production and destruction
  Real _top_production_bound;
  Real _top_destruction_bound;
};
