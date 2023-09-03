//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FunctorMaterial.h"

#include "INSFVVelocityVariable.h"
/**
 * Computes the turbulent viscosity for the k-Epsilon model.
 * Implements two near-wall treatments: equilibrium and non-equilibrium wall functions.
 */
class kEpsilonViscosityFunctorMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  kEpsilonViscosityFunctorMaterial(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  // virtual void meshChanged() override;

  // Local method to find friction velocity
  // Note: this method may be to need reimplemented for each new turbulent model
  ADReal findUStarLocalMethod(const ADReal & u, const Real & dist);

  /**
   * Obtain wall information
   * @param r the spatial argument
   * @param wall_bounded whether the spatial argument is right by a wall
   * @param min_wall_dist distance to the nearest wall
   * @param loc_normal normal of the face
   */
  void getWallData(Moose::ElemArg r,
                   bool & wall_bounded,
                   Real & min_wall_dist,
                   Point & loc_normal) const;
  void getWallData(Moose::ElemPointArg r,
                   bool & wall_bounded,
                   Real & min_wall_dist,
                   Point & loc_normal) const;
  void getWallData(Moose::ElemQpArg r,
                   bool & wall_bounded,
                   Real & min_wall_dist,
                   Point & loc_normal) const;
  void getWallData(Moose::FaceArg r,
                   bool & wall_bounded,
                   Real & min_wall_dist,
                   Point & loc_normal) const;
  void getWallData(Moose::ElemSideQpArg r,
                   bool & wall_bounded,
                   Real & min_wall_dist,
                   Point & loc_normal) const;

  /// the dimension of the simulation
  const unsigned int _dim;

  /// x-velocity
  const Moose::Functor<ADReal> & _u_var;
  /// y-velocity
  const Moose::Functor<ADReal> * const _v_var;
  /// z-velocity
  const Moose::Functor<ADReal> * const _w_var;

  /// Turbulent kinetic energy
  const Moose::Functor<ADReal> & _k;
  /// Turbulent kinetic energy dissipation rate
  const Moose::Functor<ADReal> & _epsilon;

  /// Density
  const Moose::Functor<ADReal> & _rho;
  /// Dynamic viscosity
  const Moose::Functor<ADReal> & _mu;
  /// C-mu closure coefficient
  const Moose::Functor<ADReal> & _C_mu;
  /// Wall boundaries
  std::vector<BoundaryName> _wall_boundary_names;

  /// Whether to use a linearized computation of y_plus
  const bool _linearized_yplus;

  /// Whether to constrain the time scale
  const bool _constrain_time_scale;

  /// Maximum mixing length allowed for the domain
  const Real _max_mixing_length;

  /// Whether to use a wall treatment
  const bool _wall_treatment;

  /// Whether to use a wall treatment for non-equilibrium
  const bool _non_equilibrium_treatment;

  /// Relaxation factor on the turbulent viscosity
  const Real _rf;

  /// Parameters of the wall function method
  /// Maximum number of iterations to find the friction velocity
  static constexpr int _MAX_ITERS_U_TAU{50};
  /// Relative tolerance to find the friction velocity
  static constexpr Real _REL_TOLERANCE{1e-6};

  /// Viscosity limiter - maximum viscosity must be at the wall
  Real _max_viscosity_value;

  /// Number of iterations needed to activate the computation of mu_t
  unsigned int _iters_to_activate;

  /// Initial value for mu_t
  Real _mu_t_initial;

  /// Relaxation method for production and destruction
  const MooseEnum _relaxation_method;

  /// Element Localized Damping: exponent that is increased on every use to increase damping
  std::map<const Elem *, Real> _nl_damping_map;

  /// Exponent for damping based on nonlinear iterations
  Real _damper;

  /// Previous step turbulent viscosity for damping. This must be set externally
  const Moose::Functor<ADReal> & _mu_t_old;

  /// Constants of the method
  static constexpr Real _von_karman{0.4187};
  static constexpr Real _E{9.793};
};
