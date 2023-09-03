//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "kEpsilonViscosityFunctorMaterial.h"
#include "NavierStokesMethods.h"
#include "NonlinearSystemBase.h"
#include "libmesh/nonlinear_solver.h"

registerMooseObject("NavierStokesApp", kEpsilonViscosityFunctorMaterial);

InputParameters
kEpsilonViscosityFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription(
      "Calculates the turbulent viscosity according to the k-epsilon model.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>("k", "Coupled turbulent kinetic energy.");
  params.addRequiredParam<MooseFunctorName>("epsilon",
                                            "Coupled turbulent kinetic energy dissipation rate.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Density");
  params.addRequiredParam<MooseFunctorName>("mu", "Dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("C_mu", "Coupled turbulent kinetic energy closure.");
  params.addRequiredParam<std::vector<BoundaryName>>("walls",
                                                     "Boundaries that correspond to solid walls.");
  params.addParam<bool>(
      "linearized_yplus",
      false,
      "Boolean to indicate if yplus must be estimated locally for the blending functions.");

  params.addParam<bool>(
      "constrain_time_scale", false, "Whether to limit k / epsilon to a minimum value");
  params.addParam<Real>("max_mixing_length",
                        1e10,
                        "Maximum mixing length allowed for the domain - adjust if seeking for "
                        "realizable k-epsilon answer.");
  params.addParam<bool>(
      "wall_treatment", true, "Activate wall treatment by adding wall functions.");
  params.addParam<bool>("non_equilibrium_treatment",
                        false,
                        "Use non-equilibrium wall treatment (faster than standard wall treatment)");
  params.addParam<Real>("rf", 1.0, "Relaxation factor.");
  params.addParam<unsigned int>("iters_to_activate",
                                1,
                                "Number of nonlinear iterations needed to activate the source in "
                                "the turbulent kinetic energy.");
  params.addParam<Real>("mu_t_inital", 10.0, "Initial value for the turbulent viscosity.");
  MooseEnum relaxation_method("time nl", "nl");
  params.addParam<MooseEnum>(
      "relaxation_method",
      relaxation_method,
      "The method used for relaxing the turbulent kinetic energy production. "
      "'nl' = previous nonlinear iteration and 'time' = previous timestep.");
  params.addParam<Real>(
      "damper",
      0.0,
      "Exponential damping factor for nonlinear iterations. Inactive by default with a value of 0");
  return params;
}

kEpsilonViscosityFunctorMaterial::kEpsilonViscosityFunctorMaterial(const InputParameters & params)
  : FunctorMaterial(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(getFunctor<ADReal>("v")),
    _v_var(params.isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _w_var(params.isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _k(getFunctor<ADReal>("k")),
    _epsilon(getFunctor<ADReal>("epsilon")),
    _rho(getFunctor<ADReal>(NS::density)),
    _mu(getFunctor<ADReal>("mu")),
    _C_mu(getFunctor<ADReal>("C_mu")),
    _wall_boundary_names(getParam<std::vector<BoundaryName>>("walls")),
    _linearized_yplus(getParam<bool>("linearized_yplus")),
    _constrain_time_scale(getParam<bool>("constrain_time_scale")),
    _max_mixing_length(getParam<Real>("max_mixing_length")),
    _wall_treatment(getParam<bool>("wall_treatment")),
    _non_equilibrium_treatment(getParam<bool>("non_equilibrium_treatment")),
    _rf(getParam<Real>("rf")),
    _iters_to_activate(getParam<unsigned int>("iters_to_activate")),
    _mu_t_inital(getParam<Real>("mu_t_inital")),
    _relaxation_method(getParam<MooseEnum>("relaxation_method")),
    _damper(getParam<Real>("damper")),
    _mu_t_old(getFunctor<ADReal>("mu_t_old"))
{
  _max_viscosity_value = 0.0;

  addFunctorProperty<ADReal>(
      "mu_t",
      [this](const auto & r, const auto & t) -> ADReal
      {
        unsigned int current_nl_iteration = _fe_problem.currentNonlinearSystem()
                                                .nonlinearSolver()
                                                ->get_current_nonlinear_iteration_number();

        bool wall_bounded = false;
        Real min_wall_dist = std::numeric_limits<Real>::max();
        Point loc_normal(0, 0, 0);

        // We need to know:
        // - if we are right by a wall (first cell)
        // - how far we are from the wall
        // - the normal to the wall
        getWallData(r, wall_bounded, min_wall_dist, loc_normal);

        auto time_scale = _k(r, t) / _epsilon(r, t);

        if (_constrain_time_scale)
        {
          auto min_constraint_time_scale = 0.6 * std::sqrt(_mu(r, t) / _rho(r, t) / _epsilon(r, t));
          time_scale = std::max(min_constraint_time_scale.value(), time_scale.value());
        }

        auto mu_t =
            _rho(r, t).value() * _C_mu(r, t).value() * _k(r, t).value() * time_scale.value();

        // Wall treatment
        auto mu_t_wall = mu_t;
        if (wall_bounded && _wall_treatment)
        {
          ADRealVectorValue velocity(_u_var(r, t));
          if (_v_var)
            velocity(1) = (*_v_var)(r, t);
          if (_w_var)
            velocity(2) = (*_w_var)(r, t);

          // Compute the velocity and the velocity component that is parallel to the wall
          // NOTE This is only correct if we are right by a wall
          ADReal parallel_speed = (velocity - velocity * loc_normal * loc_normal).norm();

          // Getting y_plus
          ADReal y_plus, u_tau;
          if (_non_equilibrium_treatment)
          {
            y_plus = _rho(r, t) * std::pow(_C_mu(r, t), 0.25) * std::pow(_k(r, t), 0.5) *
                     min_wall_dist / _mu(r, t);
            auto von_karman_value = (1 / _von_karman + std::log(_E * y_plus));
            u_tau = std::sqrt(std::pow(_C_mu(r, t), 0.25) * std::pow(_k(r, t), 0.5) *
                              parallel_speed / von_karman_value);
          }
          else
          {
            u_tau = this->findUStarLocalMethod(parallel_speed, min_wall_dist);
            y_plus = min_wall_dist * u_tau * _rho(r, t) / _mu(r, t);
          }

          // All y+ wall function
          if (y_plus <= 5.0) // sub-laminar layer
          {
            mu_t_wall = _mu(r, t).value();
          }
          else if (y_plus >= 30.0)
          {
            auto wall_val = Utility::pow<2>(u_tau) * _rho(r, t) * min_wall_dist / parallel_speed;
            mu_t_wall = wall_val.value();
          }
          else
          {
            auto wall_val_log =
                Utility::pow<2>(u_tau) * _rho(r, t) * min_wall_dist / parallel_speed;
            auto blending_function = (y_plus - 5.0) / 25.0;
            auto wall_val =
                blending_function * wall_val_log + (1.0 - blending_function) * _mu(r, t).value();
            mu_t_wall = wall_val.value();
          }
          mu_t = mu_t_wall;
        }

        // Dampen out mu_t using the old value
        unsigned int old_state = 2;
        Moose::SolutionIterationType type = Moose::SolutionIterationType::Time;
        if (_relaxation_method == "nl")
          type = Moose::SolutionIterationType::Nonlinear;
        const Moose::StateArg old_state_arg(old_state, type);

        auto mu_t_old = _mu_t_old(r, old_state_arg).value();
        mu_t = _rf * mu_t + (1.0 - _rf) * mu_t_old;

        ADReal return_value = 0;

        // Dampen out mu_t using the number of non-linear iterations
        if (current_nl_iteration == _iters_to_activate)
          _nl_damping_map[_current_elem] += 1.0;
        if ((_relaxation_method == "nl") && (_damper > 1e-10))
          return_value += _mu_t_inital * std::exp(-_nl_damping_map[_current_elem] / _damper);

        if ((current_nl_iteration < _iters_to_activate))
          return_value = _mu_t_inital;
        else
          return_value = std::abs(mu_t);

        if (return_value > _max_viscosity_value)
          return_value = _max_viscosity_value;

        return return_value;
      });
}

void
kEpsilonViscosityFunctorMaterial::initialSetup()
{
  // Setting up limiter for turbulent viscosity
  for (const auto & elem : _c_fe_problem.mesh().getMesh().element_ptr_range())
  {
    const auto current_argument = makeElemArg(elem);
    const auto state = determineState();
    const auto mu_t = _rho(current_argument, state).value() *
                      _C_mu(current_argument, state).value() *
                      std::pow(_k(current_argument, state).value(), 2) /
                      (std::min(1e-10, _epsilon(current_argument, state).value()));

    if (mu_t > _max_viscosity_value)
      _max_viscosity_value = std::max(mu_t, 1e3);
  }

  // Initialize damper
  for (const auto & elem : _c_fe_problem.mesh().getMesh().element_ptr_range())
    _nl_damping_map[elem] = 0.0;
}

ADReal
kEpsilonViscosityFunctorMaterial::findUStarLocalMethod(const ADReal & u, const Real & dist)
{

  /// Setting up parameters
  const auto state = determineState();
  auto rho = _rho(makeElemArg(_current_elem), state);
  auto mu = _mu(makeElemArg(_current_elem), state);
  auto nu = mu / rho;

  const ADReal a_c = 1 / _von_karman;
  const ADReal b_c = 1 / _von_karman * (std::log(_E * dist / mu) + 1.0);
  const ADReal c_c = u;
  ADReal u_tau = (-b_c + std::sqrt(std::pow(b_c, 2) + 4.0 * a_c * c_c)) / (2.0 * a_c);

  if (_linearized_yplus)
  {
    return u_tau;
  }
  else
  {
    // Newton-Raphson method to solve for u_tau
    Real rel_err;
    for (int i = 0; i < _MAX_ITERS_U_TAU; ++i)
    {
      ADReal residual = u / u_tau - 1 / _von_karman * std::log(_E * dist * u_tau / nu);

      if (residual < _REL_TOLERANCE)
        return u_tau;

      ADReal residual_derivative =
          -1 / u_tau * (u / u_tau + 1 / _von_karman * std::log(_E * dist / nu));
      ADReal new_u_tau = std::max(1e-20, u_tau - residual / residual_derivative);
      u_tau = new_u_tau;
    }

    mooseWarning("Could not find the wall friction velocity (mu: ",
                 mu,
                 " rho: ",
                 rho,
                 " velocity: ",
                 u,
                 " wall distance: ",
                 dist,
                 ") - Relative residual: ",
                 rel_err);

    return u_tau;
  }
}

void
kEpsilonViscosityFunctorMaterial::getWallData(Moose::ElemArg r,
                                              bool & wall_bounded,
                                              Real & min_wall_dist,
                                              Point & loc_normal) const
{
  const auto & elem = *r.elem;
  for (unsigned int i_side = 0; i_side < elem.n_sides(); ++i_side)
  {
    const std::vector<BoundaryID> side_bnds =
        _subproblem.mesh().getBoundaryIDs(_current_elem, i_side);

    for (const BoundaryName & name : _wall_boundary_names)
    {
      BoundaryID wall_id = _subproblem.mesh().getBoundaryID(name);
      for (BoundaryID side_id : side_bnds)
      {
        if (side_id == wall_id)
        {
          const FaceInfo * const fi = _mesh.faceInfo(&elem, i_side);
          Real dist = std::abs((fi->elemCentroid() - fi->faceCentroid()) * fi->normal());

          if (dist < min_wall_dist)
          {
            min_wall_dist = dist;
            loc_normal = fi->normal();
          }
          wall_bounded = true;
        }
      }
    }
  }
}

void
kEpsilonViscosityFunctorMaterial::getWallData(Moose::ElemPointArg, bool &, Real &, Point &) const
{
  mooseError("ElemPointArg overload not implemented");
}
void
kEpsilonViscosityFunctorMaterial::getWallData(Moose::ElemQpArg, bool &, Real &, Point &) const
{
  mooseError("ElemQpArg overload not implemented");
}
void
kEpsilonViscosityFunctorMaterial::getWallData(Moose::FaceArg r,
                                              bool & wall_bounded,
                                              Real & min_wall_dist,
                                              Point & loc_normal) const
{
  const FaceInfo * const fi = r.fi;
  mooseAssert(fi, "We should have a fi");
  const auto side_bnds = fi->boundaryIDs();
  for (const BoundaryName & name : _wall_boundary_names)
  {
    BoundaryID wall_id = _subproblem.mesh().getBoundaryID(name);
    for (BoundaryID side_id : side_bnds)
    {
      if (side_id == wall_id)
      {
        Real dist = std::abs((fi->elemCentroid() - fi->faceCentroid()) * fi->normal());

        if (dist < min_wall_dist)
        {
          min_wall_dist = dist;
          loc_normal = fi->normal();
        }
        wall_bounded = true;
      }
    }
  }
}

void
kEpsilonViscosityFunctorMaterial::getWallData(Moose::ElemSideQpArg, bool &, Real &, Point &) const
{
  mooseError("ElemSideQp overload not implemented");
}

// Talk with Mauricio about:
// - ADReal parallel_speed = (velocity - velocity * loc_normal * loc_normal).norm(); away from the
// wall
// - wall distance is not set outside of the first cell by the wall. We have to use a similar trick
// as Sterling:
//   use a distance auxvariable or some
// - the missing derivatives in mu_t and advection
// - min_wall_dist is the maximum wall distance in the auxkernel
