//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MultiSystemSolveObject.h"
#include "FEProblemBase.h"
#include "LinearSystem.h"
#include "NonlinearSystem.h"

InputParameters
MultiSystemSolveObject::validParams()
{
  InputParameters params = emptyInputParameters();
  params.addParam<std::vector<SolverSystemName>>(
      "system_names",
      "Names of the solver systems (both linear and nonlinear) that will be solved");
  // Multi-system fixed point
  // Defaults to false because of the difficulty of defining a good multi-system convergence
  // criterion, unless we add a default one to the simulation?
  params.addParam<bool>(
      "multi_system_fixed_point",
      false,
      "Whether to perform fixed point (Picard) iterations between the nonlinear systems.");
  params.addRangeCheckedParam<std::vector<Real>>(
      "multi_system_fixed_point_relaxation_factor",
      {1.0},
      "multi_system_fixed_point_relaxation_factor>0 & multi_system_fixed_point_relaxation_factor<2",
      "Relaxation factor(s) applied to system solution updates during multi-system fixed point "
      "iterations; 1 disables relaxation. If one value is provided it is applied to every system; "
      "otherwise the vector must match the number/order of systems being solved.");
  params.addParam<ConvergenceName>(
      "multi_system_fixed_point_convergence",
      "Convergence object to determine the convergence of the multi-system fixed point iteration.");
  MooseEnum multi_sys_fp_algorithm("relaxation secant", "relaxation");
  params.addParam<MooseEnum>(
      "multi_system_fixed_point_algorithm",
      multi_sys_fp_algorithm,
      "Algorithm used to accelerate the multi-system fixed point iterations. 'relaxation' applies "
      "simple solution under/over-relaxation to each system; 'secant' applies a projected secant "
      "acceleration to the coupled system solutions once per fixed point sweep. "
      "'multi_system_fixed_point_relaxation_factor' under/over-relaxes the update in both cases.");

  params.addParamNamesToGroup(
      "system_names multi_system_fixed_point multi_system_fixed_point_convergence "
      "multi_system_fixed_point_relaxation_factor multi_system_fixed_point_algorithm",
      "Multiple solver system");
  return params;
}

MultiSystemSolveObject::MultiSystemSolveObject(Executioner & ex)
  : SolveObject(ex),
    _using_multi_sys_fp_iterations(getParam<bool>("multi_system_fixed_point")),
    _multi_sys_fp_algorithm(getParam<MooseEnum>("multi_system_fixed_point_algorithm")),
    _multi_sys_fp_convergence(nullptr) // has not been created yet

{
  // Retrieve pointers to all the systems from the problem by default
  const auto & nl_sys_names = _problem.getNonlinearSystemNames();
  const auto & linear_sys_names = _problem.getLinearSystemNames();
  if (!isParamValid("system_names"))
  {
    for (const auto & sys_name : nl_sys_names)
      _systems.push_back(&_problem.getSolverSystem(_problem.solverSysNum(sys_name)));
    for (const auto & sys_name : linear_sys_names)
      _systems.push_back(&_problem.getSolverSystem(_problem.solverSysNum(sys_name)));
    _num_nl_systems = nl_sys_names.size();
  }
  else
  {
    _num_nl_systems = 0;
    // Retrieve pointers to all the user-specified systems in the order that the user specified them
    for (const auto & sys_name : getParam<std::vector<SolverSystemName>>("system_names"))
    {
      if (std::find(nl_sys_names.begin(), nl_sys_names.end(), sys_name) != nl_sys_names.end())
      {
        _systems.push_back(&_problem.getSolverSystem(_problem.solverSysNum(sys_name)));
        _num_nl_systems++;
      }
      else if (std::find(linear_sys_names.begin(), linear_sys_names.end(), sys_name) !=
               linear_sys_names.end())
        _systems.push_back(&_problem.getSolverSystem(_problem.solverSysNum(sys_name)));
      else
        paramError("system_names",
                   "System '" + sys_name +
                       "' was not found in the Problem. Did you forget to declare it in the "
                       "[Problem] block?");
    }
  }

  if (_pars.isParamSetByUser("multi_system_fixed_point_relaxation_factor") &&
      !_using_multi_sys_fp_iterations)
    paramError("Can't use relaxation factors because multisystem fixed point iteration hasn't been "
               "enabled!");
  if (_pars.isParamSetByUser("multi_system_fixed_point_algorithm") &&
      !_using_multi_sys_fp_iterations)
    paramError("multi_system_fixed_point_algorithm",
               "Can't select a multi-system fixed point acceleration algorithm because "
               "multi-system fixed point iterations haven't been enabled!");

  setupMultiSystemFixedPointRelaxationFactors();
}

void
MultiSystemSolveObject::setupMultiSystemFixedPointRelaxationFactors()
{
  _multi_sys_fp_relax_factors =
      getParam<std::vector<Real>>("multi_system_fixed_point_relaxation_factor");
  if (_multi_sys_fp_relax_factors.size() == 1)
    _multi_sys_fp_relax_factors.resize(_systems.size(), _multi_sys_fp_relax_factors[0]);
  else if (_multi_sys_fp_relax_factors.size() != _systems.size())
    paramError("multi_system_fixed_point_relaxation_factor",
               "Must provide either 1 value or " + Moose::stringify(_systems.size()) +
                   " values (one per system in the solve order).");

  // For each solver system; record whether to perform a fixed point transformation. Relaxation is
  // only needed when the factor differs from 1; the secant algorithm transforms every system's
  // solution (the factor then under/over-relaxes the secant update).
  const bool secant = _multi_sys_fp_algorithm == "secant";
  _perform_multi_sys_fp_relaxation.resize(_systems.size(), false);
  for (const auto i : make_range(_systems.size()))
    if (_using_multi_sys_fp_iterations &&
        (secant || !MooseUtils::absoluteFuzzyEqual(_multi_sys_fp_relax_factors[i], 1.0)))
      _perform_multi_sys_fp_relaxation[i] = true;
}
