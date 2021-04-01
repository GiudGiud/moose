//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SecantSolve.h"

#include "FEProblem.h"
#include "Executioner.h"
#include "MooseMesh.h"
#include "NonlinearSystem.h"
#include "AllLocalDofIndicesThread.h"
#include "Console.h"
#include "EigenExecutionerBase.h"

defineLegacyParams(SecantSolve);

InputParameters
SecantSolve::validParams()
{
  InputParameters params = emptyInputParameters();

  params.addParam<unsigned int>(
      "secant_min_its",
      1,
      "Specifies the minimum number of Secant method iterations.");
  params.addParam<unsigned int>(
      "secant_max_its",
      1,
      "Specifies the maximum number of Secant method iterations.  Mainly used when "
      "wanting to do Secant method iterations with MultiApps that are set to "
      "execute_on timestep_end or timestep_begin. Setting this parameter to 1 turns off the Secant method "
      "iterations.");
  params.addParam<bool>(
      "accept_on_max_secant_iteration",
      false,
      "True to treat reaching the maximum number of Secant method iterations as converged.");
  params.addParam<bool>("disable_secant_residual_norm_check",
                        false,
                        "Disable the Secant method residual norm evaluation thus the three parameters "
                        "secant_rel_tol, secant_abs_tol and secant_force_norms.");
  params.addParam<Real>("secant_rel_tol",
                        1e-8,
                        "The relative nonlinear residual drop to shoot for "
                        "during Secant method iterations. This check is "
                        "performed based on the Master app's nonlinear "
                        "residual.");
  params.addParam<Real>("secant_abs_tol",
                        1e-50,
                        "The absolute nonlinear residual to shoot for "
                        "during Secant method iterations. This check is "
                        "performed based on the Master app's nonlinear "
                        "residual.");
  params.addParam<PostprocessorName>("secant_custom_pp",
                                     "Postprocessor for custom secant convergence check.");
  params.addParam<Real>("custom_rel_tol",
                        1e-8,
                        "The relative nonlinear residual drop to shoot for "
                        "during Secant method iterations. This check is "
                        "performed based on the postprocessor defined by "
                        "secant_custom_pp residual.");
  params.addParam<Real>("custom_abs_tol",
                        1e-50,
                        "The absolute nonlinear residual to shoot for "
                        "during Secant method iterations. This check is "
                        "performed based on postprocessor defined by "
                        "secant_custom_pp residual.");
  params.addParam<bool>("direct_pp_value",
                        false,
                        "True to use direct postprocessor value "
                        "(scaled by value on first iteration). "
                        "False (default) to use difference in postprocessor "
                        "value between secant iterations.");
  params.addParam<bool>(
      "secant_force_norms",
      false,
      "Force the evaluation of both the TIMESTEP_BEGIN and TIMESTEP_END norms regardless of the "
      "existence of active MultiApps with those execute_on flags, default: false.");

  params.addRangeCheckedParam<Real>("max_change_fraction",
                                    std::numeric_limits<Real>::max(),
                                    "max_change_fraction>=0",
                                    "Maximum relative change in the relaxed parameters between iterations.");
  params.addRangeCheckedParam<Real>("relaxation_factor",
                                    1.0,
                                    "relaxation_factor>0 & relaxation_factor<2",
                                    "Fraction of newly computed value to keep."
                                    "Set between 0 and 2.");
  params.addParam<std::vector<std::string>>("relaxed_variables",
                                            std::vector<std::string>(),
                                            "List of master app variables to relax during Secant method Iteration");
  params.addParam<std::vector<std::string>>("relaxed_postprocessors",
                                            std::vector<std::string>(),
                                            "List of master app postprocessors to relax during Secant method Iteration");

  params.addParamNamesToGroup("secant_min_its secant_max_its accept_on_max_secant_iteration "
                              "disable_secant_residual_norm_check secant_rel_tol "
                              "secant_abs_tol secant_custom_pp secant_force_norms "
                              "relaxation_factor relaxed_variables relaxed_postprocessors",
                              "Secant method");

  params.addParam<unsigned int>(
      "max_xfem_update",
      std::numeric_limits<unsigned int>::max(),
      "Maximum number of times to update XFEM crack topology in a step due to evolving cracks");
  params.addParam<bool>("update_xfem_at_timestep_begin",
                        false,
                        "Should XFEM update the mesh at the beginning of the timestep");
  params.addParam<bool>("auto_advance",
                        "Whether to automatically advance sub-applications regardless of whether "
                        "their solve converges.");

  return params;
}

SecantSolve::SecantSolve(Executioner * ex)
  : SolveObject(ex),
    _secant_min_its(getParam<unsigned int>("secant_min_its")),
    _secant_max_its(getParam<unsigned int>("secant_max_its")),
    _has_secant_its(_secant_max_its > 1),
    _accept_max_it(getParam<bool>("accept_on_max_secant_iteration")),
    _has_secant_norm(!getParam<bool>("disable_secant_residual_norm_check")),
    _secant_rel_tol(getParam<Real>("secant_rel_tol")),
    _secant_abs_tol(getParam<Real>("secant_abs_tol")),
    _secant_custom_pp(isParamValid("secant_custom_pp") ? &getPostprocessorValue("secant_custom_pp")
                                                       : nullptr),
    _custom_rel_tol(getParam<Real>("custom_rel_tol")),
    _custom_abs_tol(getParam<Real>("custom_abs_tol")),
    _secant_force_norms(getParam<bool>("secant_force_norms")),
    _relax_factor(getParam<Real>("relaxation_factor")),
    _max_change_fraction(getParam<Real>("max_change_fraction")),
    _relaxed_vars(getParam<std::vector<std::string>>("relaxed_variables")),
    _relaxed_pps(getParam<std::vector<std::string>>("relaxed_postprocessors")),
    // this value will be set by MultiApp
    _secant_self_relaxation_factor(1.0),
    _max_xfem_update(getParam<unsigned int>("max_xfem_update")),
    _update_xfem_at_timestep_begin(getParam<bool>("update_xfem_at_timestep_begin")),
    _secant_timer(registerTimedSection("SecantSolve", 1)),
    _secant_it(0),
    _secant_status(MooseSecantConvergenceReason::UNSOLVED),
    _xfem_update_count(0),
    _xfem_repeat_step(false),
    _previous_entering_time(_problem.time() - 1),
    _solve_message(_problem.shouldSolve() ? "Solve Converged!" : "Solve Skipped!"),
    _auto_advance_set_by_user(isParamValid("auto_advance")),
    _auto_advance_user_value(_auto_advance_set_by_user ? getParam<bool>("auto_advance") : true),
    _fail_step(false)
{
  if (_secant_min_its > _secant_max_its)
    paramError("secant_min_its", "The minimum number of Secant method iterations may not exceed the maximum.");

  if (_relax_factor != 1.0)
  {
    // Store a copy of the previous solution here
    _problem.getNonlinearSystemBase().addVector("relax_previous", false, PARALLEL);
    // Store a copy of the solution before the previous solution here
    _problem.getNonlinearSystemBase().addVector("relax_before", false, PARALLEL);
  }

  if (_relaxed_vars.size() > 0 && _relaxed_pps.size() > 0)
    mooseWarning("Both variable and postprocessor relaxation are active. If the two share dofs, the "
                 "relaxation will not be correct.");

  _old_relaxed_pps_values.resize(_relaxed_pps.size());
  _older_relaxed_pps_values.resize(_relaxed_pps.size());
}

bool
SecantSolve::solve()
{
  TIME_SECTION(_secant_timer);

  Real current_dt = _problem.dt();

  _secant_timestep_begin_norm.clear();
  _secant_timestep_end_norm.clear();
  _secant_timestep_begin_norm.resize(_secant_max_its);
  _secant_timestep_end_norm.resize(_secant_max_its);

  bool converged = true;

  // need to back up multi-apps even when not doing Secant method iteration for recovering from failed
  // multiapp solve
  _problem.backupMultiApps(EXEC_TIMESTEP_BEGIN);
  _problem.backupMultiApps(EXEC_TIMESTEP_END);

  // Prepare to relax variables as a master
  std::set<dof_id_type> relaxed_dofs;
  if (_relax_factor != 1.0)
  {
    // Snag all of the local dof indices for all of these variables
    System & libmesh_nl_system = _nl.system();
    AllLocalDofIndicesThread aldit(libmesh_nl_system, _relaxed_vars);
    ConstElemRange & elem_range = *_problem.mesh().getActiveLocalElementRange();
    Threads::parallel_reduce(elem_range, aldit);

    relaxed_dofs = aldit.getDofIndices();
  }

  // Prepare to relax variables as a subapp
  std::set<dof_id_type> self_relaxed_dofs;
  if (_secant_self_relaxation_factor != 1.0)
  {
    // Snag all of the local dof indices for all of these variables
    System & libmesh_nl_system = _nl.system();
    AllLocalDofIndicesThread aldit(libmesh_nl_system, _secant_self_relaxed_variables);
    ConstElemRange & elem_range = *_problem.mesh().getActiveLocalElementRange();
    Threads::parallel_reduce(elem_range, aldit);

    self_relaxed_dofs = aldit.getDofIndices();

    if (_previous_entering_time == _problem.time())
    {
      NumericVector<Number> & solution = _nl.solution();
      NumericVector<Number> & relax_previous = _nl.getVector("self_relax_previous");
      NumericVector<Number> & relax_before = _nl.getVector("self_relax_before");
      relax_before = relax_previous;
      relax_previous = solution;
    }
  }

  Real pp_scaling = 1.0;
  std::ostringstream pp_history;

  for (_secant_it = 0; _secant_it < _secant_max_its; ++_secant_it)  //FIXME use half its
  {
    if (_has_secant_its)
    {
      if (_secant_it == 0)
      {
        if (_has_secant_norm)
        {
          // First Secant method iteration - need to save off the initial nonlinear residual
          _secant_initial_norm = _problem.computeResidualL2Norm();
          _console << COLOR_MAGENTA << "Initial Secant method Norm: " << COLOR_DEFAULT;
          if (_secant_initial_norm == std::numeric_limits<Real>::max())
            _console << " MAX ";
          else
            _console << std::scientific << _secant_initial_norm;
          _console << COLOR_DEFAULT << "\n\n";
        }
      }
      else
      {
        // For every iteration other than the first, we need to restore the state of the MultiApps
        _problem.restoreMultiApps(EXEC_TIMESTEP_BEGIN);
        _problem.restoreMultiApps(EXEC_TIMESTEP_END);
      }

      _console << COLOR_MAGENTA << "Beginning Secant method Iteration " << _secant_it << COLOR_DEFAULT
               << '\n';
    }

    // Save last postprocessor value as value before solve
    Real pp_old = 0.0;
    if (_secant_custom_pp && _secant_it > 0 && !getParam<bool>("direct_pp_value"))
      pp_old = *_secant_custom_pp;

    // FIXME Better to relax postprocessors here??
    // Other potential location for relaxation of PPs
    // Relax postprocessors for the main application
    std::cout << "Secant method iteration: " << _secant_it << std::endl;
    for (size_t i=0; i<_relaxed_pps.size(); i++)
    {
      if (_secant_it % 2 == 0)
      {
        // Get new postprocessor value
        const Real current_value = getPostprocessorValueByName(_relaxed_pps[i]);
        const Real old_value = _old_relaxed_pps_values[i];
        const Real older_value = _older_relaxed_pps_values[i];

        // Compute and set relaxed value
        Real new_value;
        if (_secant_it > 0)   //FIXME
          new_value = older_value - (old_value - older_value)*(old_value - older_value) / (current_value + older_value - 2*old_value);  // steffensen
          // new_value = old_value - (current_value - old_value) * (old_value - older_value) / (current_value - 2*old_value + older_value);  // secant
        else
          new_value = current_value;

        // Bound new value change
        if (std::abs((new_value - current_value) / current_value) < 1 - _max_change_fraction)
          new_value = (1 - _max_change_fraction) * current_value;
        else if (std::abs((new_value - current_value) / current_value) < 1 + _max_change_fraction)
          new_value = (1 + _max_change_fraction) * current_value;

        _problem.setPostprocessorValueByName(_relaxed_pps[i], new_value);

        // Print new value
        std::cout << _relaxed_pps[i] << " " << current_value << " & " << old_value << " & " << older_value << " -> " << new_value << std::endl;
      }
      // Keep track of postprocessor values for the next iteration
      _older_relaxed_pps_values[i] = _old_relaxed_pps_values[i];
      _old_relaxed_pps_values[i] = getPostprocessorValueByName(_relaxed_pps[i]);
    }

    Real begin_norm_old = (_secant_it > 0 ? _secant_timestep_begin_norm[_secant_it - 1]
                                          : std::numeric_limits<Real>::max());
    Real end_norm_old = (_secant_it > 0 ? _secant_timestep_end_norm[_secant_it - 1]
                                        : std::numeric_limits<Real>::max());
    bool relax = (_relax_factor != 1) && (_secant_it > 0) && (_secant_it % 2 == 0);
    bool solve_converged = solveStep(begin_norm_old,
                                     _secant_timestep_begin_norm[_secant_it],
                                     end_norm_old,
                                     _secant_timestep_end_norm[_secant_it],
                                     relax,
                                     relaxed_dofs);

    // Calculate error from postprocessor values
    Real pp_new = std::numeric_limits<Real>::max();
    if (_secant_custom_pp)
    {
      if ((_secant_it == 0 && getParam<bool>("direct_pp_value")) ||
          !getParam<bool>("direct_pp_value"))
        pp_scaling = *_secant_custom_pp;
      pp_new = *_secant_custom_pp;

      auto ppname = getParam<PostprocessorName>("secant_custom_pp");
      pp_history << std::setw(2) << _secant_it + 1 << " Secant method " << ppname << " = "
                 << Console::outputNorm(std::numeric_limits<Real>::max(), pp_new) << "\n";
      _console << pp_history.str();
    }

    if (solve_converged)
    {
      if (_has_secant_its)
      {
        if (_has_secant_norm)
        {
          _console << "\n 0 Secant method    |R| = "
                   << Console::outputNorm(std::numeric_limits<Real>::max(), _secant_initial_norm)
                   << '\n';

          for (unsigned int i = 0; i <= _secant_it; ++i)
          {
            Real max_norm = std::max(_secant_timestep_begin_norm[i], _secant_timestep_end_norm[i]);
            std::stringstream secant_prefix;
            if (i % 2 == 1)
              secant_prefix << " Secant half-step |R| = ";
            else
              secant_prefix << " Secant step      |R| = ";

            _console << std::setw(2) << i + 1
                     << secant_prefix.str() << Console::outputNorm(_secant_initial_norm, max_norm)
                     << '\n';
          }

          if  (_secant_it + 2 > 2 * _secant_min_its)
          {
            Real max_norm = std::max(_secant_timestep_begin_norm[_secant_it],
                                     _secant_timestep_end_norm[_secant_it]);

            Real max_relative_drop = max_norm / _secant_initial_norm;

            if (max_norm < _secant_abs_tol)
            {
              _secant_status = MooseSecantConvergenceReason::CONVERGED_ABS;
              break;
            }
            if (max_relative_drop < _secant_rel_tol)
            {
              _secant_status = MooseSecantConvergenceReason::CONVERGED_RELATIVE;
              break;
            }
          }
          if (_executioner.augmentedPicardConvergenceCheck())   //// FIXME
          {
            _secant_status = MooseSecantConvergenceReason::CONVERGED_CUSTOM;
            break;
          }
          if (std::abs(pp_new - pp_old) < _custom_abs_tol)
          {
            _secant_status = MooseSecantConvergenceReason::CONVERGED_CUSTOM;
            break;
          }
          if (std::abs((pp_new - pp_old) / pp_scaling) < _custom_rel_tol)
          {
            _secant_status = MooseSecantConvergenceReason::CONVERGED_CUSTOM;
            break;
          }
          if (_secant_it + 1 == _secant_max_its)
          {
            if (_accept_max_it)
            {
              _secant_status = MooseSecantConvergenceReason::REACH_MAX_ITS;
              converged = true;
            }
            else
            {
              _secant_status = MooseSecantConvergenceReason::DIVERGED_MAX_ITS;
              converged = false;
            }
            break;
          }
        }
      }
    }
    else
    {
      // If the last solve didn't converge then we need to exit this step completely (even in the
      // case of Secant method). So we can retry...
      converged = false;
      break;
    }

    _problem.dt() =
        current_dt; // _dt might be smaller than this at this point for multistep methods
  }

  if (converged && _secant_self_relaxation_factor != 1.0)
  {
    std::cout << "RELAXING POST SECANT" << std::endl;
    if (_previous_entering_time == _problem.time())
    {
      NumericVector<Number> & solution = _nl.solution();
      NumericVector<Number> & relax_previous = _nl.getVector("self_relax_previous");
      NumericVector<Number> & relax_before = _nl.getVector("self_relax_before");

      const Real factor = _secant_self_relaxation_factor;
      for (const auto & dof : self_relaxed_dofs)\
      {
        std::cout << "dof " << dof << " " << relax_before(dof) << " & " << relax_previous(dof) << " & " << solution(dof) << " : " <<
          relax_before(dof) - (relax_previous(dof) - relax_before(dof)) * (relax_previous(dof) - relax_before(dof))
             / (solution(dof) + relax_before(dof) - 2*relax_previous(dof)) << std::endl;
        solution.set(dof, relax_before(dof) - (relax_previous(dof) - relax_before(dof)) * (relax_previous(dof) - relax_before(dof))
             / (solution(dof) + relax_before(dof) - 2*relax_previous(dof)));
      }
      solution.close();
      _nl.update();


      // Relax the postprocessors
      std::cout << "Sub-Secant method iteration: " << _secant_it << std::endl;
      for (size_t i=0; i<_secant_self_relaxed_pps.size(); i++)
      {
        // Get new postprocessor value
        const Real current_value = getPostprocessorValueByName(_secant_self_relaxed_pps[i]);
        const Real old_value = _problem.getPostprocessorValueByName(_secant_self_relaxed_pps[i], 1);  //FIXME
        const Real older_value = _problem.getPostprocessorValueByName(_secant_self_relaxed_pps[i], 1);  //FIXME

        // Compute and set relaxed value
        Real new_value = current_value;
        const Real factor = _secant_self_relaxation_factor;
        new_value = factor * current_value + (1 - factor) * old_value;   //FIXME
        _problem.setPostprocessorValueByName(_secant_self_relaxed_pps[i], new_value);
        _problem.setPostprocessorValueByName(_secant_self_relaxed_pps[i], new_value, 1);
        // _problem.setPostprocessorValueByName(_secant_self_relaxed_pps[i], new_value, 2);

        std::cout << _secant_self_relaxed_pps[i] << " " << current_value << " & " << old_value  << " " << _problem.getPostprocessorValueByName(_secant_self_relaxed_pps[i], 1) << " -> " << new_value << std::endl;
      }
    }
    _previous_entering_time = _problem.time();
  }

  if (_has_secant_its)
  {
    _console << "Secant method converged reason: ";
    switch (_secant_status)
    {
      case MooseSecantConvergenceReason::CONVERGED_ABS:
        _console << "CONVERGED_ABS";
        break;
      case MooseSecantConvergenceReason::CONVERGED_RELATIVE:
        _console << "CONVERGED_RELATIVE";
        break;
      case MooseSecantConvergenceReason::CONVERGED_CUSTOM:
        _console << "CONVERGED_CUSTOM";
        break;
      case MooseSecantConvergenceReason::REACH_MAX_ITS:
        _console << "REACH_MAX_ITS";
        break;
      case MooseSecantConvergenceReason::DIVERGED_MAX_ITS:
        _console << "DIVERGED_MAX_ITS";
        break;
      case MooseSecantConvergenceReason::DIVERGED_NONLINEAR:
        _console << "DIVERGED_NONLINEAR";
        break;
      case MooseSecantConvergenceReason::DIVERGED_FAILED_MULTIAPP:
        _console << "DIVERGED_FAILED_MULTIAPP";
        break;
      default:
        // UNSOLVED and CONVERGED_NONLINEAR should not be hit when Secant method
        // iteration is not on here
        mooseError("Internal error: wrong Secant method status!");
        break;
    }
    _console << std::endl;
  }
  return converged;
}

bool
SecantSolve::autoAdvance() const
{
  bool auto_advance = !(_has_secant_its && _problem.isTransient());

  if (dynamic_cast<EigenExecutionerBase *>(&_executioner) && _has_secant_its)
    auto_advance = true;

  if (_auto_advance_set_by_user)
    auto_advance = _auto_advance_user_value;

  return auto_advance;
}

bool
SecantSolve::solveStep(Real begin_norm_old,
                       Real & begin_norm,
                       Real end_norm_old,
                       Real & end_norm,
                       bool relax,
                       const std::set<dof_id_type> & relaxed_dofs)
{
  bool auto_advance = autoAdvance();

  _executioner.preSolve();

  _problem.execTransfers(EXEC_TIMESTEP_BEGIN);
  if (!_problem.execMultiApps(EXEC_TIMESTEP_BEGIN, auto_advance))
  {
    _secant_status = MooseSecantConvergenceReason::DIVERGED_FAILED_MULTIAPP;
    return false;
  }

  if (_problem.haveXFEM() && _update_xfem_at_timestep_begin)
    _problem.updateMeshXFEM();

  _problem.execute(EXEC_TIMESTEP_BEGIN);

  if (_has_secant_its && _has_secant_norm)
    if (_problem.hasMultiApps(EXEC_TIMESTEP_BEGIN) || _secant_force_norms)
    {
      begin_norm = _problem.computeResidualL2Norm();

      _console << COLOR_MAGENTA << "Secant method Norm after TIMESTEP_BEGIN MultiApps: "
               << Console::outputNorm(begin_norm_old, begin_norm) << '\n';
    }

  // Perform output for timestep begin
  _problem.outputStep(EXEC_TIMESTEP_BEGIN);

  // Update warehouse active objects
  _problem.updateActiveObjects();

  if (relax)
  {
    NumericVector<Number> & solution = _nl.solution();
    NumericVector<Number> & relax_previous = _nl.getVector("relax_previous");
    NumericVector<Number> & relax_before = _nl.getVector("relax_before");

    // Save off the current and previous solutions
    relax_before = relax_previous;
    relax_previous = solution;
  }

  if (_has_secant_its)
    _console << COLOR_MAGENTA << "\nMaster solve:\n" << COLOR_DEFAULT;
  if (!_inner_solve->solve())
  {
    _secant_status = MooseSecantConvergenceReason::DIVERGED_NONLINEAR;

    _console << COLOR_RED << " Solve Did NOT Converge!" << COLOR_DEFAULT << std::endl;
    // Perform the output of the current, failed time step (this only occurs if desired)
    _problem.outputStep(EXEC_FAILED);
    return false;
  }
  else
    _secant_status = MooseSecantConvergenceReason::CONVERGED_NONLINEAR;

  _console << COLOR_GREEN << ' ' << _solve_message << COLOR_DEFAULT << std::endl;

  // Relax the "relaxed_variables" and "relaxed_postprocessors"
  if (relax)
  {
    std::cout << "Relaxing IN STEP " << std::endl;
    NumericVector<Number> & solution = _nl.solution();
    NumericVector<Number> & relax_previous = _nl.getVector("relax_previous");
    NumericVector<Number> & relax_before = _nl.getVector("relax_before");

    for (const auto & dof : relaxed_dofs)
    {
      // std::cout << "dof " << dof << " " << relax_before(dof) << " & " << relax_previous(dof) << " & " << solution(dof) << " : " <<
      //   relax_before(dof) - (relax_previous(dof) - relax_before(dof)) * (relax_previous(dof) - relax_before(dof))
      //      / (solution(dof) + relax_before(dof) - 2*relax_previous(dof)) << std::endl;
      if (std::abs(solution(dof) + relax_before(dof) - 2*relax_previous(dof)) > 1e-16)
        solution.set(dof, relax_before(dof) - (relax_previous(dof) - relax_before(dof)) * (relax_previous(dof) - relax_before(dof))
            / (solution(dof) + relax_before(dof) - 2*relax_previous(dof)));
      else
        solution.set(dof, solution(dof));
    }
    solution.close();
    _nl.update();
  }

  if (_problem.haveXFEM() && (_xfem_update_count < _max_xfem_update) && _problem.updateMeshXFEM())
  {
    _console << "\nXFEM modified mesh, repeating step" << std::endl;
    _xfem_repeat_step = true;
    ++_xfem_update_count;
  }
  else
  {
    if (_problem.haveXFEM())
    {
      _xfem_repeat_step = false;
      _xfem_update_count = 0;
      _console << "\nXFEM did not modify mesh, continuing" << std::endl;
    }

    _problem.onTimestepEnd();
    _problem.execute(EXEC_TIMESTEP_END);

    _problem.execTransfers(EXEC_TIMESTEP_END);
    if (!_problem.execMultiApps(EXEC_TIMESTEP_END, auto_advance))
    {
      _secant_status = MooseSecantConvergenceReason::DIVERGED_FAILED_MULTIAPP;
      return false;
    }
  }

  if (_fail_step)
  {
    _fail_step = false;
    return false;
  }

  _executioner.postSolve();

  if (_has_secant_its && _has_secant_norm)
    if (_problem.hasMultiApps(EXEC_TIMESTEP_END) || _secant_force_norms)
    {
      end_norm = _problem.computeResidualL2Norm();

      _console << COLOR_MAGENTA << "Secant method Norm after TIMESTEP_END MultiApps: "
               << Console::outputNorm(end_norm_old, end_norm) << '\n';
    }

  return true;
}
