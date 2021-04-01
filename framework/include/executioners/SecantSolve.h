//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SolveObject.h"

// System includes
#include <string>

class SecantSolve;

template <>
InputParameters validParams<SecantSolve>();

class SecantSolve : public SolveObject
{
public:
  SecantSolve(Executioner * ex);

  static InputParameters validParams();

  /**
   * Secant solve the FEProblem.
   * @return True if solver is converged.
   */
  virtual bool solve() override;

  /// Enumeration for Secant convergence reasons
  enum class MooseSecantConvergenceReason
  {
    UNSOLVED = 0,
    CONVERGED_NONLINEAR = 1,
    CONVERGED_ABS = 2,
    CONVERGED_RELATIVE = 3,
    CONVERGED_CUSTOM = 4,
    REACH_MAX_ITS = 5,
    DIVERGED_MAX_ITS = -1,
    DIVERGED_NONLINEAR = -2,
    DIVERGED_FAILED_MULTIAPP = -3
  };

  /**
   * Get the number of Secant iterations performed
   * Because this returns the number of Secant iterations, rather than the current
   * iteration count (which starts at 0), increment by 1.
   *
   * @return Number of Secant iterations performed
   */
  unsigned int numSecantIts() const { return _secant_it + 1; }

  /// Check the solver status
  MooseSecantConvergenceReason checkConvergence() const { return _secant_status; }

  /// This function checks the _xfem_repeat_step flag set by solve.
  bool XFEMRepeatStep() const { return _xfem_repeat_step; }

  /// Clear Secant status
  void clearSecantStatus() { _secant_status = MooseSecantConvergenceReason::UNSOLVED; }

  /// Whether or not this has Secant iterations
  bool hasSecantIteration() { return _has_secant_its; }

  /// Set relaxation factor for the current solve as a MultiApp
  void setMultiAppRelaxationFactor(Real factor) { _secant_self_relaxation_factor = factor; }

  /// Set relaxation variables for the current solve as a MultiApp
  void setMultiAppRelaxationVariables(const std::vector<std::string> & vars)
  {
    _secant_self_relaxed_variables = vars;
  }

  /// Set relaxation postprocessors for the current solve as a MultiApp
  void setMultiAppRelaxationPostprocessors(const std::vector<std::string> & pps)
  {
    _secant_self_relaxed_pps = pps;
  }

  /**
   * Whether sub-applications are automatically advanced no matter what happens during their solves
   */
  bool autoAdvance() const;

  /// mark the current solve as failed due to external conditions
  void failStep() { _fail_step = true; }

protected:
  /**
   * Perform one Secant iteration or a full solve.
   *
   * @param begin_norm_old Residual norm after timestep_begin execution of previous Secant
   * iteration
   * @param begin_norm     Residual norm after timestep_begin execution
   * @param end_norm_old   Residual norm after timestep_end execution of previous Secant iteration
   * @param end_norm       Residual norm after timestep_end execution
   * @param relax          Whether or not we do relaxation in this iteration
   * @param relaxed_dofs   DoFs to be relaxed
   *
   * @return True if both nonlinear solve and the execution of multiapps are successful.
   *
   * Note: this function also set _xfem_repeat_step flag for XFEM. It tracks _xfem_update_count
   * state.
   * FIXME: The proper design will be to let XFEM use Secant iteration to control the execution.
   */
  bool solveStep(Real begin_norm_old,
                 Real & begin_norm,
                 Real end_norm_old,
                 Real & end_norm,
                 bool relax,
                 const std::set<dof_id_type> & relaxed_dofs);

   /// Miniumum Secant iterations
   unsigned int _secant_min_its;
  /// Maximum Secant iterations
  unsigned int _secant_max_its;
  /// Whether or not we activate Secant iteration
  bool _has_secant_its;
  /// Whether or not to treat reaching maximum number of Secant iteration as converged
  bool _accept_max_it;
  /// Whether or not to use residual norm to check the Secant convergence
  bool _has_secant_norm;
  /// Relative tolerance on residual norm
  Real _secant_rel_tol;
  /// Absolute tolerance on residual norm
  Real _secant_abs_tol;
  /// Postprocessor value for user-defined secant convergence check
  const PostprocessorValue * const _secant_custom_pp;
  /// Relative tolerance on postprocessor value
  Real _custom_rel_tol;
  /// Absolute tolerance on postprocessor value
  Real _custom_abs_tol;
  /// Whether or not we force evaluation of residual norms even without multiapps
  bool _secant_force_norms;
  const Real _max_change_fraction;

  /// Relaxation factor for Secant Iteration
  const Real _relax_factor;
  /// The variables (transferred or not) that are going to be relaxed
  std::vector<std::string> _relaxed_vars;
  /// The postprocessors (transferred or not) that are going to be relaxed
  std::vector<std::string> _relaxed_pps;
  /// The relaxed postprocessors value at the previous iteration
  std::vector<PostprocessorValue> _old_relaxed_pps_values;
  /// The relaxed postprocessors value from two iterations before
  std::vector<PostprocessorValue> _older_relaxed_pps_values;

  /// Relaxation factor outside of Secant iteration (used as a subapp)
  Real _secant_self_relaxation_factor;
  /// Variables to be relaxed outside of Secant iteration (used as a subapp)
  std::vector<std::string> _secant_self_relaxed_variables;
  /// Postprocessors to be relaxed outside of Secant iteration (used as a subapp)
  std::vector<std::string> _secant_self_relaxed_pps;

  /// Maximum number of xfem updates per step
  unsigned int _max_xfem_update;
  /// Controls whether xfem should update the mesh at the beginning of the time step
  bool _update_xfem_at_timestep_begin;

private:
  /// Timer for Secant iteration
  const PerfID _secant_timer;

  ///@{ Variables used by the Secant iteration
  /// Secant iteration counter
  unsigned int _secant_it;
  /// Initial residual norm
  Real _secant_initial_norm;
  /// Full history of residual norm after evaluation of timestep_begin
  std::vector<Real> _secant_timestep_begin_norm;
  /// Full history of residual norm after evaluation of timestep_end
  std::vector<Real> _secant_timestep_end_norm;
  /// Status of Secant solve
  MooseSecantConvergenceReason _secant_status;
  ///@}

  /// Counter for number of xfem updates that have been performed in the current step
  unsigned int _xfem_update_count;
  /// Whether step should be repeated due to xfem modifying the mesh
  bool _xfem_repeat_step;

  /// Time of previous Secant solve as a subapp
  Real _previous_entering_time;

  const std::string _solve_message;

  /// Whether the user has set the auto_advance parameter for handling advancement of
  /// sub-applications in multi-app contexts
  const bool _auto_advance_set_by_user;

  /// The value of auto_advance set by the user for handling advancement of sub-applications in
  /// multi-app contexts
  const bool _auto_advance_user_value;

  /// force the current step to fail, triggering are repeat with a cut dt
  bool _fail_step;
};
