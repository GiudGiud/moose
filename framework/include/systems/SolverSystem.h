//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "SystemBase.h"
#include "MooseTypes.h"

#include "libmesh/system.h"

#include <string>

class SubProblem;
class FEProblemBase;

namespace SparsityPattern = libMesh::SparsityPattern;

class SolverSystem : public SystemBase
{
public:
  SolverSystem(SubProblem & subproblem,
               FEProblemBase & fe_problem,
               const std::string & name,
               Moose::VarKindType var_kind);
  virtual ~SolverSystem();

  virtual void preInit() override;
  virtual void restoreSolutions() override final;

  void serializeSolution();

  /**
   * Quit the current solve as soon as possible.
   */
  virtual void stopSolve(const ExecFlagType & exec_flag,
                         const std::set<TagID> & vector_tags_to_close) = 0;

  /**
   * Returns the convergence state
   * @return true if converged, otherwise false
   */
  virtual bool converged() = 0;

  /**
   * If the system has a kernel that corresponds to a time derivative
   */
  virtual bool containsTimeKernel() = 0;

  /**
   * Returns the names of the variables that have time derivative kernels
   * in the system.
   */
  virtual std::vector<std::string> timeKernelVariableNames() = 0;

  /**
   * Set the solution to a given vector.
   * @param soln The vector which should be treated as the solution.
   */
  void setSolution(const NumericVector<Number> & soln);

  /**
   * Apply solution under/over-relaxation for fixed point iterations.
   *
   * The relaxed update is:
   *   u <- relaxation_factor * u_new + (1 - relaxation_factor) * u_old
   *
   * @param[in] relaxation_factor  The factor applied to the new solution
   * @param[in] iteration_type  Type of iteration; which "previous" value to use
   */
  void applyFixedPointRelaxation(const Real relaxation_factor,
                                 const Moose::SolutionIterationType iteration_type);

  /**
   * Accumulate this system's contribution to a projected (Aitken) secant step size for a coupled
   * multi-system fixed point iteration.
   *
   * Using the fixed point residual g(x) = f(x) - x, this adds the inner products
   *   dx . dg   and   dg . dg,   with dx = x_n - x_{n-1},  dg = g(x_n) - g(x_{n-1})
   * to the running totals. The caller forms a single step size mu = sum(dx.dg) / sum(dg.dg) shared
   * by all systems, so the secant slope captures the coupling between systems rather than treating
   * each system as a self-map. Requires two previous solution states and the recorded evaluation
   * f(x_{n-1}); does not modify the solution.
   *
   * @param[in] prev_eval_tag  Tag of the vector storing the previous evaluation f(x_{n-1})
   * @param[in] iteration_type  Type of iteration; which "previous" solution states to use
   * @param[in,out] dx_dot_dg  Running total of dx . dg, incremented by this system's contribution
   * @param[in,out] dg_dot_dg  Running total of dg . dg, incremented by this system's contribution
   */
  void secantSlopeContribution(const TagID prev_eval_tag,
                               const Moose::SolutionIterationType iteration_type,
                               Real & dx_dot_dg,
                               Real & dg_dot_dg);

  /**
   * Apply a projected (Aitken) secant update to this system's solution using a step size mu
   * shared across the coupled systems:
   *   u <- x_n - mu * g(x_n),   g(x) = f(x) - x
   * which is then under/over-relaxed by relaxation_factor. The current evaluation f(x_n) is
   * recorded for the next iteration's secant slope. When use_secant is false (first iteration, or
   * a vanishing coupled residual change), no secant step is taken: the evaluation is recorded and
   * an (optionally relaxed) Picard step is applied instead.
   *
   * @param[in] use_secant  Whether a valid secant step size is available; if not, fall back to a
   *                        relaxed Picard step
   * @param[in] mu  The secant step size shared across the coupled systems
   * @param[in] relaxation_factor  The under/over-relaxation factor applied to the update
   * @param[in] prev_eval_tag  Tag of the vector storing the previous evaluation f(x_{n-1})
   * @param[in] iteration_type  Type of iteration; which "previous" solution states to use
   */
  void applyFixedPointSecant(const bool use_secant,
                             const Real mu,
                             const Real relaxation_factor,
                             const TagID prev_eval_tag,
                             const Moose::SolutionIterationType iteration_type);

  /**
   * Set the side on which the preconditioner is applied to.
   * @param pcs The required preconditioning side
   */
  void setPCSide(MooseEnum pcs);

  /**
   * Get the current preconditioner side.
   */
  Moose::PCSideType getPCSide() { return _pc_side; }

  /**
   * Set the norm in which the linear convergence will be measured.
   * @param kspnorm The required norm
   */
  void setMooseKSPNormType(MooseEnum kspnorm);

  /**
   * Get the norm in which the linear convergence is measured.
   */
  Moose::MooseKSPNormType getMooseKSPNormType() { return _ksp_norm; }

  virtual const NumericVector<Number> * const & currentSolution() const override final;

  virtual void compute(ExecFlagType type) override;

protected:
  void checkInvalidSolution();

  virtual NumericVector<Number> & solutionInternal() const override final;

  /**
   * Whether a system matrix is formed from coloring. This influences things like when to compute
   * time derivatives
   */
  virtual bool matrixFromColoring() const { return false; }

  /// solution vector from solver
  const NumericVector<Number> * _current_solution;

  /// Preconditioning side
  Moose::PCSideType _pc_side;
  /// KSP norm type
  Moose::MooseKSPNormType _ksp_norm;

  /// Boolean to see if solution is invalid
  bool _solution_is_invalid;
};

inline const NumericVector<Number> * const &
SolverSystem::currentSolution() const
{
  return _current_solution;
}

inline NumericVector<Number> &
SolverSystem::solutionInternal() const
{
  return *system().solution;
}
