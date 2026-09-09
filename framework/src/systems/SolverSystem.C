//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SolverSystem.h"
#include "SolutionInvalidity.h"
#include "FEProblemBase.h"
#include "TimeIntegrator.h"
#include "MooseUtils.h"

SolverSystem::SolverSystem(SubProblem & subproblem,
                           FEProblemBase & fe_problem,
                           const std::string & name,
                           Moose::VarKindType var_kind)
  : SystemBase(subproblem, fe_problem, name, var_kind),
    _current_solution(nullptr),
    _pc_side(Moose::PCS_DEFAULT),
    _ksp_norm(Moose::KSPN_UNPRECONDITIONED)
{
}

SolverSystem::~SolverSystem() = default;

void
SolverSystem::preInit()
{
  SystemBase::preInit();

  _current_solution = system().current_local_solution.get();

  if (_serialized_solution.get())
    _serialized_solution->init(system().n_dofs(), false, SERIAL);
}

void
SolverSystem::restoreSolutions()
{
  // call parent
  SystemBase::restoreSolutions();
  // and update _current_solution
  _current_solution = system().current_local_solution.get();
}

void
SolverSystem::serializeSolution()
{
  if (_serialized_solution.get())
  {
    if (!_serialized_solution->initialized() || _serialized_solution->size() != system().n_dofs())
    {
      _serialized_solution->clear();
      _serialized_solution->init(system().n_dofs(), false, SERIAL);
    }

    _current_solution->localize(*_serialized_solution);
  }
}

void
SolverSystem::setSolution(const NumericVector<Number> & soln)
{
  _current_solution = &soln;

  auto tag = _subproblem.getVectorTagID(Moose::SOLUTION_TAG);
  associateVectorToTag(const_cast<NumericVector<Number> &>(soln), tag);

  if (_serialized_solution.get())
    serializeSolution();
}

void
SolverSystem::applyFixedPointRelaxation(const Real relaxation_factor,
                                        const Moose::SolutionIterationType iteration_type)
{
  if (MooseUtils::absoluteFuzzyEqual(relaxation_factor, 1.0))
    return;

  mooseAssert(hasSolutionState(1, iteration_type),
              "Fixed point relaxation was requested but the old fixed point solution was not "
              "saved.");

  // This might be paranoid but who knows, maybe someone requests nonghosted
  mooseAssert(solutionStateParallelType(1, iteration_type) == solution().type(),
              "Fixed point relaxation requires the previous fixed point solution state to have "
              "the same parallel type as the system solution.");

  auto & sol = solution();
  sol.scale(relaxation_factor);
  sol.add(1.0 - relaxation_factor, solutionState(1, iteration_type));
  sol.close();
  update();
}

void
SolverSystem::secantSlopeContribution(const TagID prev_eval_tag,
                                      const Moose::SolutionIterationType iteration_type,
                                      Real & dx_dot_dg,
                                      Real & dg_dot_dg)
{
  mooseAssert(hasVector(prev_eval_tag),
              "Secant fixed point update requested but the previous-evaluation vector was not "
              "allocated.");
  mooseAssert(hasSolutionState(2, iteration_type),
              "Secant fixed point update requires two previous solution states.");

  auto & sol = solution();                     // f(x_n)
  auto & prev_eval = getVector(prev_eval_tag); // f(x_{n-1})
  const auto & xn = solutionState(1, iteration_type);
  const auto & xnm1 = solutionState(2, iteration_type);

  // Fixed point residuals g(x) = f(x) - x at the two most recent iterates
  auto g_n = sol.clone();
  g_n->add(-1.0, xn); // g(x_n) = f(x_n) - x_n
  auto g_nm1 = prev_eval.clone();
  g_nm1->add(-1.0, xnm1); // g(x_{n-1}) = f(x_{n-1}) - x_{n-1}

  // Iterate and residual changes over the last sweep
  auto dg = g_n->clone();
  dg->add(-1.0, *g_nm1); // dg = g(x_n) - g(x_{n-1})
  auto dx = xn.clone();
  dx->add(-1.0, xnm1); // dx = x_n - x_{n-1}

  dx_dot_dg += dx->dot(*dg);
  dg_dot_dg += dg->dot(*dg);
}

void
SolverSystem::applyFixedPointSecant(const bool use_secant,
                                    const Real mu,
                                    const Real relaxation_factor,
                                    const TagID prev_eval_tag,
                                    const Moose::SolutionIterationType iteration_type)
{
  mooseAssert(hasVector(prev_eval_tag),
              "Secant fixed point update requested but the previous-evaluation vector was not "
              "allocated.");

  // f(x_n): the raw system solve output at the current iterate
  auto & sol = solution();
  // f(x_{n-1}): the raw system solve output recorded during the previous iteration
  auto & prev_eval = getVector(prev_eval_tag);

  // Without a valid secant step size (first sweep, or a vanishing coupled residual change), record
  // the evaluation and take the usual (optionally relaxed) Picard step.
  if (!use_secant)
  {
    prev_eval = sol;
    prev_eval.close();
    applyFixedPointRelaxation(relaxation_factor, iteration_type);
    return;
  }

  mooseAssert(solutionStateParallelType(1, iteration_type) == sol.type(),
              "Secant fixed point update requires the previous solution state to have the same "
              "parallel type as the system solution.");

  const auto & xn = solutionState(1, iteration_type); // x_n produced f(x_n)

  // g(x_n) = f(x_n) - x_n, computed before overwriting the solution
  auto g_n = sol.clone();
  g_n->add(-1.0, xn);

  // Record f(x_n) for the next iteration's secant slope
  prev_eval = sol;
  prev_eval.close();

  // x_{n+1} = x_n - mu * g(x_n), then under/over-relaxed by relaxation_factor
  sol = xn;
  sol.add(-mu, *g_n);
  if (!MooseUtils::absoluteFuzzyEqual(relaxation_factor, 1.0))
  {
    sol.scale(relaxation_factor);
    sol.add(1.0 - relaxation_factor, xn);
  }
  sol.close();
  update();
}

void
SolverSystem::setPCSide(MooseEnum pcs)
{
  if (pcs == "left")
    _pc_side = Moose::PCS_LEFT;
  else if (pcs == "right")
    _pc_side = Moose::PCS_RIGHT;
  else if (pcs == "symmetric")
    _pc_side = Moose::PCS_SYMMETRIC;
  else if (pcs == "default")
    _pc_side = Moose::PCS_DEFAULT;
  else
    mooseError("Unknown PC side specified.");
}

void
SolverSystem::setMooseKSPNormType(MooseEnum kspnorm)
{
  if (kspnorm == "none")
    _ksp_norm = Moose::KSPN_NONE;
  else if (kspnorm == "preconditioned")
    _ksp_norm = Moose::KSPN_PRECONDITIONED;
  else if (kspnorm == "unpreconditioned")
    _ksp_norm = Moose::KSPN_UNPRECONDITIONED;
  else if (kspnorm == "natural")
    _ksp_norm = Moose::KSPN_NATURAL;
  else if (kspnorm == "default")
    _ksp_norm = Moose::KSPN_DEFAULT;
  else
    mooseError("Unknown ksp norm type specified.");
}

void
SolverSystem::checkInvalidSolution()
{
  auto & solution_invalidity = _app.solutionInvalidity();

  // sync all solution invalid counts to rank 0 process
  solution_invalidity.syncIteration();

  if (solution_invalidity.hasInvalidSolution())
  {
    if (_fe_problem.acceptInvalidSolution())
      if (_fe_problem.showInvalidSolutionConsole())
        solution_invalidity.print(_console);
      else
        mooseWarning("The Solution Invalidity warnings are detected but silenced! "
                     "Use Problem/show_invalid_solution_console=true to show solution counts");
    else
      // output the occurrence of solution invalid in a summary table
      if (_fe_problem.showInvalidSolutionConsole())
        solution_invalidity.print(_console);
  }
}

void
SolverSystem::compute(const ExecFlagType type)
{
  // Let's try not to overcompute
  bool compute_tds = false;
  if (type == EXEC_LINEAR)
    compute_tds = true;
  else if (type == EXEC_NONLINEAR)
  {
    if (_fe_problem.computingScalingJacobian() || matrixFromColoring())
      compute_tds = true;
  }
  else if ((type == EXEC_TIMESTEP_END) || (type == EXEC_FINAL))
  {
    if (_fe_problem.solverParams(number())._type == Moose::ST_LINEAR)
      // We likely don't have a final residual evaluation upon which we compute the time derivatives
      // so we need to do so now
      compute_tds = true;
  }

  // avoid division by dt which might be zero.
  if (compute_tds && _fe_problem.dt() > 0.)
    for (auto & ti : _time_integrators)
    {
      // Do things like compute integration weights
      ti->preStep();
      ti->computeTimeDerivatives();
    }
}
