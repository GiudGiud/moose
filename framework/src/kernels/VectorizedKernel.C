//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Kernel.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "MooseVariableScalar.h"
#include "SubProblem.h"
#include "NonlinearSystem.h"

#include "libmesh/threads.h"
#include "libmesh/quadrature.h"

InputParameters
VectorizedKernel::validParams()
{
  InputParameters params = KernelBase::validParams();
  params.registerBase("Kernel");
  return params;
}

VectorizedKernel::VectorizedKernel(const InputParameters & parameters)
  : KernelBase(parameters),
    MooseVariableInterface<Real>(this,
                                 false,
                                 "variable",
                                 Moose::VarKindType::VAR_SOLVER,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _var(*mooseVariable()),
    _test(_var.phi()),
    _grad_test(_var.gradPhi()),
    _phi(_assembly.phi(_var)),
    _grad_phi(_assembly.gradPhi(_var)),
    _u(_is_implicit ? _var.sln() : _var.slnOld()),
    _grad_u(_is_implicit ? _var.gradSln() : _var.gradSlnOld())
{
  addMooseVariableDependency(mooseVariable());

  // TODO: initialize with max NQp
  _qp_residuals.resize(4);
  _qp_jacobians.resize(4);
}

void
VectorizedKernel::computeResidual()
{
  if (!hasVectorTags())
    return;

  prepareVectorTag(_assembly, _var.number());
  const unsigned int nqp = _qrule->n_points();

  precalculateResidual();
  for (_i = 0; _i < _test.size(); _i++)
  {
    computeQpResiduals();
    MooseArray<Real> qp_residuals(_nqp);
    qp_residuals = _qp_residuals;
    Real local_re = 0.;
    const auto * const A = _JxW.data();
    const auto * const B = _coord.data();
    const auto * const C = qp_residuals.data();

#pragma clang loop vectorize(enable) interleave(enable)
    for (unsigned int qp = 0; qp < 4; ++qp)
    {
      Real temp = A[qp]* B[qp];
      local_re +=  temp * C[qp];
    }
    _local_re(_i) += local_re;
  }

  accumulateTaggedLocalResidual();
}

void
VectorizedKernel::computeJacobian()
{
  prepareMatrixTag(_assembly, _var.number(), _var.number());
  _nqp = _qrule->n_points();

  precalculateJacobian();
  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < _phi.size(); _j++)
    {
      computeQpJacobians();
      Real local_ke = 0.;
      #pragma omp simd reduction(+: local_ke)
      for (unsigned int qp = 0; qp < 4; ++qp)
        local_ke += _JxW[qp] * _coord[qp] * _qp_jacobians[qp];
      _local_ke(_i, _j) += local_ke;
    }

  accumulateTaggedLocalMatrix();
}

void
VectorizedKernel::computeOffDiagJacobian(const unsigned int jvar_num)
{
  const auto & jvar = getVariable(jvar_num);

  if (jvar_num == _var.number())
    computeJacobian();
  else
  {
    prepareMatrixTag(_assembly, _var.number(), jvar_num);
    _nqp = _qrule->n_points();

    // This (undisplaced) jvar could potentially yield the wrong phi size if this object is acting
    // on the displaced mesh
    auto phi_size = jvar.dofIndices().size();
    mooseAssert(
        phi_size * jvar.count() == _local_ke.n(),
        "The size of the phi container does not match the number of local Jacobian columns");

    if (_local_ke.m() != _test.size())
      return;

    precalculateOffDiagJacobian(jvar_num);
    if (jvar.count() == 1)
    {
      for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < phi_size; _j++)
        {
          computeQpOffDiagJacobians(jvar.number());
          #pragma omp simd reduction(+: _local_ke) private(_qp_jacobians)
          for (unsigned int qp = 0; qp < _nqp; qp++)
            _local_ke(_i, _j) += _JxW[qp] * _coord[qp] * _qp_jacobians[qp];
        }
    }
    else
    {
      unsigned int n = phi_size;
      for (_i = 0; _i < _test.size(); _i++)
        for (_j = 0; _j < n; _j++)
          for (unsigned int qp = 0; qp < _nqp; qp++)
          {
            // TODO: Qp-vectorize this
            RealEigenVector v = _JxW[qp] * _coord[qp] *
                                computeQpOffDiagJacobianArray(static_cast<ArrayMooseVariable &>(
                                    const_cast<MooseVariableFieldBase &>(jvar)));
            for (unsigned int k = 0; k < v.size(); ++k)
              _local_ke(_i, _j + k * n) += v(k);
          }
    }

    accumulateTaggedLocalMatrix();
  }
}

void
VectorizedKernel::computeOffDiagJacobianScalar(const unsigned int jvar)
{
  MooseVariableScalar & jv = _sys.getScalarVariable(_tid, jvar);
  prepareMatrixTag(_assembly, _var.number(), jvar);
  _nqp = _qrule->n_points();

  for (_i = 0; _i < _test.size(); _i++)
    for (_j = 0; _j < jv.order(); _j++)
    {
      computeQpOffDiagJacobianScalars(jvar);
  #pragma omp simd reduction(+: _local_ke) private(_qp_jacobians)
      for (_qp = 0; _qp < _nqp; _qp++)
        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * _qp_jacobians[_qp];
    }
  accumulateTaggedLocalMatrix();
}

void
VectorizedKernel::computeResidualAndJacobian()
{
  computeResidual();

  for (const auto & [ivariable, jvariable] : _fe_problem.couplingEntries(_tid, _sys.number()))
  {
    const unsigned int ivar = ivariable->number();
    const unsigned int jvar = jvariable->number();

    if (ivar != _var.number())
      continue;

    if (_is_implicit)
    {
      prepareShapes(jvar);
      computeOffDiagJacobian(jvar);
    }
  }

  /// TODO: add nonlocal Jacobians and scalar Jacobians
}
