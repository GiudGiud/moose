//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Diffusion.h"

registerMooseObject("MooseApp", Diffusion);

InputParameters
Diffusion::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

Diffusion::Diffusion(const InputParameters & parameters) : Kernel(parameters) {}

void
Diffusion::computeResidual()
{
  prepareVectorTag(_assembly, _var.number());
  const unsigned int n_test = _test.size();
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
#pragma omp simd
    for (_i = 0; _i < n_test; _i++) // target for auto vectorization
      _local_re(_i) += _grad_u[_qp] * _JxW[_qp] * _coord[_qp] * _grad_test[_i][_qp];
  accumulateTaggedLocalResidual();
}

void
Diffusion::computeJacobian()
{
  prepareMatrixTag(_assembly, _var.number(), _var.number());

  const unsigned int n_test = _test.size();
  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    for (_j = 0; _j < _phi.size(); _j++)
#pragma omp simd
      for (_i = 0; _i < n_test; _i++) // target for auto vectorization
        _local_ke(_i, _j) += _grad_phi[_j][_qp] * _JxW[_qp] * _coord[_qp] * _grad_test[_i][_qp];

  accumulateTaggedLocalMatrix();
}

Real
Diffusion::computeQpResidual()
{
  mooseError("AJJAJA residual");
  return _grad_u[_qp] * _grad_test[_i][_qp];
}

Real
Diffusion::computeQpJacobian()
{
  mooseError("HAHAHAJA");
  return _grad_phi[_j][_qp] * _grad_test[_i][_qp];
}
