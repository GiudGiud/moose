//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorizedDiffusion.h"

registerMooseObject("MooseApp", VectorizedDiffusion);

InputParameters
VectorizedDiffusion::validParams()
{
  InputParameters params = VectorizedKernel::validParams();
  params.addClassDescription("The Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
                             "form of $(\\nabla \\phi_i, \\nabla u_h)$, using a SIMD-vectorized implementation.");
  return params;
}

VectorizedDiffusion::VectorizedDiffusion(const InputParameters & parameters) : VectorizedKernel(parameters) {}

void
VectorizedDiffusion::computeQpResiduals()
{
#pragma omp simd //private(_qp_residuals)
  for (unsigned int qp = 0; qp < _nqp; qp++)
     _qp_residuals[qp] = 0.;
  for (unsigned int qp = 0; qp < _nqp; qp++)
    // It tends to unroll this loop instead of vectorizing it
    // #pragma omp simd private(_qp_residuals)
      for (unsigned int i = 0; i < LIBMESH_DIM; i++)
        _qp_residuals[qp] += _grad_u[qp]._coords[i] * _grad_test[_i][qp]._coords[i];
  // std::cout << 0 << " " << _qp_residuals[0] << " " << _grad_u[0] << " " << _grad_test[_i][0] << std::endl;
  // std::cout << 1 << " " << _qp_residuals[1] << " " << _grad_u[1] << " " << _grad_test[_i][1] << std::endl;
  // std::cout << 2 << " " << _qp_residuals[2] << " " << _grad_u[2] << " " << _grad_test[_i][2] << std::endl;
  // std::cout << 3 << " " << _qp_residuals[3] << " " << _grad_u[3] << " " << _grad_test[_i][3] << std::endl;
}

void
VectorizedDiffusion::computeQpJacobians()
{
#pragma omp simd //private(_qp_jacobians)
  for (unsigned int qp = 0; qp < _nqp; qp++)
    _qp_jacobians[qp] = 0.;

  for (unsigned int qp = 0; qp < _nqp; qp++)
  // It tends to unroll this loop instead of vectorizing it
    #pragma omp simd //private(_qp_jacobians)
    for (unsigned int i = 0; i < LIBMESH_DIM; i++)
      _qp_jacobians[qp] += _grad_phi[_j][qp]._coords[i] * _grad_test[_i][qp]._coords[i];
}
