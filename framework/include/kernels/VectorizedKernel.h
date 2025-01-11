//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "KernelBase.h"
#include "MooseVariableInterface.h"

/**
 * A base class for kernels that use SIMD vectorization over the quadrature point vector
 */
class VectorizedKernel : public KernelBase, public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  VectorizedKernel(const InputParameters & parameters);

  /// Compute this VectorizedKernel's contribution to the residual
  virtual void computeResidual() override;

  /// Compute this VectorizedKernel's contribution to the diagonal Jacobian entries
  virtual void computeJacobian() override;

  /// Computes d-residual / d-jvar... storing the result in Ke.
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

  /// Compute the residual and Jacobian together
  virtual void computeResidualAndJacobian() override;

  /**
   * Computes jacobian block with respect to a scalar variable
   * @param jvar The number of the scalar variable
   */
  virtual void computeOffDiagJacobianScalar(unsigned int jvar) override;

  virtual const MooseVariable & variable() const override { return _var; }

protected:
  /**
   * Compute this VectorizedKernel's contribution to the residual at the current quadrature point
   */
  virtual void computeQpResiduals() = 0;

  /**
   * Compute this VectorizedKernel's contribution to the Jacobian at the current quadrature point
   */
  virtual void computeQpJacobians() = 0;

  /**
   * For coupling standard variables
   */
  virtual void computeQpOffDiagJacobians(unsigned int /*jvar*/) { }

  /**
   * For coupling scalar variables
   */
  virtual void computeQpOffDiagJacobianScalars(unsigned int /*jvar*/) { }

  /**
   * For coupling array variables
   */
  virtual RealEigenVector computeQpOffDiagJacobianArray(const ArrayMooseVariable & jvar)
  {
    return RealEigenVector::Zero(jvar.count());
  }

  /// This is a regular VectorizedKernel so we cast to a regular MooseVariable
  MooseVariable & _var;

  /// the current test function
  const VariableTestValue & _test;

  /// gradient of the test function
  const VariableTestGradient & _grad_test;

  /// the current shape functions
  const VariablePhiValue & _phi;

  /// gradient of the shape function
  const VariablePhiGradient & _grad_phi;

  /// Holds the solution at current quadrature points
  const VariableValue & _u;

  /// Holds the solution gradient at the current quadrature points
  const VariableGradient & _grad_u;

  /// Holds the number of quadrature points in the current element
  unsigned int _nqp;

  /// Holds the residual at the current quadrature points
  VariableValue _qp_residuals;

  /// Holds a Jacobian coefficient (diagonal or not) at the current quadrature points
  VariableValue _qp_jacobians;
};
