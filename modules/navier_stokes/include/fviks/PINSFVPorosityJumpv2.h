//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVInterfaceKernel.h"

class PINSFVPorosityJumpv2 : public FVInterfaceKernel
{
public:
  static InputParameters validParams();
  PINSFVPorosityJumpv2(const InputParameters & params);

protected:
  ADReal computeQpResidual() override {return 0;};
  void computeResidual(const FaceInfo & fi) override;
  void computeJacobian(const FaceInfo & fi) override;

private:

  const unsigned int _dim;

  /// Pressure on one side
  const ADVariableValue & _p_elem;

  /// Pressure on the other side
  const ADVariableValue & _p_neighbor;

  /// Porosity functor
  const Moose::Functor<ADReal> & _eps;

  /// Density functor
  const Moose::Functor<ADReal> & _rho;

  /// Velocity X-component functor
  const Moose::Functor<ADReal> & _u;

  /// Velocity Y-component functor
  const Moose::Functor<ADReal> * const _v;

  /// Velocity Z-component functor
  const Moose::Functor<ADReal> * const _w;

  /// which momentum component this kernel applies to
  const int _index;
};
