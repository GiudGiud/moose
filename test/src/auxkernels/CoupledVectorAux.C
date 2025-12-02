//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CoupledVectorAux.h"

registerMooseObject("MooseTestApp", CoupledVectorAux);

InputParameters
CoupledVectorAux::validParams()
{
  InputParameters params = VectorAuxKernel::validParams();

  MooseEnum operators("value gradient laplacian dot");

  params.addRequiredCoupledVar("coupled", "Coupled value to apply the operator on");
  params.addRequiredParam<MooseEnum>(
      "operator", operators, "The operator to use in the calculation");

  params.addParam<unsigned int>(
      "component_1",
      0,
      "Row index to return for gradient operation, and first slice index for laplacian operation");
  params.addParam<unsigned int>(
      "component_2", 0, "Second slice index to return for laplacian operation");

  return params;
}

CoupledVectorAux::CoupledVectorAux(const InputParameters & parameters)
  : VectorAuxKernel(parameters),
    _operator(getParam<MooseEnum>("operator")),
    _coupled_val(coupledVectorValue("coupled")),
    _coupled_grad(coupledVectorGradient("coupled")),
    _coupled_lapl(coupledVectorSecond("coupled")),
    _coupled_dot(coupledVectorDot("coupled")),
    _comp_1(getParam<unsigned int>("component_1")),
    _comp_2(getParam<unsigned int>("component_2"))
{
}

RealVectorValue
CoupledVectorAux::computeValue()
{
  switch (_operator)
  {
    case 0:
      return _coupled_val[_qp];
    case 1:
      return RealVectorValue(_coupled_grad[_qp](_comp_1, 0),
                             _coupled_grad[_qp](_comp_1, 1),
                             _coupled_grad[_qp](_comp_1, 2));
    case 2:
      return _coupled_lapl[_qp].slice(_comp_1).slice(_comp_2);
    case 3:
      return _coupled_dot[_qp];
    default:
      mooseError("Should not reach here");
  }
}
