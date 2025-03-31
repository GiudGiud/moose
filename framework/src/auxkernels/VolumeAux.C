//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VolumeAux.h"

registerMooseObject("MooseApp", VolumeAux);

InputParameters
VolumeAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription(
      "Auxiliary Kernel that samples element, side or nodal volumes. Nodal volumes are defined as "
      "the sum of the neighbor element volumes divided by their number of nodes");
  return params;
}

VolumeAux::VolumeAux(const InputParameters & parameters) : AuxKernel(parameters)
{
  if (mooseVariableBase()->feType() != libMesh::FEType(CONSTANT, MONOMIAL) &&
      mooseVariableBase()->feType() != libMesh::FEType(FIRST, LAGRANGE))
    paramError("variable", "Must be of type CONSTANT MONOMIAL");

  // Nodal volumes are fairly simple right now
  if (_mesh.hasSecondOrderElements())
    paramError("variable", "Nodal volumes cannot currently be computed for a second order mesh");
  if (_mesh.getUniqueCoordSystem() != Moose::COORD_XYZ)
    paramError("variable",
               "Nodal variables cannot currently be computed for non-Cartesian coordinates");
}

Real
VolumeAux::computeValue()
{
  if (isNodal())
  {
    if (_bnd)
      paramError("variable",
                 "Nodal side volumes are not implemented. Consider using NodalArea from the "
                 "Contact module");
    // Loop over the node neighbor elements and add their contribution to the nodal volume
    const auto & node_to_elems = _mesh.nodeToElemMap();
    Real volume = 0;
    for (const auto elem_id : node_to_elems.at(_current_node->id()))
    {
      const auto elem = _mesh.elemPtr(elem_id);
      // Alternatively, we would need to use the quadrature and the coordinate system weights
      if (elem->active())
        volume += elem->volume() / elem->n_nodes();
    }
    return volume;
  }
  else
    return _bnd ? _current_side_volume : _current_elem_volume;
}

void
VolumeAux::compute()
{
  _var.setDofValue(computeValue(), 0);
}
