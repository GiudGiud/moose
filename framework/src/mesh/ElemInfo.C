//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElemInfo.h"

ElemInfo::ElemInfo(const Elem * const elem)
  : _elem(elem),
    _volume(_elem->volume()),
    _centroid(_elem->vertex_average()),
    _coord_transform_factor(1.0),
    _dof_indices(std::vector<std::vector<dof_id_type>>())
{
}

std::string
ElemInfo::print(bool print_full_elem_info) const
{
  std::string out = "";
  if (print_full_elem_info)
    _elem->print_info();
  else
  {
    out += "id:       " + std::to_string(_elem->id()) + "\n";
    out += "centroid: " + Moose::stringify(_centroid) + "\n";
    out += "dofs:     " + Moose::stringify(_dof_indices) + "\n";
  }
  return out;
}
