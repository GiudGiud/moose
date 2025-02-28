//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "2DExpansionConeMeshGenerator.h"
#include "CastUniquePointer.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "MooseMesh.h"

#include <cmath>
#include <vector>

/// Register the new mesh generator with MOOSE
registerMooseObject("MooseApp", QuadrilateralExpansionConeMeshGenerator);

//
// The QuadrilateralExpansionConeMeshGenerator creates a 2D mesh for an annular sector (or cone)
// by generating nodes on (num_radial_layers+1) radial levels (with optionally geometric spacing)
// and (num_angular_divisions+1) equally spaced angular positions between theta_min and theta_max.
// It then creates QUAD4 elements between adjacent nodes.
//
InputParameters
QuadrilateralExpansionConeMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();

  params.addParam<Point>("apex", Point(0.0, 0.0, 0.0), "The apex (origin) of the expansion cone");
  params.addRequiredParam<Real>("inner_radius", "The inner radius of the cone (must be > 0)");
  params.addRequiredParam<Real>("outer_radius",
                                "The outer radius of the cone (must be > inner_radius)");
  params.addRequiredParam<Real>("theta_min", "The minimum angle (in radians) of the cone sector");
  params.addRequiredParam<Real>("theta_max", "The maximum angle (in radians) of the cone sector");
  params.addRequiredParam<unsigned int>("num_radial_layers",
                                        "Number of layers in the radial direction");
  params.addRequiredParam<unsigned int>("num_angular_divisions",
                                        "Number of divisions in the angular direction");
  params.addParam<Real>(
      "radial_expansion",
      1.0,
      "Radial expansion factor for non-uniform spacing (1.0 for uniform spacing)");
  params.addParam<subdomain_id_type>(
      "block_id", 0, "Subdomain (block) id for the generated elements");

  params.addClassDescription(
      "Generates a 2D quadrilateral mesh for an expansion cone (annular sector). "
      "Nodes are placed on radial and angular grids with optional geometric spacing in the radial "
      "direction.");
  return params;
}

QuadrilateralExpansionConeMeshGenerator::QuadrilateralExpansionConeMeshGenerator(
    const InputParameters & parameters)
  : MeshGenerator(parameters),
    _apex(getParam<Point>("apex")),
    _inner_radius(getParam<Real>("inner_radius")),
    _outer_radius(getParam<Real>("outer_radius")),
    _theta_min(getParam<Real>("theta_min")),
    _theta_max(getParam<Real>("theta_max")),
    _num_radial_layers(getParam<unsigned int>("num_radial_layers")),
    _num_angular_divisions(getParam<unsigned int>("num_angular_divisions")),
    _radial_expansion(getParam<Real>("radial_expansion")),
    _block_id(getParam<subdomain_id_type>("block_id"))
{
  if (_inner_radius <= 0)
    mooseError("inner_radius must be > 0");
  if (_outer_radius <= _inner_radius)
    mooseError("outer_radius must be greater than inner_radius");
  if (_theta_max <= _theta_min)
    mooseError("theta_max must be greater than theta_min");
  if (_radial_expansion <= 0)
    mooseError("radial_expansion must be > 0");
}

std::unique_ptr<MeshBase>
QuadrilateralExpansionConeMeshGenerator::generate()
{
  // Build a 2D replicated mesh
  auto mesh = buildMeshBaseObject();

  // Compute radial coordinates (layer radii)
  // If radial_expansion == 1, use uniform spacing; otherwise, use a geometric progression.
  std::vector<Real> radial_coords(_num_radial_layers + 1);
  if (std::fabs(_radial_expansion - 1.0) < 1e-8)
  {
    Real dr = (_outer_radius - _inner_radius) / _num_radial_layers;
    for (unsigned int i = 0; i <= _num_radial_layers; ++i)
      radial_coords[i] = _inner_radius + i * dr;
  }
  else
  {
    // Geometric progression: find the initial step such that
    // sum_{i=0}^{N-1} step * (radial_expansion^i) = outer_radius - inner_radius.
    Real denom =
        (std::pow(_radial_expansion, _num_radial_layers) - 1.0) / (_radial_expansion - 1.0);
    Real step = (_outer_radius - _inner_radius) / denom;
    radial_coords[0] = _inner_radius;
    for (unsigned int i = 1; i <= _num_radial_layers; ++i)
      radial_coords[i] = radial_coords[i - 1] + step * std::pow(_radial_expansion, i - 1);
  }

  // Compute angular coordinates (in radians)
  std::vector<Real> angular_coords(_num_angular_divisions + 1);
  Real dtheta = (_theta_max - _theta_min) / _num_angular_divisions;
  for (unsigned int j = 0; j <= _num_angular_divisions; ++j)
    angular_coords[j] = _theta_min + j * dtheta;

  // Create a grid of nodes. We store node ids in a 2D array:
  std::vector<std::vector<dof_id_type>> node_ids(
      _num_radial_layers + 1, std::vector<dof_id_type>(_num_angular_divisions + 1));

  // Loop over each radial level and angular division to add nodes.
  for (unsigned int i = 0; i <= _num_radial_layers; ++i)
  {
    for (unsigned int j = 0; j <= _num_angular_divisions; ++j)
    {
      Real r = radial_coords[i];
      Real theta = angular_coords[j];
      // Map polar coordinates to Cartesian: p = apex + (r*cos(theta), r*sin(theta))
      Point p(_apex(0) + r * std::cos(theta),
              _apex(1) + r * std::sin(theta),
              _apex(2)); // For 2D, z remains unchanged
      node_ids[i][j] = mesh->add_point(p)->id();
    }
  }

  // Create quadrilateral elements.
  // For each cell in the grid (i = 0..num_radial_layers-1, j = 0..num_angular_divisions-1),
  // define a QUAD4 element with nodes ordered counter-clockwise.
  for (unsigned int i = 0; i < _num_radial_layers; ++i)
  {
    for (unsigned int j = 0; j < _num_angular_divisions; ++j)
    {
      // Define the connectivity: lower left, lower right, upper right, upper left.
      std::vector<dof_id_type> connectivity = {
          node_ids[i][j], node_ids[i][j + 1], node_ids[i + 1][j + 1], node_ids[i + 1][j]};

      // Create a new QUAD4 element.
      Elem * elem = Elem::build(QUAD4).release();
      for (unsigned int k = 0; k < connectivity.size(); ++k)
        elem->set_node(k) = mesh->node_ptr(connectivity[k]);
      elem->subdomain_id() = _block_id;
      mesh->add_elem(elem);
    }
  }

  // Finalize and prepare the mesh
  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
