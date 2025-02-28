//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//* Licensed under LGPL 2.1, please see LICENSE for details

#include "MeshGenerator.h"
#include "CastUniquePointer.h"

#include "GlobalSketchMeshGenerator.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include <cmath>
#include <vector>

registerMooseObject("MooseApp", GlobalSketchMeshGenerator);

InputParameters
GlobalSketchMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<Real>("width", "Total width of the EDT domain");
  params.addRequiredParam<Real>("height", "Total height of the EDT domain");
  params.addRequiredParam<Real>("chimney_radius", "Radius of the central air chimney");
  params.addParam<Point>("chimney_center",
                         Point(0.5, 0.5, 0.0),
                         "Center of the air chimney in normalized coordinates");
  params.addParam<subdomain_id_type>("solid_block_id", 1, "Block ID for the solid region");
  params.addParam<subdomain_id_type>("chimney_block_id", 2, "Block ID for the air chimney region");
  params.addRequiredParam<unsigned int>("nx", "Number of elements in x direction");
  params.addRequiredParam<unsigned int>("ny", "Number of elements in y direction");
  params.addClassDescription(
      "Generates a 2D global sketch of the EDT with a central circular air chimney.");
  return params;
}

GlobalSketchMeshGenerator::GlobalSketchMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _width(getParam<Real>("width")),
    _height(getParam<Real>("height")),
    _chimney_radius(getParam<Real>("chimney_radius")),
    _chimney_center(getParam<Point>("chimney_center")),
    _solid_block_id(getParam<subdomain_id_type>("solid_block_id")),
    _chimney_block_id(getParam<subdomain_id_type>("chimney_block_id")),
    _nx(getParam<unsigned int>("nx")),
    _ny(getParam<unsigned int>("ny"))
{
}

std::unique_ptr<MeshBase>
GlobalSketchMeshGenerator::generate()
{
  auto mesh = buildReplicatedMesh(2);
  Real dx = _width / _nx;
  Real dy = _height / _ny;

  // Build a structured grid of nodes.
  std::vector<std::vector<dof_id_type>> node_ids(_ny + 1, std::vector<dof_id_type>(_nx + 1));
  for (unsigned int j = 0; j <= _ny; ++j)
    for (unsigned int i = 0; i <= _nx; ++i)
    {
      Real x = i * dx;
      Real y = j * dy;
      node_ids[j][i] = mesh->add_point(Point(x, y, 0.0))->id();
    }

  // Chimney center in physical coordinates.
  Real cx = _chimney_center(0) * _width;
  Real cy = _chimney_center(1) * _height;

  // Create elements and assign block id based on distance from chimney center.
  for (unsigned int j = 0; j < _ny; ++j)
    for (unsigned int i = 0; i < _nx; ++i)
    {
      std::vector<dof_id_type> conn = {
          node_ids[j][i], node_ids[j][i + 1], node_ids[j + 1][i + 1], node_ids[j + 1][i]};
      Elem * elem = Elem::build(QUAD4).release();
      for (unsigned int k = 0; k < 4; ++k)
        elem->set_node(k) = mesh->node_ptr(conn[k]);

      // Compute centroid.
      Point centroid(0.0, 0.0, 0.0);
      for (unsigned int k = 0; k < 4; ++k)
        centroid += *mesh->node_ptr(conn[k]);
      centroid /= 4.0;

      // Use the distance from (cx,cy) to decide the block.
      Real dist = std::sqrt(std::pow(centroid(0) - cx, 2) + std::pow(centroid(1) - cy, 2));
      elem->subdomain_id() = (dist < _chimney_radius) ? _chimney_block_id : _solid_block_id;

      mesh->add_elem(elem);
    }

  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
