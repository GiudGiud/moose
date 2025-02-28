//* This file is part of the MOOSE framework
#include "MeshGenerator.h"
#include "CastUniquePointer.h"

#include "HexagonalAssemblyMeshGenerator.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include <cmath>
#include <vector>

registerMooseObject("MooseApp", HexagonalAssemblyMeshGenerator);

InputParameters
HexagonalAssemblyMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<Real>("side_length", "Side length of a single hexagon");
  params.addRequiredParam<unsigned int>("num_rows", "Number of rows of hexagons");
  params.addRequiredParam<unsigned int>("num_cols", "Number of columns of hexagons");
  params.addParam<subdomain_id_type>("assembly_block_id", 1, "Block ID for the assemblies");
  params.addClassDescription(
      "Generates a 2D quadrilateral mesh composed of hexagonal assemblies (each hexagon subdivided "
      "into 6 quads) arranged in a staggered grid.");
  return params;
}

HexagonalAssemblyMeshGenerator::HexagonalAssemblyMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _side_length(getParam<Real>("side_length")),
    _num_rows(getParam<unsigned int>("num_rows")),
    _num_cols(getParam<unsigned int>("num_cols")),
    _assembly_block_id(getParam<subdomain_id_type>("assembly_block_id"))
{
}

std::unique_ptr<MeshBase>
HexagonalAssemblyMeshGenerator::generate()
{
  auto mesh = buildReplicatedMesh(2);

  // Pre-calculate constants.
  Real s = _side_length;
  Real r = s; // For a regular hexagon, the distance from center to vertex equals the side length.
  Real hex_height = std::sqrt(3.0) * s;
  Real angle_offset = 0.0;

  // Loop over rows and columns (staggered grid).
  for (unsigned int row = 0; row < _num_rows; ++row)
  {
    for (unsigned int col = 0; col < _num_cols; ++col)
    {
      // Compute hexagon center; stagger columns.
      Real cx = col * 1.5 * s;
      Real cy = row * hex_height + ((col % 2) ? (hex_height / 2.0) : 0.0);
      Point center(cx, cy, 0.0);

      // Compute vertices (6 vertices).
      std::vector<Point> vertices;
      for (unsigned int k = 0; k < 6; ++k)
      {
        Real theta = angle_offset + k * M_PI / 3.0;
        vertices.push_back(Point(cx + r * std::cos(theta), cy + r * std::sin(theta), 0.0));
      }

      // Compute midpoints.
      std::vector<Point> midpoints;
      for (unsigned int k = 0; k < 6; ++k)
      {
        Point m = (vertices[k] + vertices[(k + 1) % 6]) * 0.5;
        midpoints.push_back(m);
      }

      // Add nodes to the mesh.
      dof_id_type center_id = mesh->add_point(center)->id();
      std::vector<dof_id_type> vertex_ids(6), midpoint_ids(6);
      for (unsigned int k = 0; k < 6; ++k)
      {
        vertex_ids[k] = mesh->add_point(vertices[k])->id();
        midpoint_ids[k] = mesh->add_point(midpoints[k])->id();
      }

      // Create 6 quadrilaterals covering the hexagon.
      for (unsigned int k = 0; k < 6; ++k)
      {
        std::vector<dof_id_type> conn = {
            vertex_ids[k], midpoint_ids[k], center_id, midpoint_ids[(k + 5) % 6]};
        Elem * elem = Elem::build(QUAD4).release();
        for (unsigned int n = 0; n < 4; ++n)
          elem->set_node(n) = mesh->node_ptr(conn[n]);
        elem->subdomain_id() = _assembly_block_id;
        mesh->add_elem(elem);
      }
    }
  }

  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
