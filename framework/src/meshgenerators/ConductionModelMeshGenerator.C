//* This file is part of the MOOSE framework
#include "MeshGenerator.h"
#include "CastUniquePointer.h"

#include "ConductionModelMeshGenerator.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include <cmath>
#include <vector>

registerMooseObject("MooseApp", ConductionModelMeshGenerator);

InputParameters
ConductionModelMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<unsigned int>("nx", "Number of elements in x direction");
  params.addRequiredParam<unsigned int>("ny", "Number of elements in y direction");
  params.addRequiredParam<Real>("width", "Width of the domain");
  params.addRequiredParam<Real>("height", "Height of the domain");
  params.addParam<Real>(
      "curvature_amplitude", 0.0, "Amplitude of sinusoidal curvature on the top boundary");
  params.addParam<subdomain_id_type>("block_id", 1, "Block ID for the domain");
  params.addClassDescription("Generates a 2D quadrilateral mesh for a conduction model, with "
                             "optional top boundary curvature.");
  return params;
}

ConductionModelMeshGenerator::ConductionModelMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _nx(getParam<unsigned int>("nx")),
    _ny(getParam<unsigned int>("ny")),
    _width(getParam<Real>("width")),
    _height(getParam<Real>("height")),
    _curvature_amplitude(getParam<Real>("curvature_amplitude")),
    _block_id(getParam<subdomain_id_type>("block_id"))
{
}

std::unique_ptr<MeshBase>
ConductionModelMeshGenerator::generate()
{
  auto mesh = buildReplicatedMesh(2);
  Real dx = _width / _nx;
  Real dy = _height / _ny;

  // Create nodes (apply curvature to top row if specified)
  std::vector<std::vector<dof_id_type>> node_ids(_ny + 1, std::vector<dof_id_type>(_nx + 1));
  for (unsigned int j = 0; j <= _ny; ++j)
  {
    for (unsigned int i = 0; i <= _nx; ++i)
    {
      Real x = i * dx;
      Real y = j * dy;
      if (j == _ny && _curvature_amplitude != 0.0)
        y += _curvature_amplitude * std::sin(M_PI * x / _width);
      node_ids[j][i] = mesh->add_point(Point(x, y, 0.0))->id();
    }
  }

  // Create QUAD4 elements.
  for (unsigned int j = 0; j < _ny; ++j)
  {
    for (unsigned int i = 0; i < _nx; ++i)
    {
      std::vector<dof_id_type> conn = {
          node_ids[j][i], node_ids[j][i + 1], node_ids[j + 1][i + 1], node_ids[j + 1][i]};
      Elem * elem = Elem::build(QUAD4).release();
      for (unsigned int k = 0; k < 4; ++k)
        elem->set_node(k) = mesh->node_ptr(conn[k]);
      elem->subdomain_id() = _block_id;
      mesh->add_elem(elem);
    }
  }

  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
