//* This file is part of the MOOSE framework
#include "MeshGenerator.h"
#include "CastUniquePointer.h"
#include "SolidConjugateMeshGenerator.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include <cmath>
#include <vector>

registerMooseObject("MooseApp", SolidConjugateMeshGenerator);

InputParameters
SolidConjugateMeshGenerator::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<unsigned int>("nx", "Number of elements in x direction");
  params.addRequiredParam<unsigned int>("ny", "Number of elements in y direction");
  params.addRequiredParam<unsigned int>("nz", "Number of elements in z direction");
  params.addRequiredParam<Real>("x_length", "Length of domain in x direction");
  params.addRequiredParam<Real>("y_length", "Length of domain in y direction");
  params.addRequiredParam<Real>("z_length", "Length of domain in z direction");
  params.addParam<subdomain_id_type>("block_id", 1, "Block ID for the solid domain");
  params.addClassDescription("Generates a structured 3D hexahedral mesh for the solid domain.");
  return params;
}

SolidConjugateMeshGenerator::SolidConjugateMeshGenerator(const InputParameters & parameters)
  : MeshGenerator(parameters),
    _nx(getParam<unsigned int>("nx")),
    _ny(getParam<unsigned int>("ny")),
    _nz(getParam<unsigned int>("nz")),
    _x_length(getParam<Real>("x_length")),
    _y_length(getParam<Real>("y_length")),
    _z_length(getParam<Real>("z_length")),
    _block_id(getParam<subdomain_id_type>("block_id"))
{
}

std::unique_ptr<MeshBase>
SolidConjugateMeshGenerator::generate()
{
  auto mesh = buildReplicatedMesh(3);
  Real dx = _x_length / _nx;
  Real dy = _y_length / _ny;
  Real dz = _z_length / _nz;

  std::vector<std::vector<std::vector<dof_id_type>>> node_ids(
      _nz + 1, std::vector<std::vector<dof_id_type>>(_ny + 1, std::vector<dof_id_type>(_nx + 1)));

  for (unsigned int k = 0; k <= _nz; ++k)
    for (unsigned int j = 0; j <= _ny; ++j)
      for (unsigned int i = 0; i <= _nx; ++i)
      {
        Real x = i * dx;
        Real y = j * dy;
        Real z = k * dz;
        node_ids[k][j][i] = mesh->add_point(Point(x, y, z))->id();
      }

  // Create HEX8 elements.
  for (unsigned int k = 0; k < _nz; ++k)
    for (unsigned int j = 0; j < _ny; ++j)
      for (unsigned int i = 0; i < _nx; ++i)
      {
        std::vector<dof_id_type> conn = {node_ids[k][j][i],
                                         node_ids[k][j][i + 1],
                                         node_ids[k][j + 1][i + 1],
                                         node_ids[k][j + 1][i],
                                         node_ids[k + 1][j][i],
                                         node_ids[k + 1][j][i + 1],
                                         node_ids[k + 1][j + 1][i + 1],
                                         node_ids[k + 1][j + 1][i]};
        Elem * elem = Elem::build(HEX8).release();
        for (unsigned int n = 0; n < 8; ++n)
          elem->set_node(n) = mesh->node_ptr(conn[n]);
        elem->subdomain_id() = _block_id;
        mesh->add_elem(elem);
      }

  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
