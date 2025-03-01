//* FVUNSMeshReader.C
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//* Licensed under LGPL 2.1, please see LICENSE for details
//*
//* This mesh generator reads a FieldView Unstructured mesh file (.fvuns)
//* and constructs a MeshBase object by reading all the nodes, then creating
//* the elements with proper connectivity, and finally setting the subdomain IDs.
//*
//* Assumptions on the .fvuns file format (binary):
//*   1. A 5-character signature "FVUNS".
//*   2. An integer (4 bytes) indicating the mesh dimension (2 or 3).
//*   3. An integer (4 bytes) for the number of nodes.
//*   4. For each node: three doubles (x,y,z).
//*   5. An integer (4 bytes) for the number of elements.
//*   6. For each element:
//*         - An integer (4 bytes): number of nodes in the element.
//*         - That many integers (4 bytes each): the connectivity (node indices, 0-indexed).
//*         - An integer (4 bytes): the subdomain ID for that element.
//*
//* Depending on your actual .fvuns specification you may need to adjust the reading.
#include "MeshGenerator.h"
#include "CastUniquePointer.h"
#include "FVUNSMeshReader.h"

#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/elem.h"
#include "libmesh/point.h"
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

registerMooseObject("MooseApp", FVUNSMeshReader);

InputParameters
FVUNSMeshReader::validParams()
{
  InputParameters params = MeshGenerator::validParams();
  params.addRequiredParam<std::string>("filename",
                                       "The name of the FVUNS mesh file (.fvuns) to read");
  params.addClassDescription(
      "Reads a FieldView Unstructured mesh file (.fvuns) into a MOOSE MeshBase object.");
  return params;
}

FVUNSMeshReader::FVUNSMeshReader(const InputParameters & parameters)
  : MeshGenerator(parameters), _filename(getParam<std::string>("filename"))
{
}

std::unique_ptr<MeshBase>
FVUNSMeshReader::generate()
{
  // Open the .fvuns file in binary mode.
  std::ifstream in(_filename.c_str(), std::ios::binary);
  if (!in)
    mooseError("Failed to open FVUNS file: ", _filename);

  // Read and verify the file signature.
  char signature[6];
  in.read(signature, 5);
  signature[5] = '\0';
  if (std::string(signature) != "FIELDVIEW")
    mooseError("Invalid FVUNS file signature in file: ", _filename);

  // Read the mesh dimension (assumed to be stored as an int).
  int dim = 0;
  in.read(reinterpret_cast<char *>(&dim), sizeof(int));
  if (dim < 2 || dim > 3)
    mooseError("Unsupported mesh dimension in FVUNS file: ", dim);

  // Create a replicated mesh based on the dimension.
  std::unique_ptr<MeshBase> mesh;
  if (dim == 2)
    mesh = buildReplicatedMesh(2);
  else
    mesh = buildReplicatedMesh(3);

  // Read number of nodes.
  int num_nodes = 0;
  in.read(reinterpret_cast<char *>(&num_nodes), sizeof(int));

  // Read each node and add to the mesh.
  for (int i = 0; i < num_nodes; i++)
  {
    double x, y, z;
    in.read(reinterpret_cast<char *>(&x), sizeof(double));
    in.read(reinterpret_cast<char *>(&y), sizeof(double));
    in.read(reinterpret_cast<char *>(&z), sizeof(double));
    Point p(x, y, z);
    mesh->add_point(p);
  }

  // Read number of elements.
  int num_elems = 0;
  in.read(reinterpret_cast<char *>(&num_elems), sizeof(int));

  // For each element, read its connectivity and subdomain id.
  for (int i = 0; i < num_elems; i++)
  {
    int elem_node_count = 0;
    in.read(reinterpret_cast<char *>(&elem_node_count), sizeof(int));
    std::vector<int> connectivity(elem_node_count);
    for (int j = 0; j < elem_node_count; j++)
    {
      in.read(reinterpret_cast<char *>(&connectivity[j]), sizeof(int));
    }
    int subdomain_id = 0;
    in.read(reinterpret_cast<char *>(&subdomain_id), sizeof(int));

    // Create an element based on the number of nodes.
    Elem * elem = nullptr;
    if (elem_node_count == 3)
      elem = Elem::build(TRI3).release();
    else if (elem_node_count == 4)
      elem = Elem::build(QUAD4).release();
    else if (elem_node_count == 8)
      elem = Elem::build(HEX8).release();
    else
      mooseError("Unsupported element type with ", elem_node_count, " nodes in FVUNS file.");

    // Set the nodes of the element.
    for (unsigned int j = 0; j < connectivity.size(); j++)
      elem->set_node(j) = mesh->node_ptr(connectivity[j]);
    elem->subdomain_id() = subdomain_id;
    mesh->add_elem(elem);
  }

  // Finalize mesh.
  mesh->prepare_for_use();
  return dynamic_pointer_cast<MeshBase>(mesh);
}
