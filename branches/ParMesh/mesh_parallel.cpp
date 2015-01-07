/*
 * mesh_parallel.cpp
 *
 *  Created on: Jan 7, 2015
 *      Author: rogerio
 */

#include "mesh.h"

namespace MeshDB{

	int Mesh::get_rank() const{
		return rank;
	}

	int Mesh::get_nproc() const{
		return nproc;
	}


	// The routine ParMETIS V3 PartMeshKway takes a distributed mesh and computes its partitioning, while
	// ParMETIS V3 Mesh2Dual takes a distributed mesh and constructs a distributed dual graph. Both of these rou-
	// tines require an elmdist array that specifies the distribution of the mesh elements, but that is otherwise identical
	// to the vtxdist array. They also require a pair of arrays called eptr and eind, as well as the integer parameter
	// ncommonnodes.
	// The eptr and eind arrays are similar in nature to the xadj and adjncy arrays used to specify the adjacency
	// list of a graph but now for each element they specify the set of nodes that make up each element. Specifically, the set
	// of nodes that belong to element i is stored in array eind starting at index eptr[i] and ending at (but not including)
	// index eptr[i + 1] (in other words, eind[eptr[i]] up through and including eind[eptr[i + 1]-1]). Hence,
	// the node lists for each element are stored consecutively in the array eind. This format allows the specification of
	// meshes that contain elements of mixed type.
	// The ncommonnodes parameter specifies the degree of connectivity that is desired between the vertices of the
	// dual graph. Specifically, an edge is placed between two vertices if their corresponding mesh elements share at least
	// g nodes, where g is the ncommonnodes parameter. Hence, this parameter can be set to result in a traditional dual
	// graph (e.g., a value of two for a triangle mesh or a value of four for a hexahedral mesh). However, it can also be set
	// higher or lower for increased or decreased connectivity.
	// Additionally, ParMETIS V3 PartMeshKway requires an elmwgt array that is analogous to the vwgt array.
	void Mesh::mesh_partition(){

	}

	void Mesh::mesh_distribution(){

	}

	void Mesh::bdry_linkSetup(){

	}
}

