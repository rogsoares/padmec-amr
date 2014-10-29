/*
 * MeshData2.cpp
 *
 *  Created on: Oct 27, 2014
 *      Author: rogerio
 */

#include "MeshData.h"

namespace PRS{

	void MeshData::mapNodeID(pMesh theMesh){
		pEntity node;
		int i,ID;

		i = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			ID = EN_id(node);
			nodeID_map[ID] = i++;
		}
		VIter_delete(vit);
	}

	void MeshData::getMappedNodeID(int ID, int& mappedID){
		mappedID = nodeID_map[ID];
	}

	void MeshData::getMappedNodeID(int ID0, int ID1, int ID2, int& mappedID0, int& mappedID1, int& mappedID2){
		mappedID0 = nodeID_map[ID0];
		mappedID1 = nodeID_map[ID1];
		mappedID2 = nodeID_map[ID2];
	}

	void MeshData::getRowsAndCols(int id0, int id1, int* idxm_IJ, int* idxn_IJ){

	}

	IS MeshData::get_DirichletISRows() const{

	}

	IS MeshData::get_DirichletISCols() const{

	}

	int MeshData::nrows() const{

	}

	int MeshData::ncols() const{

	}

	int MeshData::nfreen() const{

	}
}
