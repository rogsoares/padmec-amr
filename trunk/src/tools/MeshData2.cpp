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
		msize = M_numVertices(theMesh);
	}

	void MeshData::findDirichletNodes(){
		std::set<int> dirichletnodes_set;
		std::set<int> freenodes_set;
		std::set<int>::iterator iter_set;
		int i, nidxm, nidxn, ID, idx;
		// find dirichlet nodes
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			flag = GEN_tag( node->getClassification() );
			iter_set = dirichletflag_set.find(flag);
			if (iter_set!=dirichlet_set.end()){
				dirichletnodes_set.insert(EN_id(node));
			}
			else{
				freenodes_set.insert(EN_id(node));
			}
		}
		VIter_delete(vit);

		// defines variables and allocate memory for vectors
		nnfreen = (int)dirichletnodes_set.size();
		nfreen = msize - nnfreen;
		dirichlet_idx = new int[nnfreen];
		dirichlet_data = new double[nnfreen];
		i = 0;
		for(iter_set=dirichletnodes_set.begin(); iter_set!=dirichletnodes_set.end(); iter_set++){
			node = *iter_set
			ID = EN_id(node);
			flag = GEN_tag( node->getClassification() );
			getMappedNodeID(ID,idx);
			dirichlet_idx[i] = idx;
			dirichlet_data[i] = pSimPar->getBC_Value(flag);
			i++;
		}
		dirichletnodes_set.clear();
	}

	void MeshData::getMappedIndices_free(int i_th_element, int i_th_edge, int* idxm_IJ, int* idxn_IJ){
		int row = i_th_element;
		int col = i_th_edge;
		idxm_IJ[0] = rows_cols_indices_free_mat.getValue(row,5*col);
		idxm_IJ[1] = rows_cols_indices_free_mat.getValue(row,5*col+1);
		idxn_IJ[0] = rows_cols_indices_free_mat.getValue(row,5*col+2);
		idxn_IJ[1] = rows_cols_indices_free_mat.getValue(row,5*col+3);
		idxn_IJ[2] = rows_cols_indices_free_mat.getValue(row,5*col+4);
	}

	void MeshData::getMappedIndices_dirichlet(int i_th_element, int* idx_IJ, int* idx_JK, int idx_IK){
		int row = i_th_element;
		int col = i_th_edge;
		idxm_IJ[0] = rows_cols_indices_nfree_mat.getValue(row,5*col);
		idxm_IJ[1] = rows_cols_indices_nfree_mat.getValue(row,5*col+1);
		idxn_IJ[0] = rows_cols_indices_nfree_mat.getValue(row,5*col+2);
		idxn_IJ[1] = rows_cols_indices_nfree_mat.getValue(row,5*col+3);
		idxn_IJ[2] = rows_cols_indices_nfree_mat.getValue(row,5*col+4);
	}

	int MeshData::nrows() const{
		return msize;
	}

	int MeshData::ncols() const{
		return msize;
	}

	int MeshData::numFreeNodes() const{
		return nfreen;
	}

	const int* MeshData::getDirichlet_idx() const{
		return dirichlet_idx;
	}

	const double* MeshData::getDirichlet_data() const{
		return dirichlet_data;
	}

}
