/*
 * MEBFV_Assembly.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "MEBFV_elliptic.h"

namespace PRS{
	void MEBFV_elliptic::Assembly_A(){
		int numElements = pGCData->getNumElements();
		double Aij[6], Ajk[6], Aik[6];
		int idxm_IJ[2], idxn_IJ[3];
		int idxm_JK[2], idxn_JK[3];
		int idxm_IK[2], idxn_IK[3];
		int id0, id1, id2;

		for (int i=0; i<numElements; i++){
			pGCData->getIJK_IDs(i,id0,id1,id2);
			pGCData->getElementMatrices(i,Aij,Ajk,Aik);
			pMData->getRowsAndCols(id0,id1,idxm_IJ,idxn_IJ);
			pMData->getRowsAndCols(id1,id2,idxm_JK,idxn_JK);
			pMData->getRowsAndCols(id0,id2,idxm_IK,idxn_IK);
			MatSetValues(A,2,idxm_IJ,3,idxn_IJ,Aij,ADD_VALUES);
			MatSetValues(A,2,idxm_JK,3,idxn_JK,Ajk,ADD_VALUES);
			MatSetValues(A,2,idxm_IK,3,idxn_IK,Aik,ADD_VALUES);
		}
		MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	}

	void MEBFV_elliptic::Assembly_b(){
		// get matrix contribution from prescribed nodes (dirichletMat)
		// pMData->get_isrow(): all matrix rows of free nodes.         (size=FN, FN=Free Nodes)
		// pMData->get_iscol(): all matrix columns of dirichlet nodes. (size=PN, PN=Prescribed Nodes)
		// dirichletMat size: FN x PN
		MatGetSubMatrix(A,pMData->get_DirichletISRows(),pMData->get_DirichletISCols(),PETSC_DECIDE,MAT_INITIAL_MATRIX,&dirichletMat);
		MatMult(dirichletMat,dirichletVec,b);	// dirichletVec is a vector containing prescribed values. dirichletVec size: PN
		VecScale(b,-1.0);						// define correction for dirichlet in the rhs
		VecAXPY(b,1.0,SST);						// set source/sink Terms (SST) to RHS vector
		VecAssemblyBegin(b);
		VecAssemblyEnd(b);
	}
}

