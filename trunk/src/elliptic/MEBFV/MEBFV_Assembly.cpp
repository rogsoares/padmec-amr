/*
 * MEBFV_Assembly.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "MEBFV_elliptic.h"

namespace PRS{
	void MEBFV_elliptic::Assembly_A(){
		int i,j;
		int numElements = pGCData->getNumElements();
		double Aij[6], Ajk[6], Aik[6];
		int idxm_IJ[2], idxn_IJ[3];
		int idxm_JK[2], idxn_JK[3];
		int idxm_IK[2], idxn_IK[3];
		int id0, id1, id2;

		for (i=0; i<numElements; i++){
			pGCData->getElementMatrices(i,Aij,Ajk,Aik);
			pMData->getIndices(i,0,idxm_IJ,idxn_IJ);
			pMData->getIndices(i,1,idxm_JK,idxn_JK);
			pMData->getIndices(i,2,idxm_IK,idxn_IK);
			MatSetValues(A,2,idxm_IJ,3,idxn_IJ,Aij,ADD_VALUES);
			MatSetValues(A,2,idxm_JK,3,idxn_JK,Ajk,ADD_VALUES);
			MatSetValues(A,2,idxm_IK,3,idxn_IK,Aik,ADD_VALUES);
		}
		MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

	}

	void MEBFV_elliptic::Assembly_b(){
		MatMult(dirichletMat,dirichletVec,b);	// dirichletVec is a vector containing prescribed values. dirichletVec size: PN
		VecScale(b,-1.0);						// define correction for dirichlet in the rhs
		VecAXPY(b,1.0,SST);						// set source/sink Terms (SST) to RHS vector
		VecAssemblyBegin(b);
		VecAssemblyEnd(b);
	}
}

