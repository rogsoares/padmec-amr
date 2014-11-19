/*
 * MEBFV_Initialize.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "MEBFV_elliptic.h"

namespace PRS{
	void MEBFV_elliptic::Initialize(){

		if (initialize){
			return;
		}

		int nrows = pMData->nrows();	// number of mesh nodes
		int ncols = pMData->ncols();	// number of mesh nodes
		int nfnodes = pMData->numFreeNodes(); // number of mesh free nodes
		int nnfnodes = nrows - nfnodes;	// number of dirichlet nodes (prescribed)

		// define A matrix
		MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nrows,ncols,80,PETSC_NULL,80,PETSC_NULL,&A);
		MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,nrows,ncols,80,PETSC_NULL,80,PETSC_NULL,&dirichletMat);

		// define b vector
		VecCreate(PETSC_COMM_WORLD,&b);
		VecSetSizes(b,PETSC_DECIDE,nfnodes);
		VecSetFromOptions(b);

		VecDuplicate(b,&x);
		VecDuplicate(b,&dirichletVec);
		VecDuplicate(b,&SST);

		// setup rhs vector
		VecZeroEntries(dirichletVec);
		VecSetValues(dirichletVec,nnfnodes,pMData->getDirichlet_idx(),pMData->getDirichlet_data(),INSERT_VALUES);
		setSST();

		// do not initialize this variable again
		initialize = true;
	}

	void MEBFV_elliptic::setSST(){
		int nwells, nnodes, i, j, row;
		double Qi;
		nwells = pPPData->getNumWells();			// number of wells prescribed (Neumann boundary value)
		for (i=0; i<nwells; i++){
			nnodes = pPPData->getNumWellNodes(i);	// number of nodes located on well
			for (j=0; j<nnodes; j++){
				pPPData->getFlowRate(i,j,row,Qi);
				VecSetValue(SST,row,Qi,ADD_VALUES);
			}
		}
		VecAssemblyBegin(SST);
		VecAssemblyEnd(SST);
	}

	void MEBFV_elliptic::freeMemory(){
		MatDestroy(A);
		MatDestroy(dirichletMat);
		VecDestroy(b);
		VecDestroy(x);
		VecDestroy(SST);
		VecDestroy(dirichletVec);
	}
}
