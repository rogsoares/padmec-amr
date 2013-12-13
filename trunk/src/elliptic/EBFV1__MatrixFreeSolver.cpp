/*
 * EBFV1__MatrixFreeSolver.cpp
 *
 *  Created on: 25/08/2012
 *      Author: rogsoares
 */

#include "EBFV1_elliptic.h"
#include <time.h>

namespace PRS           // PRS: Petroleum Reservoir Simulator
{

	double EBFV1_elliptic::setMatrixFreeOperation(pMesh theMesh){
	double startt = MPI_Wtime();
	PetscLogDouble flops1;
	PetscGetFlops(&flops1);
	//	cout << "Mflops: " << flops1/1.0e6 << endl;
	PetscInt m, n, its;
	ierr = VecGetOwnershipRange(output,&m,&n); CHKERRQ(ierr);
	PetscInt nLRows = n - m;
	PetscInt numGF = pMData->getNum_GF_Nodes();
	ierr = MatCreateShell(PETSC_COMM_WORLD,nLRows,nLRows,numGF,numGF,matvec_struct,&matrix);CHKERRQ(ierr);
	ierr = MatSetFromOptions(matrix); CHKERRQ(ierr);
	ierr = MatShellSetOperation(matrix, MATOP_MULT,(void(*)(void))&EBFV1_elliptic::MatMultUser); CHKERRQ(ierr);

	/*
	 * Solve system of equation using a matrix-free approach
	 */
	//PetscInt m, n;
//	MatGetSize(matvec_struct->G,&m,&n); printf("G  : %d x %d\n",m,n);
//	VecGetSize(matvec_struct->RHS,&m);  printf("RHS: %d\n",m);
//	VecGetSize(output,&m);              printf("Sol: %d\n",m);

	// create auxiliary vector: used during matrix-vector operations
	ierr = VecCreate(PETSC_COMM_WORLD,&matvec_struct->z); CHKERRQ(ierr);
	ierr = VecSetSizes(matvec_struct->z,PETSC_DECIDE,matvec_struct->F_nrows); CHKERRQ(ierr);
	ierr = VecSetFromOptions(matvec_struct->z); CHKERRQ(ierr);
	
	// use last pressure solution as guess solution for iterative solver and then reduce number of iterations
	static bool guess_sol = false;
	PetscTruth guessNonZero = PETSC_FALSE;
	if (guess_sol){
		int row = 0;
		guessNonZero = PETSC_TRUE;
		// get solution from mesh nodes and set output vector for guess solution
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			if ( pSimPar->isNodeFree( GEN_tag( node->getClassification() )) ){
				ierr = VecSetValue(output,row++,pPPData->getPressure(node),INSERT_VALUES); CHKERRQ(ierr);
			}
		}
		VIter_delete(vit);
	}
	guess_sol = true;

	// call PETSc to solver system of equation
	KSP_solver(matrix,matvec_struct->G,matvec_struct->RHS,output,pSimPar,guessNonZero,KSPBCGS,PCASM,its);
	cout << "KSP_Solver: Num. iter =  " << its << endl;
	
	//printVectorToFile(output,"SOLUTION.txt");
	//VecView(matvec_struct->RHS,PETSC_VIEWER_STDOUT_WORLD); //STOP();
	//MatView(matvec_struct->G,PETSC_VIEWER_STDOUT_WORLD); //STOP();
//	VecView(output,PETSC_VIEWER_STDOUT_WORLD);
//				PetscReal val;
//				VecNorm(output,NORM_2,&val);
//				PetscPrintf(PETSC_COMM_WORLD,"norm: %f",val);
//				STOP();

	PetscLogDouble flops2;
	PetscGetFlops(&flops2);
	//	cout << "Mflops: " << flops2/1.0e6 << endl;
	//cout << "Mflops: " << (flops2-flops1)/1.0e6 << endl;


	double endt = MPI_Wtime();
	//cout << "Elaprse-time: " << endt-startt << endl;
	return endt-startt;
}
}
