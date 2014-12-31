/*
 * residualCalculator.cpp
 *
 *  Created on: 27/12/2014
 *      Author: Julio Cezar
 */

#include "includes.h"

residualCalculator::residualCalculator(int n,Vec *y,Mat *A){

    //create A*u Vector
    VecCreate(PETSC_COMM_WORLD,&Au);
    VecSetSizes(Au,PETSC_DECIDE,n);
    VecSetFromOptions(Au);

    //copy RHS for the member field b
    VecDuplicate(*y,&b);
    VecCopy(*y,b);

    //copy original MATRIX for the member field matrix
    MatDuplicate(*A,MAT_COPY_VALUES,&matrix);

}

void residualCalculator::Calculate(Vec *r,Vec u){
    //calculate the residual
    VecCopy(b,*r);
    MatMult(matrix,u,Au);
    VecAXPY(*r,-1,Au); //r=y-Au
    VecNorm(*r,NORM_2,&rNorm);
    cout<<"\nRESIDUALNORM ->"<<rNorm<<endl;
}

PetscReal residualCalculator::get_residualNorm(){
//get the latest residual norm
    return rNorm;
}

int residualCalculator::printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}


