/*
 * smoother.cpp
 *
 *  Created on: 20/12/2014
 *      Author: Julio Cezar
 */

#include "includes.h"

smoother::smoother(){

int value; //negligible variable

     ifstream opt("MG_parameters/options.dat", ios::binary);

        if (opt) {

            char buffer[32] ; //offsets columns

                    opt.read (buffer,11);
                    opt>>value;
                    opt.read (buffer,23);
                    opt>>value;
                    opt.read (buffer,9);
                    opt>>PreIts;
                    opt.read (buffer,13); //why is 13?..I don't know,the buffer behaves strangely with respect to the last parameter,
                                          //which has 10 columns ("PostIts = ")
                    opt>>PostIts;

            opt.close();

        }

}

int smoother::PreSmoothing(Mat A, Vec y, Vec u,bool* firstuse){

    //cout<<"\n\npre smoothing...";
	KSP ksp;
    PetscInt its;
    PetscErrorCode ierr;
	PC preconditioner;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);                     //DIFFERENT_NON_ZERO_PATTERN is no longer a parameter
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCJACOBI);CHKERRQ(ierr);

        if(*firstuse==true){
            ierr = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);
            *firstuse=false;
      //      cout<<"first use...";
        }
                else
                    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PreIts);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,y,u);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
    cout<<"\n\nPre Smoothing current iteration : "<<its<<endl;
   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);                      //variable adress is now required on KSPDDestroy Funcion

        return 0;
}

int smoother::PostSmoothing(Mat A, Vec y, Vec u){
 //cout<<"\n\npost smoothing...";
	KSP ksp;
    PetscInt its;
    PetscErrorCode ierr;
	PC preconditioner;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);                     //DIFFERENT_NON_ZERO_PATTERN is no longer a parameter
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCJACOBI);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PostIts);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,y,u);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
    cout<<"\n\nPost smoothing current iteration : "<<its<<endl;
    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);                      //variable adress is now required on KSPDDestroy Funcion

        return 0;

}
