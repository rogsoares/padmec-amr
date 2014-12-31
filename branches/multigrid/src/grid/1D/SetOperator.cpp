/*
 * SetOperator.cpp
 *
 *  Created on: 14/12/2014
 *      Author: Julio Cezar
 */

#include "includes.h"
#include <sstream> //output file1.txt,file2.txt...


SetOperator::SetOperator(){

//get boundary conditions
    ifstream bound("MG_parameters/boundaries.dat", ios::binary);
    char buffer[23]; //offset 13 columns("boundaries = ") in the file,then get the value
    int value;//used to get Restriction type

        if (bound) {

            bound.read (buffer,13);// read data as a block
            bound>>D0;//get first condition
            bound.read (buffer,13);// read data as a block:
            bound>>DL;//get second condition
            bound.close();

        }

//get restriction operator
    ifstream rest("MG_parameters/options.dat", ios::binary);

            if (rest) {

                rest.read (buffer,11);// read data as a block ...number is 15,because this line is trivial
                rest>>value;
                rest.read (buffer,23);// read data as a block:
                rest>>value;
                rest.close();

            }
            switch(value){
                case 0:
                    Rop=INJECTION;
                    break;
                case 1:
                    Rop=FW;
                    break;
                default :
                    cout<<"\n\nMG Error in MG_parameters/options.dat->RESTRICTION OPERATOR = "<<value<<"\n\n";
                    cout<<"use 0 (INJECTION) or 1 (FULL-WEIGHTING)\n\n";
                    exit(0);
            }

}

///SetOperator
            /*
            inputs:
                int n->current number of internal points of the nth grid / (R for SetOperator_Interpolation)
            output:
                Mat R/I ->nth grid (RESTRICTION/INTERPOLATION) operator

                    creates the restriction/interpolation operator R/I,then returns it back to a Grid object
            */

int SetOperator::SetOperator_Restriction(int n_fine,int n_coarse,Mat* R){
    int correction=0;   //used for R assembling
    PetscScalar value[3];
    PetscInt col[3];

        switch (Rop){

            case INJECTION :
                cout<<"\n using injection operator...";
                ierr=MatCreate(PETSC_COMM_WORLD,R);CHKERRQ(ierr);
                ierr=MatSetSizes(*R,PETSC_DECIDE,PETSC_DECIDE,n_coarse,n_fine);CHKERRQ(ierr);
                ierr=MatSetFromOptions(*R);CHKERRQ(ierr);
                ierr=MatSetUp(*R);CHKERRQ(ierr);

                //R
                //cout<<"\nR...";
                value[0] = 1.0 ; value[1]=0 ;value [2]=0;
                correction=2; //the first one goes in the second column
                col[2]=correction;
                    for (int i=0;i<n_coarse;i++){
                        col[0] = col[2]-1;++correction;
                        col[1] = correction;++correction;
                        col[2] = correction;
                        ierr = MatSetValues(*R,1,&i,1,col,value,INSERT_VALUES);CHKERRQ(ierr);
                    }
                ierr = MatAssemblyBegin(*R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = MatAssemblyEnd(*R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

                break;

            case FW :
                cout<<"\nusing FW operator...";
                ierr=MatCreate(PETSC_COMM_WORLD,R);CHKERRQ(ierr);
                ierr=MatSetSizes(*R,PETSC_DECIDE,PETSC_DECIDE,n_coarse,n_fine);CHKERRQ(ierr);
                ierr=MatSetFromOptions(*R);CHKERRQ(ierr);
                ierr=MatSetUp(*R);CHKERRQ(ierr);

                //R
                //cout<<"\nR...";
                value[0] = 1.0/4.0 ; value[1] =2.0/4.0; value[2] =1.0/4.0;
                    for (int i=0;i<n_coarse;i++){
                        col[0] = correction;++correction;
                        col[1] = correction;++correction;
                        col[2] = correction;
                        ierr = MatSetValues(*R,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
                    }
                ierr = MatAssemblyBegin(*R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = MatAssemblyEnd(*R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                break;
            }

            return 0;
}

int SetOperator::SetOperator_Interpolation(int n_fine,int n_coarse,Mat R,Mat* I){

    int correction=0;

    PetscScalar value[3];
    PetscInt line[3];

    ierr=MatCreate(PETSC_COMM_WORLD,I);CHKERRQ(ierr);  //Mat Create AIJ?
    ierr=MatSetSizes(*I,PETSC_DECIDE,PETSC_DECIDE,n_fine,n_coarse);CHKERRQ(ierr);
    ierr=MatSetFromOptions(*I);CHKERRQ(ierr);
    ierr=MatSetUp(*I);CHKERRQ(ierr);
    //cout<<"I...";

      value[0] = 1.0/2.0 ; value[1] =2.0/2.0; value[2] =1.0/2.0;
                    for (int i=0;i<n_coarse;i++){
                        line[0] = correction;++correction;
                        line[1] = correction;++correction;
                        line[2] = correction;
                        ierr = MatSetValues(*I,3,line,1,&i,value,INSERT_VALUES);CHKERRQ(ierr);
                    }
                ierr = MatAssemblyBegin(*I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
                ierr = MatAssemblyEnd(*I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    return 0;
}

///SetOperator_coarseA
        /*
        inputs:
            int n->current number of internal points of the nth grid,R,I operators
        output:
            Mat A ->nth grid

                assemble the coarse grid operator A
        */

int SetOperator::SetOperator_CoarseA(int n_coarse,int n_fine,Mat R,Mat A_fine,Mat I,Mat *A){
    ierr=MatCreate(PETSC_COMM_WORLD,A);CHKERRQ(ierr);
    ierr=MatSetSizes(*A,PETSC_DECIDE,PETSC_DECIDE,n_coarse,n_coarse);CHKERRQ(ierr);
    ierr=MatSetFromOptions(*A);CHKERRQ(ierr);
    ierr=MatSetUp(*A);CHKERRQ(ierr);
    //A_coarse=RAI
    //cout<<"Acoarse...";
    MatMatMatMult(R,A_fine,I,MAT_INITIAL_MATRIX ,PETSC_DEFAULT,A);
    ierr = MatAssemblyBegin(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    return 0;
}

//for test purposes
int SetOperator::printMatrixToFile(Mat &A,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(A,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}


