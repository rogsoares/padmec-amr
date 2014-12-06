#include "includes.h"
#include "libincludes.h"
#include <sstream> //output file1.txt,file2.txt...

SetOperator::SetOperator(){

    ifstream bound("parameters/boundaries.dat", ios::binary);

        if (bound) {

        double value; //store boundaries values and passes them to the variable D0,DL
        char buffer[13] ; //offset 13 columns("boundaries = ") in the file,then get the value

        bound.read (buffer,13);// read data as a block
        bound>>value;
        D0=value;value=0; //get first condition
        bound.read (buffer,13);// read data as a block:
        bound>>value;
        DL=value; //get second condition

        bound.close();

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

    PetscScalar value[3];
    PetscInt col[3];
    int correction=0;   //used for R assembling

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
    return 0;
}

int SetOperator::SetOperator_Interpolation(int n_fine,int n_coarse,Mat R,Mat* I){

    int nonzeros;//MatGetRow parameter
    const PetscInt* colsIndex;//MatGetRow parameter
    const PetscScalar* values;//MatGetRow parameter

    ierr=MatCreate(PETSC_COMM_WORLD,I);CHKERRQ(ierr);  //Mat Create AIJ?
    ierr=MatSetSizes(*I,PETSC_DECIDE,PETSC_DECIDE,n_fine,n_coarse);CHKERRQ(ierr);
    ierr=MatSetFromOptions(*I);CHKERRQ(ierr);
    ierr=MatSetUp(*I);CHKERRQ(ierr);

    //Transpose R
    //cout<<"I...";
        for(int i=0;i<n_coarse;i++){
            MatGetRow(R,i,&nonzeros,&colsIndex,&values);//get R row's i nonzeros
            MatSetValues(*I,nonzeros,colsIndex,1,&i,values,INSERT_VALUES); //we take R's i row,and transform into a I column i
        }
        for(int i=0;i<n_coarse;i++){
            MatRestoreRow(R,i,&nonzeros,&colsIndex,&values);
        }
    ierr = MatAssemblyBegin(*I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(*I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    MatScale(*I,2.0);

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


