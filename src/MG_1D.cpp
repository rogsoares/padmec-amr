#include "includes.h"
#include "libincludes.h"



int MG_1D::createComponents(int *n,int gridlevel){

    bool evendetected=true;
    n_fine=*n;

        while(evendetected==true){
            if(*n%2==0){
                            *n-=1;
            }
                else{   //add grid
                    n_coarse=(*n-1)/2;
                    evendetected=false;
                }
        }
            *n=n_coarse; //going up ...

    //cout<<"\nGrid "<<gridlevel<<":\n"<<"e-"<<n_coarse<<"x1"<<"\t Acoarse-"<<n_coarse<<"x"<<n_coarse<<endl<<endl; //for test purposes

    //create error vec
    ierr=VecCreate(PETSC_COMM_WORLD,&e);CHKERRQ(ierr);
    ierr=VecSetSizes(e,PETSC_DECIDE,n_coarse);CHKERRQ(ierr);
    ierr=VecSetFromOptions(e);CHKERRQ(ierr);
    VecSet(e,0); //initial value for e
    ierr=VecDuplicate(e,&r);CHKERRQ(ierr);
    //A_coarse
    ierr=MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr=MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n_coarse,n_coarse);CHKERRQ(ierr);
    ierr=MatSetFromOptions(A);CHKERRQ(ierr);
    ierr=MatSetUp(A);CHKERRQ(ierr);
    //Restriction       //alternative way ->avarage between consecutive points   //Mat Create AIJ?
    ierr=MatCreate(PETSC_COMM_WORLD,&R);CHKERRQ(ierr);
    ierr=MatSetSizes(R,PETSC_DECIDE,PETSC_DECIDE,n_coarse,n_fine);CHKERRQ(ierr);
    ierr=MatSetFromOptions(R);CHKERRQ(ierr);
    ierr=MatSetUp(R);CHKERRQ(ierr);
    //Interpolation
    ierr=MatCreate(PETSC_COMM_WORLD,&I);CHKERRQ(ierr);  //Mat Create AIJ?
    ierr=MatSetSizes(I,PETSC_DECIDE,PETSC_DECIDE,n_fine,n_coarse);CHKERRQ(ierr);
    ierr=MatSetFromOptions(I);CHKERRQ(ierr);
    ierr=MatSetUp(I);CHKERRQ(ierr);

    //cout<<"\nGrid "<<gridlevel<<"'s components...ok\n\n";
    return 0;

}

int MG_1D::assembleRAI(int n,Mat A_fine){

   // cout<<"\n\nAssembling RAI...\n\n";

    PetscScalar value[3];
    PetscInt col[3];
    int nonzeros;//MatGetRow parameter
    const PetscInt* colsIndex;//MatGetRow parameter
    const PetscScalar* values;//MatGetRow parameter
    int correction=0;   //used for R assembling
    int i;

    //R
    //cout<<"\nR...";
    value[0] = 1.0/4.0 ; value[1] =2.0/4.0; value[2] =1.0/4.0;
        for (i=0;i<n_coarse;i++){
            col[0] = correction;++correction;
            col[1] = correction;++correction;
            col[2] = correction;
            ierr = MatSetValues(R,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
        }
    ierr = MatAssemblyBegin(R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(R,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    //Transpose R to form I
    //cout<<"I...";
        for(int i=0;i<n_coarse;i++){
            MatGetRow(R,i,&nonzeros,&colsIndex,&values);//get R row's i nonzeros
            MatSetValues(I,nonzeros,colsIndex,1,&i,values,INSERT_VALUES); //we take R's i row,and transform into a I column i
        }
        for(int i=0;i<n_coarse;i++){
            MatRestoreRow(R,i,&nonzeros,&colsIndex,&values);
        }
    ierr = MatAssemblyBegin(I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    MatScale(I,2.0);

    //A_coarse=RAI
    //cout<<"Acoarse...";
    MatMatMatMult(R,A_fine,I,MAT_INITIAL_MATRIX ,PETSC_DEFAULT,&A);
    ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    /*
    printMatrixToFile(R,"Restrictionmat");
    printMatrixToFile(I,"Interpolationmat");
    printMatrixToFile(A,"Coarsediscretmat");*/
    return 0;
}

void MG_1D::Restrict(Vec r_fine,int i){

       // if(i==0)
         //   cout<<"\n->restrict from  base grid to level "<<i; //for test purposes
            //    else
             //       cout<<"\n->restrict from level "<<i-1<<" to "<<i; //for test purposes
    MatMult(R,r_fine,r);
    //cout<<"...ok\n";
        printVectorToFile(r,"r_restricted");

}

void MG_1D::Interpolate(Vec e_fine,int i){

                //printVectorToFile(e,"e_coarse");
                //printVectorToFile(e_fine,"e_fine");


    //if(i==0)
      //  cout<<"\n->interpolate from level "<<i<<" to base grid"; //for test purposes ;//for test purposes
      //      else
        //        cout<<"\n->interpolate from level "<<i<<" to "<<i-1; //for test purposes ;//for test purposes
    MatMult(I,e,e_fine);
    //cout<<"...ok\n";
}


int MG_1D::solver(){

	KSP ksp;
    PetscErrorCode ierr;
	PC preconditioner;
	PetscInt its;
	int itnum=5; //iteration in a grid defined by the program
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCJACOBI);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,itnum);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,r,e);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
    //cout<<"\n\ncurrent iteration(smoothing -> e) : "<<its<<endl;
   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

        return 0;
}

MG_1D::~MG_1D(){
    //cout<<"\nDestroying components (from coarse grid)...\n";   //for test purposes
    VecDestroy(&e);
    VecDestroy(&r);
    MatDestroy(&A);
    MatDestroy(&R);
    MatDestroy(&I);
}

int MG_1D::printMatrixToFile(Mat &A,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(A,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

int MG_1D::printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}




