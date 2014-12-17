#include "includes.h"
#include "libincludes.h"
#include <sstream> //output file1.txt,file2.txt...


int grid1D::createComponents(int n,int gridlevel,Mat A_fine){

    int indexes[2]={0,0}; //used to build contribution vector (where a porcentage of the boundaries are)

    firstuse=true;
    n_fine=n;
    SetOperator op;   //build R and I

    //cout<<"\nGrid "<<gridlevel<<":\n"<<"e-"<<n_coarse<<"x1"<<"\t Acoarse-"<<n_coarse<<"x"<<n_coarse<<endl<<endl; //for test purposes

    //initialize iteration counter
    currentitnum=0;

    //create error vec
    ierr=VecCreate(PETSC_COMM_WORLD,&e);CHKERRQ(ierr);
    ierr=VecSetSizes(e,PETSC_DECIDE,n_coarse);CHKERRQ(ierr);
    ierr=VecSetFromOptions(e);CHKERRQ(ierr);
    VecSet(e,0); //e=0
    //create residual vec
    ierr=VecDuplicate(e,&r);CHKERRQ(ierr);
    //create and set up contribution vec
    ierr=VecCreate(PETSC_COMM_WORLD,&boundCtbt);CHKERRQ(ierr);
    ierr=VecSetSizes(boundCtbt,PETSC_DECIDE,n_fine);CHKERRQ(ierr);
    ierr=VecSetFromOptions(boundCtbt);CHKERRQ(ierr);
    VecSet(boundCtbt,0); //e=0
    //get Dirichlet conditions
    D0=op.getD0();
    DL=op.getDL();
    PetscScalar vals[2]={D0,DL};
    indexes[1]=n_fine-1;
    VecSetValues(boundCtbt,2,indexes,vals,INSERT_VALUES);
    VecAssemblyBegin(boundCtbt);VecAssemblyEnd(boundCtbt);

    //Restriction operator
    op.SetOperator_Restriction(n_fine,n_coarse,&R);
    //Interpolation operator
    op.SetOperator_Interpolation(n_fine,n_coarse,R,&I);
    //coarse discretization matrix A
    op.SetOperator_CoarseA(n_fine,n_coarse,R,A_fine,I,&A);
    // cout<<"\nGrid "<<gridlevel<<"'s components...ok\n\n";

    printMatrixToFile(R,"Restrictionmat");
    printMatrixToFile(I,"Interpolationmat");
    printMatrixToFile(A,"Coarsediscretmat");

    return 0;

}

void grid1D::Restrict(Vec r_fine,int i){
       // PetscReal ResNorm; //Multigrid Convergence Condition
   // VecNorm(r_fine,NORM_2,&ResNorm);
    //    cout<<"\nRESIDUALNORM_atMG_1D.cpp  before restriction(at level:"<<i<<")->"<<ResNorm<<endl;

        //if(i==0)
          //  cout<<"\n->restrict from  base grid to level "<<i+1; //for test purposes
            //    else
              //      cout<<"\n->restrict from level "<<i<<" to "<<i+1; //for test purposes
    MatMult(R,r_fine,r);
    //cout<<"...ok\n";

       // VecNorm(r,NORM_2,&ResNorm);
       // cout<<"\nRESIDUALNORM_atMG_1D.cpp  after restriction(at level:"<<i+1<<")->"<<ResNorm<<endl;


    //printVectorToFile(r,thisgoesout);
    stringstream ss;
    ss<<"coarse/OutputR_ucrs"<<i<<".txt";
    printVectorToFile(r,ss.str().c_str());
}

void grid1D::Interpolate(Vec e_fine,int i){
                //printVectorToFile(e,"e_coarse");
                //printVectorToFile(e_fine,"e_fine");


    //if(i==0)
      //  cout<<"\n->interpolate from level "<<i+1<<" to base grid"; //for test purposes ;//for test purposes
        //    else
          //      cout<<"\n->interpolate from level "<<i+1<<" to "<<i; //for test purposes ;//for test purposes

       // PetscReal ResNorm; //Multigrid Convergence Condition
       // VecNorm(e,NORM_2,&ResNorm);
        //cout<<"\nERRORNORM_atMG_1D.cpp  before interpolation(at level:"<<i+1<<")->"<<ResNorm<<endl;

    MatMult(I,r,e_fine);
    VecAXPY(e_fine,0.5,boundCtbt); //adding boundary contributions
      //  VecNorm(e_fine,NORM_2,&ResNorm);
    //    cout<<"\nERRORNORM_atMG_1D.cpp  after interpolation(at level:"<<i<<")->"<<ResNorm<<endl;
    //cout<<"...ok\n";

	 stringstream ss;
    ss<<"coarse/OutputI_ucrs"<<i+1<<".txt";
    printVectorToFile(e_fine,ss.str().c_str());
}


int grid1D::solver(){
        currentitnum+=15; //iteration variation...(number of iterations per cycle)

	KSP ksp;
    PetscErrorCode ierr;
	PC preconditioner;
	PetscInt its;
	  /*PetscReal ResNorm; //Multigrid Convergence Condition
        VecNorm(r,NORM_2,&ResNorm);
        cout<<"\nRESIDUALNORM_atMG_1D.cpp  before error relaxation ->"<<ResNorm<<endl;
        VecNorm(e,NORM_2,&ResNorm);
        cout<<"\nERRONORM_atMG_1D.cpp  before error relaxation ->"<<ResNorm<<endl;*/
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCJACOBI);CHKERRQ(ierr);
      if(firstuse==true){
                ierr = KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);
                firstuse=false;
            }
                    else
                        ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,currentitnum);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,r,e);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
    //cout<<"\n\ncurrent iteration(smoothing -> e) : "<<its<<endl;
   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    //VecNorm(e,NORM_2,&ResNorm);
        //cout<<"\nERRORNORM_atMG_1D.cpp  after  relaxation ->"<<ResNorm<<endl;

        return 0;
}


grid1D::~grid1D(){
   // cout<<"\nDestroying components (from coarse grid)...\n";   //for test purposes
    VecDestroy(&e);
    VecDestroy(&r);
    VecDestroy(&boundCtbt);
    MatDestroy(&A);
    MatDestroy(&R);
    MatDestroy(&I);
}

int grid1D::printMatrixToFile(Mat &A,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(A,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

int grid1D::printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}




