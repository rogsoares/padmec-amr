#include "includes.h"
#include "libincludes.h"

Controller_V::Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int*counter){

    MAXGRIDS=0;
    *counter=0;//counts multigrid iterations
    grid1D* pGrid; //pointer to  grid
    gridGenerator gridOpt(1);
    PetscInt its=0;//,oldits=0;
    PetscReal eNorm; //Multigrid Convergence Condition
    Vec Au; //residual factor
    VecCreate(PETSC_COMM_WORLD,&Au);
    VecSetSizes(Au,PETSC_DECIDE,n);
    VecSetFromOptions(Au);
    VecSet(e,0);
    firstuse=true; //alternates :PETSC INTIAL GUESS NON-ZERO

    //VecNorm(e,NORM_2,&ResNorm);
    //cout<<"\nERRORNORM_atController_V.cppbefore Multigrid->"<<ResNorm<<endl;
    currentitnum=0; //iteration counter initialization
    double CPUTIME=MPI_Wtime();

    pGrid=gridOpt.generateGrids(n,&MAXGRIDS); //allocate memmory
    gridOpt.getCoarsePoints(pGrid);//get generated points...

    //pGrid=generateGrids(n,&MAXGRIDS);//return an array of grids
         for(int i=0;i<MAXGRIDS;i++){

               if(i>0)
                    pGrid[i].createComponents(pGrid[i-1].n_coarse,i+1,pGrid[i-1].A);
                        else
                            pGrid[i].createComponents(n,i+1,A);

         }

        //do{

            its=0;
            //smooth the error
            KSPJacobi(A,y,u,its);
            //calculate the residual
            VecCopy(y,r);
            MatMult(A,u,Au);
            VecAXPY(r,-1,Au); //r=y-Au
      //         VecNorm(r,NORM_2,&ResNorm);
    //cout<<"\nRESIDUALNORM_atController_V.cppbefore Multigrid->"<<ResNorm<<endl;

            Multigrid(n,A,u,e,pGrid,MAXGRIDS);

            VecAXPY(u,1.0,e);//u'=u+e
            //VecNorm(r,NORM_2,&ResNorm);
            //cout<<"\nRESIDUALNORM_atController_V.cppafter Multigrid->"<<ResNorm<<endl;
            VecNorm(e,NORM_2,&eNorm);
            //cout<<"\nERRORNORM_atController_V.cppafter Multigrid->"<<ResNorm<<endl;
           // cout<<endl<<endl<<"\t\t====end of cycle "<<*counter+1<<"===="<<endl<<endl<<endl<<endl;
            ++*counter;

//      }while(eNorm>1);

        delete[] pGrid;
        cout<<"\n = = = = Multigrid Type : V-cycle =  =  =  =\n\n";
        cout<<"\nSmoother : Jacobi(PRECONDITIONING)+RICHARDSON(KSP)\n";
        cout<<"\nNumber of coarse grids used : "<<MAXGRIDS<<endl;
        cout<<"\nTotal number of cycles : "<<*counter<<endl;
        cout<<"\nCPU Time : "<<MPI_Wtime() - CPUTIME<<endl<<endl;
}



int Controller_V::Multigrid(int n,Mat A_fine,Vec r_fine,Vec e_fine,grid1D* pGrid,int MAXGRIDS){

///GO UPWARDS

        for (int i=0;i<MAXGRIDS;i++){

            //restrict r
                    if (i==0)
                        pGrid[i].Restrict(r_fine,i);
                            else
                                pGrid[i].Restrict(pGrid[i-1].r,i);

            //smooth the error
            pGrid[i].solver();

        }

///GO DOWNWARDS

        for (int i=MAXGRIDS-1;i>=0;i--){
            if(i>0){
                pGrid[i].Interpolate(pGrid[i-1].r,i);//test
                //pGrid[i].Interpolate(pGrid[i-1].e,i);//original

            }
                else{
                    pGrid[i].Interpolate(r_fine,i);//original
                    //pGrid[i].Interpolate(e_fine,i);//original
                    cout<<"\nend of multigrid...\n\n";
                                return 0;// check this

                }
        }

}


int Controller_V::KSPJacobi(Mat A, Vec y, Vec v, PetscInt &its){
    currentitnum+=10; //iteration variation...(number of iterations per cycle)
    //cout<<"\nCONTROLLER'S current itnum : "<<currentitnum;
	KSP ksp;
    PetscErrorCode ierr;
	PC preconditioner;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);                     //DIFFERENT_NON_ZERO_PATTERN is no longer a parameter
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
    ierr = KSPSolve(ksp,y,v);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
    //cout<<"\n\ncurrent iteration : "<<its<<endl;
   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);                      //variable adress is now required on KSPDDestroy Funcion

        return 0;
}

Controller_V::~Controller_V(){
    //cout<<"\n\nshutting down multigrid process...";
}

int Controller_V::printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}
