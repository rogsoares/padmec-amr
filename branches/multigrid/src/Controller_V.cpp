#include "includes.h"
#include "libincludes.h"

Controller_V::Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int*counter){

    int level=0; //grid level
    int n_temp=n; //this variable is used in grid construction
    MAXGRIDS=0;
    *counter=0;//counts multigrid iterations
    grid1D* pGrid; //pointer to  grid
    PetscInt its=0,oldits=0;
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
    pGrid=generateGrids(n,&MAXGRIDS);//return an array of grids
         for(int i=0;i<MAXGRIDS;i++){

               if(i>0)
                    pGrid[i].createComponents(&n_temp,i+1,pGrid[i-1].A);
                        else
                            pGrid[i].createComponents(&n_temp,i+1,A);



            /*set restriction,interpolation and coarse grid matrices
             A(coarse)=RAI // A(coarse) is an aproximation of the original matrix from the fine grid
             I=R transposed
            */
                //if(i==0)
                  //  pGrid[i].assembleRAI(n,A);   //if this is the first jump,then take original data
                  //      else
                    //        pGrid[i].assembleRAI(pGrid[i-1].n_coarse,pGrid[i-1].A);//if this is jump n,then take past coarse grid's data

         }

        //do{

            its=0;
            //smooth the error
            KSPJacobi(A,y,u,its);
             printVectorToFile(u,"coarse/true_u2");
            //calculate the residual
            VecCopy(y,r);
            MatMult(A,u,Au);
            VecAXPY(r,-1,Au); //r=y-Au
      //         VecNorm(r,NORM_2,&ResNorm);
    //cout<<"\nRESIDUALNORM_atController_V.cppbefore Multigrid->"<<ResNorm<<endl;

          /* recursive function here...*/
            Multigrid(n,A,u,e,pGrid,&level,MAXGRIDS,u);

            VecAXPY(u,1.0,e);//u'=u+e
            //VecNorm(r,NORM_2,&ResNorm);
            //cout<<"\nRESIDUALNORM_atController_V.cppafter Multigrid->"<<ResNorm<<endl;
            VecNorm(e,NORM_2,&eNorm);
            //cout<<"\nERRORNORM_atController_V.cppafter Multigrid->"<<ResNorm<<endl;
           // cout<<endl<<endl<<"\t\t====end of cycle "<<*counter+1<<"===="<<endl<<endl<<endl<<endl;
            ++*counter;

//      }while(eNorm>1);

            for(int i=(MAXGRIDS-1);i--;i>=0){
                    pGrid[i].~grid1D();
            }
        delete[] pGrid;
        cout<<"\n = = = = Multigrid Type : V-cycle =  =  =  =\n\n";
        cout<<"\nSmoother : Jacobi(PRECONDITIONING)+RICHARDSON(KSP)\n";
        cout<<"\nNumber of coarse grids used : "<<MAXGRIDS<<endl;
        cout<<"\nTotal number of cycles : "<<*counter<<endl;
        cout<<"\nCPU Time : "<<MPI_Wtime() - CPUTIME<<endl<<endl;
}

grid1D* Controller_V::generateGrids(int n,int*MAXGRIDS){
    grid1D* newgrids;

    //cout<<"\n\ncreating grids...\n\n";
                   while(n>=5){        //think about the last level matrix A=R(2x5)A(5x5)I(5x2)=A(2x2)->if n is smaller ,then the coarse_matrix will become a vector
                if(n%2==0){
                                n-=1;
                }
                    else{   //add grid
                        n=(n-1)/2;
                        ++*MAXGRIDS;
                      //  cout<<"\ngrid "<<*MAXGRIDS<<" created("<<n<<" points) :: level->"<<*MAXGRIDS-1<<"\n";    //for test purposes
                    }
            }
    newgrids=new grid1D[*MAXGRIDS];
    return newgrids;

}

int Controller_V::Multigrid(int n,Mat A_fine,Vec r_fine,Vec e_fine,grid1D* pGrid,int *i,int MAXGRIDS,Vec errorF){

    /* recursive function here...*/

    if(n>=5){

        //restrict r
        pGrid[*i].Restrict(r_fine,*i);
        //smooth the error
        pGrid[*i].solver();
        ++*i; ///next grid level
                return Multigrid(pGrid[*i-1].n_coarse,pGrid[*i-1].A,pGrid[*i-1].r,pGrid[*i-1].e,pGrid,i,MAXGRIDS,errorF);
    }

        else {
            --*i;  ///if the expansion is cancelled,then ++*i(called before "return Multigrid()...") must come back to *i...
            if(*i>0){
                pGrid[*i].Interpolate(pGrid[*i-1].r,*i);//test
                //pGrid[*i].Interpolate(pGrid[*i-1].e,*i);//original
                return Multigrid(0,pGrid[*i-1].A,pGrid[*i-1].r,pGrid[*i-1].e,pGrid,i,MAXGRIDS,errorF);

            }
                else{
                    pGrid[*i].Interpolate(errorF,*i);
                   // cout<<"\nend of multigrid...\n\n";
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
