#include "includes.h"
#include "libincludes.h"

Controller_V::Controller_V(int n,Vec r,Vec e,Vec u,Vec y,Mat A,int*counter){

    int level=0; //grid level
    int n_temp=n; //this variable is used in grid const ruction
    MAXGRIDS=0;
    MG_1D* pGrid; //pointer to  grid
    PetscInt its,oldits;
    Vec Ap; //residual factor
    VecCreate(PETSC_COMM_WORLD,&Ap);
    VecSetSizes(Ap,PETSC_DECIDE,n);
    VecSetFromOptions(Ap);

    pGrid=generateGrids(n,&MAXGRIDS);//return an array of grids
         for(int i=0;i<MAXGRIDS;i++){
            pGrid[i].createComponents(&n_temp,i+1);
         }

        do{
            oldits=its;            //its remains static ,when the convergence condition is reached

            //smooth the error
            KSPJacobi(A,y,u,its,3);

            //calculate the residual
            VecCopy(y,r);
            MatMult(A,u,Ap);
            VecAXPY(r,-1,Ap); //r=y-Ap

          /* recursive function here...*/
            Multigrid(n,A,r,e,pGrid,&level,MAXGRIDS,e);

            VecAXPY(u,1.0,e);//u'=u+e
            ++*counter;

        }while(oldits!=its);

            for(int i=(MAXGRIDS-1);i--;i>=0){
                    pGrid[i].~MG_1D();
            }
        delete[] pGrid;
}

MG_1D* Controller_V::generateGrids(int n,int*MAXGRIDS){
    MG_1D* newgrids;

   // cout<<"\n\ncreating grids...\n\n";
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
    newgrids=new MG_1D [*MAXGRIDS];
    return newgrids;

}

int Controller_V::Multigrid(int n,Mat A_fine,Vec r_fine,Vec e_fine,MG_1D* pGrid,int *i,int MAXGRIDS,Vec errorF){

    /* recursive function here...*/

    if(n>=5){
        /*set restriction,interpolation and coarse grid matrices
         A(coarse)=RAI // A(coarse) is an aproximation of the original matrix from the fine grid
         I=R transposed
        */
        pGrid[*i].assembleRAI(n,A_fine);
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
                pGrid[*i].Interpolate(pGrid[*i-1].e,*i);
                return Multigrid(0,pGrid[*i-1].A,pGrid[*i-1].r,pGrid[*i-1].e,pGrid,i,MAXGRIDS,errorF);
                /*n is considered to be zero because we jump "if" at the top (of Multigrid function)
                        and we are only using "else" for interpoltation means...(only grids_1d and i are needed)
                */

            }
                else{
                    pGrid[*i].Interpolate(errorF,*i);
                   // cout<<"\nend of multigrid...\n\n";
                                return 0;

                }
        }


}


int Controller_V::KSPJacobi(Mat A, Vec y, Vec v, PetscInt &its,int itnum){
	KSP ksp;
    PetscErrorCode ierr;
	PC preconditioner;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);                     //DIFFERENT_NON_ZERO_PATTERN is no longer a parameter
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCJACOBI);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,itnum);CHKERRQ(ierr);
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
