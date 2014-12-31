/*
 * Controller_V.cpp
 *
 *  Created on: 10/10/2014
 *      Author: Julio Cezar
 */

#include "includes.h"

Controller_V::Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int*counter){

    MAXGRIDS=0;
    *counter=0;//counts multigrid iterations
    PetscInt its=0;//,oldits=0;
    PetscReal rNorm; //residual norm
    ///--------------------------------------
    VecSet(e,0);
    firstuse=false; /// false ,because we already have an aproximation for u(in the future ,this must be a constructor's parameter)

    ///read some data...
    grid1D* pGrid; //pointer to  grid
    gridGenerator gridOpt(1); //int parameter : 1-read data and create a file containing specific grid points
    smoother solver;

    ///custom classes+tools initialization
    residualCalculator residual(n,&y,&A);

    currentitnum=0; //iteration counter initialization

    ///=======================GENERATE GRIDS=========================================
    pGrid=gridOpt.generateGrids(n,&MAXGRIDS); //allocate memmory
    gridOpt.getCoarsePoints(pGrid);//get generated points...
    //pGrid=generateGrids(n,&MAXGRIDS);//return an array of grids
         for(int i=0;i<MAXGRIDS;i++){

               if(i>0)
                    pGrid[i].createComponents(pGrid[i-1].n_coarse,i+1,pGrid[i-1].A);
                        else
                            pGrid[i].createComponents(n,i+1,A);

         }
    ///===============================================================================

    //START V-CYCLE

        double CPUTIME=MPI_Wtime();//start performance test

      //  do{

            its=0;
            cout<<"\npre-smoothing in fine grid\n";
            //smooth the error(pre-smoothing)
            solver.PreSmoothing(A,y,u,&firstuse);
             //calculate the residual
            residual.Calculate(&r,u);
            //restrict->solve(->restrict->solve->restrict->solve....)->interpolate->solve(->interpolate->solve->interpolate->solve->interpolate...)
            InterGridOperations(A,r,e,pGrid,MAXGRIDS,solver);
            //u'=u+e
            VecAXPY(u,1.0,e);
            //post-smoothing
            cout<<"\npost-smoothing in fine grid\n";
            solver.PostSmoothing(A,y,u);
            residual.Calculate(&r,u);

            cout<<endl<<endl<<"\t\t====end of cycle "<<*counter+1<<"===="<<endl<<endl<<endl<<endl;
            ++*counter;

     // }while(residual.get_residualNorm()>0.01);

        delete[] pGrid;
        cout<<"\n = = = = Multigrid Type : V-cycle =  =  =  =\n\n";
        cout<<"\nSmoother : Jacobi(PRECONDITIONING)+RICHARDSON(KSP)\n";
        cout<<"\nNumber of coarse grids used : "<<MAXGRIDS<<endl;
        cout<<"\nTotal number of cycles : "<<*counter<<endl;
        cout<<"\nCPU Time : "<<MPI_Wtime() - CPUTIME<<endl<<endl;
}



int Controller_V::InterGridOperations(Mat A_fine,Vec r_fine,Vec e_fine,grid1D* pGrid,int MAXGRIDS,smoother solver){

///GO UPWARDS     from fine grid to kth grid

        for (int i=0;i<MAXGRIDS;i++){
            //restrict r
                    if (i==0)
                        pGrid[i].Restrict(r_fine,i);
                            else
                                pGrid[i].Restrict(pGrid[i-1].r,i);
            //smooth the error
            solver.PreSmoothing(pGrid[i].A,pGrid[i].r,pGrid[i].e,&pGrid[i].firstuse);
        }

///GO DOWNWARDS      from kth grid to fine grid

        for (int i=MAXGRIDS-1;i>=0;i--){
            if(i>0){
            //smooth the error before interpolation,ensuring it'll have a good interpolation
                solver.PostSmoothing(pGrid[i].A,pGrid[i].r,pGrid[i].e);
                pGrid[i].Interpolate(pGrid[i-1].e,i);//interpolate

            }
                    else{
                        solver.PostSmoothing(pGrid[i].A,pGrid[i].r,pGrid[i].e);
                        pGrid[i].Interpolate(e_fine,i);//interpolate to the fine grid
                        cout<<"\nend of multigrid...\n\n";
                        return 0;// check this
                    }
        }

}

int Controller_V::printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

Controller_V::~Controller_V(){
    //cout<<"\n\nshutting down multigrid (V-cycle) process...";
}


