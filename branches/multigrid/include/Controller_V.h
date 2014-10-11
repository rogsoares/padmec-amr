#ifndef _CONTROLLER_
#define _CONTROLLER_
#include "MG_1D.h"

///V-CYCLE MG PATTERN

class Controller_V {
    friend class MG_1D;
    public:
        Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int* counter);//start multigrid
        ~Controller_V();
        int Multigrid(int n,Mat A_fine,Vec r_fine,Vec e_fine,MG_1D* pGrid,int *i,int MAXGRIDS,Vec errorF);//erroF =original error
        int MAXGRIDS;
        MG_1D* generateGrids(int n,int*MAXGRIDS);
        int KSPJacobi(Mat A, Vec y, Vec v, PetscInt &its,int itnum);


};



#endif
