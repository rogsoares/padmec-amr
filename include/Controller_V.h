#ifndef _CONTROLLER_
#define _CONTROLLER_
#include "grid1D.h"

///V-CYCLE MG PATTERN

class Controller_V {
    //friend class MG_1D;
    public:
        int currentitnum; //number of iterations
        Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int* counter);//start multigrid
        ~Controller_V();
        int Multigrid(int n,Mat A_fine,Vec r_fine,Vec e_fine,grid1D* pGrid,int *i,int MAXGRIDS,Vec errorF);//erroF =original error
        int MAXGRIDS;
        grid1D* generateGrids(int n,int*MAXGRIDS);
        int KSPJacobi(Mat A, Vec y, Vec v, PetscInt &its);
        bool firstuse;
        int printVectorToFile(Vec& v,const char* filename);


};



#endif
