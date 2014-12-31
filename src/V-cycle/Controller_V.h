
/*
 * V-CYCLE  CLASS
 *
 * Controller_V.h
 *
 *  Created on: 10/10/2014
 *      Author: Julio Cezar
 */

#ifndef _CONTROLLER_
#define _CONTROLLER_
#include "../grid/1D/grid1D.h"
#include "../smoother/smoother.h"


///V-CYCLE MG PATTERN

class Controller_V {
    public:
        int currentitnum; //number of iterations
        Controller_V(int n,Vec r,Vec e,Vec y,Vec u,Mat A,int* counter);//set up and launch multigrid
        ~Controller_V();
        int InterGridOperations(Mat A_fine,Vec r_fine,Vec e_fine,grid1D* pGrid,int MAXGRIDS,smoother solver);//restriction+presmoothing+interpolation+postsmoothing
        int MAXGRIDS;
        bool firstuse;
        int printVectorToFile(Vec& v,const char* filename);

};

#endif
