/*
 * GRID (1D TYPE) CLASS
 *
 * grid1D.h
 *
 *  Created on: 16/12/2014
 *      Author: Julio Cezar
 */


#ifndef _G1D_H
#define _G1D_H
#include <string.h>


///->obs if there's any problem when adding a variable to this class,compile without the destructor,then compile again with the destructor

class grid1D{
    private :
        PetscErrorCode ierr;
        double D0,DL; //boundaries
        Mat R,I; // m. R-Restriction(n_cXn) m. I-Interpolation(nXn_c) AR-A*R(n_cxn
    public :
        grid1D(){};
        ~grid1D();
        Mat A;//A-coarse discretization
        Vec e;//coarse errors (must be smoothed)
        Vec r;  // r-residual(rhs)
        Vec boundCtbt;
        int n_coarse;//coarse grid's internal points
        int n_fine;//fine grid's internal points
        int currentitnum;
        bool firstuse;
        int createComponents(int n,int gridlevel,Mat A_fine);
        void Restrict(Vec r_fine,int i);        //restriction proc.(from fine to coarse)
        void Interpolate(Vec e_fine,int i);     //interpolation proc.(from coarse to fine)
        int printVectorToFile(Vec& v,const char* filename);
        int printMatrixToFile(Mat &A,const char *filename);

};

#endif  //so..if G1D.h is already included, ifndef will jump to endif,avoiding double inclusion
