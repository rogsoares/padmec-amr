#ifndef _SETOP_H
#define _SETOP_H

class SetOperator{
    private :
        PetscErrorCode ierr;
         //**Dirichlet boundaries**
        double D0;
        double DL;
        //************************
    public :
        SetOperator();
        double getD0(){return D0;}
        double getDL(){return DL;}
        int SetOperator_Restriction(int n_fine,int n_coarse,Mat *R);
        int SetOperator_Interpolation(int n_fine,int n_coarse,Mat R,Mat *I);
        int SetOperator_CoarseA(int n_fine,int n_coarse,Mat R,Mat A_fine,Mat I,Mat *A);
        int printMatrixToFile(Mat &A,const char* filename);
};
#endif
