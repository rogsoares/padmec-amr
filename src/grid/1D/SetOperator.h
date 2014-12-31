/*
 * OPERATOR'S CLASS (ASSEMBLY+SETUP)
 *
 * SetOperator.h
 *
 *  Created on: 14/12/2014
 *      Author: Julio Cezar
 */

#ifndef _SETOP_H
#define _SETOP_H

enum R_operator {INJECTION,FW};


class SetOperator{
    private :
        PetscErrorCode ierr;
         //**Dirichlet boundaries**
        double D0;
        double DL;
        //injection=0 full-weighting=1
        R_operator Rop;

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
