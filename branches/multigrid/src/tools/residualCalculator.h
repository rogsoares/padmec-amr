/*
 * RESIDUAL CALCULATOR  CLASS
 *
 * residualCalcualtor.h
 *
 *  Created on: 27/12/2014
 *      Author: Julio Cezar
 */

#ifndef _RESCALC_
#define _RESCALC_

class residualCalculator{
    private :
        Vec Au; //A(matrix)*u(aproximate value of the solution - variable)
        Vec b; //RHS of the problem
        Mat matrix; //original matrix of the problem
        PetscReal rNorm;

    public :
        residualCalculator(int n,Vec *y,Mat *A);
        void Calculate(Vec *r,Vec u);
        PetscReal get_residualNorm();
        int printVectorToFile(Vec& v,const char* filename);
};


#endif // _RESCALC_


