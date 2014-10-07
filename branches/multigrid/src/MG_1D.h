#ifndef _MG1D_H
#define _MG1D_H

class MG_1D{
    private :
        PetscErrorCode ierr;
        Mat R,I; // m. R-Restriction(n_cXn) m. I-Interpolation(nXn_c) AR-A*R(n_cxn)

    public :
        Mat A;//A-coarse discretization
        Vec e;
        Vec r;   //e-coarse error(must be smoothed) r-residual(rhs)
        int n_coarse;//use for restriction and interpolation ...basis of a new coarse grid
        int n_fine;//use for restriction and interpolation
        int printMatrixToFile(Mat &A,const char *filename);
        MG_1D(){};
        ~MG_1D();
        int createComponents(int *n,int gridlevel);
        int assembleRAI(int n,Mat A_fine);
        void Restrict(Vec r_fine,int i);        //restriction proc.(from fine to coarse)
        void Interpolate(Vec e_fine,int i);     //interpolation proc.(from coarse to fine)
        int solver();
        int printVectorToFile(Vec& v,const char* filename);


};

#endif  //so..if MG_1D.h is already included, ifndef will jump to endif,avoiding double inclusion
