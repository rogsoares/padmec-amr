#include "includes.h"
/*
solving example 5.11 from Ertkein Turcay book...
    Solve the problem described in Example 5.8 using the implicit backward-difference formulation.Use a timestep size of 15 days.
*/
static char help[] = "Solve the problem described in Example 5.8 using the implicit backward-difference formulation.Use a timestep size of 15 days..\n\n";


int printMatrixToFile(Mat& A,const char* filename);
int printVectorToFile(Vec& v,const char* filename);
int KSP_linearsolver(Mat A, Vec y, Vec v, PetscInt &its,int itnum);
int KSP_Jacobi(Mat A, Vec y, Vec v, PetscInt &its,int itnum);
MG_1D* creategrids(int n,int*MAXGRIDS);
void deleteCoarseGrids(MG_1D* grids,int MAXGRIDS);

int main(int argc,char **args){

    PetscMPIInt size;
    PetscInt n,i,N;
    PetscInt its,oldits;
    PetscScalar value[3];
    PetscScalar neg_one = -1.0,one = 1.0;
    PetscInt col[3];
    PetscErrorCode ierr;

    PetscInitialize(&argc,&args,(char *)0,help);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);


///               ***********GRID TYPE > BLOCK-CENTERED *************
/*

  GB1   GB2        GBi-1  GBi  GBi+1     GBn-1  GBn        -->Flow Direction
_____________     ___________________    _____________
|  .  |  .  | ... |  .  |  .  |  .  |....|  .  |  .  |
|     |     |     |     |     |     |    |     |     |
-------------     -------------------    -------------
---X1----X2--------Xi-1----Xi---Xi+1-------Xn-1---Xn------------->x    (Lx=DeltaX1+DeltaX2....+DeltaXn)
X1-1/2    X1+3/2 Xi-3/2     Xi+1/2      Xn-3/2      Xn+1/2             DeltaX=Xi+1/2-Xi-1/2
    X1+1/2            Xi-1/2     Xi+3/2       Xn-1/2

*as we have mid-terms(xi-1/2,xi+1/2) we need more than n discrete points,we need N points

 \\ N= 1(first point)+ 2n(points+boundaries)
____
| . -> pressure point
|
-----
*/
///*******************************************************************************************************


    ///                                              Set Grid + problem

    n=50;N=1+(2*n);          //NUMBER OF PRESSURE POINTS/TOTAL GRID POINTS
    double Lx=1000;          //GRID LENGHT(x ft)
    double deltaX=Lx/N;      //DELTA(distance between points)
    double Ly=1000;          //GRID WIDTH(y ft)
    double Lz=75;            //GRID HEIGHT(z ft)
    double Ax=Ly*Lz;         //cross-sectional area normal to the x axis
    double Vb=Lx*Ly*Lz;       //bulk volume (ft³)
    ///more about constants on pages 34 and 104(Ertkein Tucay)
    double cl=0.0000035;     //compressibility of phase (in psi-¹)
    double kx=0.015;         //permeability in the direction of the x axis (in darcy)
    double por=0.18;         //porosity
    double bl=1.0;           //formation volume factor (phase l)
    double deltaT=15;        //time variation(days)
    double Mil=10;           //viscosity of phase l(in cp)
    double qlsci=-150;        //production rate of phase l at standartconditions for gridblock i(for problem 5.8->i=4)
    int target=4 ;         //thats the block with production rate
    ///more about alphac and betac in page 63(Ertkein Tucay)
    double alphal=5.614583;   //volume conversion factor(for costumary units)
    double betac=1.127;       //transmissibility conversion factor(for costumary units)
    ///transmissibility for uniformblocks Txi-1/2=Txi+1/2
    double Tlx=betac*((Ax*kx)/(Mil*bl*Lx)); //transmissibility =0 (for no-flow boundary condition) ->first and last element of the grid)
    ///backward-difference aproximation to the flow equation
    double BDA=(Vb*por*cl)/(alphal*bl*deltaT); //then,there's pressure variation (Pin+1-Pin) multiplying the equation
    ///program variables////////////////////////////////
    PetscScalar *P ; //individual pressure of a vector
    PetscScalar val;
    PetscReal r_norm;//residual norm
    Vec p,y,u,r,e;
    Mat A; //Discretization matrix
    double startt;
    int itnum=0;
    double bool_Prate=0; //0->there's no production rate 1.0-production rate detected
    int MAXGRIDS=0;
    int gridindex=0;//used in grid selection
    int n_temp=n;//used to creategrids
    PetscReal ResNorm;




    /*
    ==Vectors==
    p-pressure values to be calculated(Pstart=6000 psia)
    y-rhs
    u-exact solution
    r-residual
    e-error
    ==Matrices==
    A-Fine grid Matrix
    R-Restriction Matrix(fine->coarse )
    I-Interpolation Matrix(R transposed) (coarse->fine)

    */

    //*REMEMBER : only calculate pressure in X1,X2...Xn points!!!(size n vector)
    ierr=VecCreate(PETSC_COMM_WORLD,&p);CHKERRQ(ierr);
    ierr=VecSetSizes(p,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr=VecSetFromOptions(p);CHKERRQ(ierr);
    ierr=VecDuplicate(p,&y);CHKERRQ(ierr);
    ierr=VecDuplicate(p,&u);CHKERRQ(ierr);
    ierr=VecDuplicate(p,&r);CHKERRQ(ierr);
    ierr=VecDuplicate(p,&e);CHKERRQ(ierr);


    VecSet(p,6000);//initial pressure in psia

    ///                                              Discretize Problem


    ierr=MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
    ierr=MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,n,n);CHKERRQ(ierr);
    ierr=MatSetFromOptions(A);CHKERRQ(ierr);
    ierr=MatSetUp(A);CHKERRQ(ierr);                   //MatSetUp obrigatory!!

        value[0] =Tlx ; value[1] =-1.0*((2.0*Tlx)+BDA); value[2] =Tlx;
                     for (i=1; i<n-1; i++) {
                       col[0] = i-1; col[1] = i; col[2] = i+1;
                       ierr = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
                                                                                  }
        value[0] = Tlx; value[1]=-1.0*(Tlx+BDA);//(for no-flow boundary condition)
        i = n - 1; col[0] = n - 2; col[1] = n - 1;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
        i = 0; col[0] = 0; col[1] = 1;value[0] = -1.0*(Tlx+BDA); value[1] =Tlx;
        ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    printMatrixToFile(A,"matrix");


    /// ***************************************************************************************************

    ///set Exact solution
    /*
    Uxx=F(x)
    r=F(x)-A*v //residual(what we want)
    U=(integrated twice)F(x)
    */

    ///**********************************************************
    ///             set RHS and solve simultaneously
    ///**********************************************************

    /*solve for p
    variables : pressure(p),error(e),residual(r),rhs(y)
    constants :

    */

    //for test purposes
  //  cout<<"\nenter itnum : ";
    itnum=1;
    startt = MPI_Wtime();
    int *multigridits=new int[itnum]; //counts preconditioning


        for(int it=1;it<=itnum;it++){

            VecGetArray(p,&P);
            ///WE must update the rhs elements before each iteration...(problem varies on time)
                for(int index=0;index<n;index++){

                        if (index==target-1) bool_Prate=1.0; //here we assign target=target-1 because gridblock 4 has index=3,since the counter starts from 0
                            else  bool_Prate=0;//notice the bool variable ,if (0),then qlsci(production rate) is equal to zero

                    val=-((bool_Prate*qlsci)+(BDA*P[index]));// individual pressure value
                    ierr=VecSetValues(y,1,&index,&val,INSERT_VALUES);CHKERRQ(ierr);

                }
            VecRestoreArray(p,&P);

            KSP_linearsolver(A,y,p,its,1);
            ///Multigrid (type : V-cycle )
                  // do{
                      //  oldits=its;            //its remains static ,when the convergence condition is reached

                        ///TEST COMPARISON RICHARDSON+JACOBI x MULTIGRID

                        //Richardson Relaxation
                        //KSP_Jacobi(A,y,p,its,200); ///RICHARDSON+JACOBI
                        VecNorm(e,NORM_2,&ResNorm);
                      //  cout<<"\nERRORNORM_atMain.cppbefore Multigrid->"<<ResNorm<<endl;
                        Controller_V MG(n,r,e,p,y,A,&multigridits[it-1]);///MULTIGRID
                         VecNorm(e,NORM_2,&ResNorm);
                     //   cout<<"\nERRORNORM_atMain.cppafter Multigrid->"<<ResNorm<<endl;

                    //x }while(its!=oldits);

        }
    cout<<"\n\ndays :  "<<deltaT*itnum<<endl<<endl;
    cout << "Total Time elapsed[s]: " <<  MPI_Wtime() - startt << endl;
    cout<<"\n***MULTIGRID ITS *** "<<endl;
        for(int i=0;i<itnum;i++){
            cout<<"\nfor "<<deltaT*(i+1)<<" days :"<<multigridits[i]<<endl;
        }
        delete[] multigridits;
    cout<<endl;


    printVectorToFile(y,"rhs");
    printVectorToFile(p,"pressure");
    printVectorToFile(r,"residual");

    return 0;
}

int KSP_linearsolver(Mat A, Vec y, Vec v, PetscInt &its,int itnum){
	KSP ksp;
    PetscErrorCode ierr;
	PC preconditioner;
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
	ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);                     //DIFFERENT_NON_ZERO_PATTERN is no longer a parameter
	ierr = KSPSetType(ksp,KSPRICHARDSON);CHKERRQ(ierr);
	ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
	ierr = PCSetType(preconditioner,PCKSP);CHKERRQ(ierr);
    ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,itnum);CHKERRQ(ierr);
    ierr = KSPSolve(ksp,y,v);CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
   // cout<<"\n\ncurrent iteration : "<<its<<endl;
   // ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);                      //variable adress is now required on KSPDDestroy Funcion

        return 0;
}

int KSP_Jacobi(Mat A, Vec y, Vec v, PetscInt &its,int itnum){
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
   // cout<<"\n\ncurrent iteration : "<<its<<endl;
    //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);                      //variable adress is now required on KSPDDestroy Funcion

        return 0;
}



int printVectorToFile(Vec& v,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	VecView(v,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}

int printMatrixToFile(Mat& A,const char* filename){
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	MatView(A,viewer);
	PetscViewerDestroy(&viewer);
	return 0;
}



