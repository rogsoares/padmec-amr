#ifndef _EBFV1_ELLIPTIC_EQUATION_H_
#define _EBFV1_ELLIPTIC_EQUATION_H_

#include "Elliptic_equation.h"

namespace PRS           // PRS: Petroleum Reservoir Simulator
{
	/**
	 * For EBFV1 elliptic formulation, a matrix-free procedure is applied to solve pressure field. A pointer-function is used to tell 
	 * Petsc how to make a matrix-vector product and all matrices and vector are passed to this function as a structure.
	 */
	struct Data_struct{
		Mat *E;
		Mat G;
		Mat *F;
		Vec RHS;
		int ndom;
		int nrows;
		int *rows;

		Vec z;
		int F_nrows;
		int F_ncols;
	};

	struct Vectors{
		Vec u_old;
		Vec v_tmp;
		Vec resvec;
		Vec rhs;
		bool key;
	};


	int MatMultUser(Mat mat, Vec u, Vec y);

	class EBFV1_elliptic : public Elliptic_equation{
	public:

		EBFV1_elliptic();
		EBFV1_elliptic(pMesh, PhysicPropData *, SimulatorParameters *, GeomData *, MeshData *);
		~EBFV1_elliptic();
		double solver(pMesh);
		double pressureGradient(pMesh theMesh);

	private:
		bool DF_key;
		dblarray Lij;
		Vectors* pVec;
		Mat *EF_multiDom;

		// After fill the matrix, it must be assembled (PETSC)
		int assemblyMatrix(Mat);
		
		// Fill matrices E,G and F for a specific edge and domain. Theses function are called inside a loop of domains
		int divergence_E(Mat, pEntity, const int&, int, double*);
		int divergence_G(Mat, pEntity, const int&, int, double*);
		int gradient_F_edges(Mat, pEntity, const int&, int, double*);
		int gradient_F_bdry(pMesh, Mat, const int&);
		int F_bdryFaces(pMesh, Mat, const int&);
		int F_bdryEdges(pMesh, Mat, const int&);

		 // Solves (EF + G)u = q as showed below:
		 // 			G * u^k+1 = q - EF * u^k		 
		double solveIteratively();

		// Associate to mesh nodes new pressure values computed
		double updatePressure(pMesh theMesh);
		
		// Compute pressure gradient and associate them to mesh nodes for all domains. 
		// Note: nodes between two or more domains can store a vector for each one
		int pressureGradient(int, int);
		int resetPressureGradient(pMesh theMesh, int, char*);
		int updatePressureGradient(int, int);

		// PETSc matrix for matrix-free procedure
		Mat matrix;

		// Pointer to be used inside matrix-free procedure
		Data_struct *matvec_struct;
		double setMatrixFreeOperation(pMesh);

		/* Group set of calls to compute and assembly matrices E,F and G*/
		double assembly_EFG_RHS(pMesh);

		/* MatMultUser is called several time by Petsc during KSP solver until convergence be reached. */
		static int MatMultUser(Mat mat, Vec u, Vec y){
			void *ctx;
			PetscErrorCode ierr = MatShellGetContext(mat,&ctx);  CHKERRQ(ierr);
			Data_struct* mats = (Data_struct*)(ctx);

			// the desired multiplication => A*u = y, where A = EF+G
			// step 1: y = G*u;
			ierr = MatMult(mats->G,u,y); CHKERRQ(ierr);

			// step 2:  y = y + (E*F*u)_dom1 + (E*F*u)_dom2 + ... + (E*F*u)_domN
			int m,n;
			for (int i=0; i<mats->ndom; i++){
				// z = [F]*u
				ierr = VecZeroEntries(mats->z); CHKERRQ(ierr);
				ierr = MatMult(mats->F[i],u,mats->z); CHKERRQ(ierr);
				// y = y + ([E]*z)_dom_i
				ierr = MatMultAdd(mats->E[i],mats->z,y,y); CHKERRQ(ierr);
			}
			return ierr;
		}

		/*set well contribution to right hand side*/
		int wellsContributionToRHS(pMesh theMesh, Vec&);

		/* Matrices E,F and G are assembled for all nodes.
		 * After that, sub matrices are extracted to form the system of equations.
		 * It is not seem to be the wisest way to assembly the final system of equation
		 * but EBFV1 formulation is a little bit complex to make something efficient.*/
		int set_SOE(pMesh theMesh, Mat, Mat&, bool, Vec&, bool, bool);

		/*Free memory from matrices related to Data_struct*/
		double freeMemory();


		PhysicPropData *pPPData;
		SimulatorParameters *pSimPar;
		GeomData *pGCData;
		MeshData *pMData;
		pMesh theMesh;
	};
}
#endif
