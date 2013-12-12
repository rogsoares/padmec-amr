#ifndef _ELLIPTIC_EQUATION_H_
#define _ELLIPTIC_EQUATION_H_

#include "Restart.h"
#include "exportVTK.h"

/**
 * Base class for elliptic solvers.
 */
namespace PRS           // PRS: Petroleum Reservoir Simulator
{
	class Elliptic_equation{
	public:

		Elliptic_equation(){}
		virtual ~Elliptic_equation(){}
		virtual double solver(pMesh)=0;
		virtual double pressureGradient(pMesh)=0;

	protected:

		PetscErrorCode ierr;
		Vec output;

		// solve system of equations: Ax=y
		double KSP_solver(Mat A, Mat pcMatrix, Vec y, Vec x, SimulatorParameters *pSimPar, PetscTruth guessNonZero, KSPType ksptype, PCType pctype,PetscInt &its){
			double startt = MPI_Wtime();
			KSP ksp;
			PC preconditioner;
			ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRQ(ierr);
			ierr = KSPSetOperators(ksp,A,pcMatrix,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
			ierr = KSPSetType(ksp,ksptype);CHKERRQ(ierr);
			ierr = KSPGetPC(ksp,&preconditioner);CHKERRQ(ierr);
			ierr = PCSetType(preconditioner,pctype);CHKERRQ(ierr);
			ierr = KSPSetInitialGuessNonzero(ksp,guessNonZero);CHKERRQ(ierr);
			ierr = KSPGMRESSetRestart(ksp,pSimPar->getKrylov_restart()); CHKERRQ(ierr);
			ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
			ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
			ierr = KSPSolve(ksp,y,x);CHKERRQ(ierr);
			ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
			ierr = KSPDestroy(ksp); CHKERRQ(ierr);
			double endt = MPI_Wtime();
			return endt-startt;
		}
	};
}
#endif
