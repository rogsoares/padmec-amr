/*
 * MEBFV_Elliptic_main.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "MEBFV_elliptic.h"

namespace PRS{
	MEBFV_elliptic::MEBFV_elliptic(){
	}

	MEBFV_elliptic::MEBFV_elliptic(pMesh mesh, PhysicPropData *pdata, SimulatorParameters *psimpar, GeomData *pgdata, MeshData *pmdata){
		theMesh = mesh;
		pPPData = pdata;
		pSimPar = psimpar;
		pGCData = pgdata;
		pMData = pmdata;
		this->initialize = false;
	}

	MEBFV_elliptic::~MEBFV_elliptic(){
	}

	// solves system of equation for pressure field
	double MEBFV_elliptic::solver(pMesh theMesh){
	#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tIN\n";
	#endif

		Initialize();
		Assembly_A();
		Assembly_b();
		Solve();

	#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tOUT\n";
	#endif
		return 0;
	}

	double MEBFV_elliptic::Solve(){
		double startt = MPI_Wtime();
		KSP ksp;
		PC preconditioner;
		KSPCreate(PETSC_COMM_WORLD,&ksp);
		KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);
		KSPSetType(ksp,KSPCG);
		KSPGetPC(ksp,&preconditioner);
		PCSetType(preconditioner,PCJACOBI);
		KSPSetInitialGuessNonzero(ksp,PETSC_FALSE);
		KSPSetFromOptions(ksp);
		KSPSolve(ksp,b,x);
		//KSPGetIterationNumber(ksp,&its);
		KSPDestroy(ksp);
		double endt = MPI_Wtime();
		return endt-startt;
	}
}
