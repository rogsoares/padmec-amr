/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"
#include <sstream>

namespace PRS{

	int SIMULATION_core::solver(){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
		LogFiles(OPENLG,0,0,0,0,pSimPar->getOutputPathName(),pSimPar->useRestart(),pSimPar->getTStepNumber(),pSimPar->getCPU_time());
		double timeStep, t1, t2, t3, CPU_cum;
		if (simFlag==STEADY_STATE){
			bool adapt;
			int count = 0;
			do{
				pElliptic_eq->solver(theMesh);
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);
				#ifndef NOADAPTATION
					adapt = adaptation();
				#endif
				count++;
			}while (adapt);
			cout<< "Loops :"<< count <<endl;
		}
		else if (simFlag==TRANSIENT){
			CPU_cum = .0;
			while ( !pSimPar->finishSimulation() ){
				t1 = MPI_Wtime();
				pElliptic_eq->solver(theMesh);
				pHyperbolic_eq->solver(theMesh,timeStep);
				#ifndef NOADAPTATION
					adaptation();
				#endif
				t2 = MPI_Wtime();
				CPU_cum += t2 - t1;
				cout << setprecision(2) << fixed << "CPU time (step/cumulated)[s] " << t2 - t1 << "\t" << CPU_cum << endl;

			}
		}
		PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
		return 0;
	}
}
