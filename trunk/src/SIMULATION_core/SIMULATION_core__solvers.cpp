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
		LogFiles(0,0,0,0,0,0,OPENLG,pSimPar->getOutputPathName(),pSimPar->useRestart(),pSimPar->getStepOutputFile(),pSimPar->getCumulativeSimulationTime(),pSimPar->getCPU_time());
		double timeStep, t1, t2, t3, t4, CPU_cum;
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
			double assemblyT,hsolverT,psolverT,gradT;
			int KSPiter;
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(theMesh);
				pElliptic_eq->getCPUtime(assemblyT,psolverT,gradT,KSPiter);
				hsolverT = pHyperbolic_eq->solver(theMesh,timeStep);
				#ifndef NOADAPTATION
					adaptation();
				#endif
				LogFiles(timeStep,assemblyT,psolverT,gradT,KSPiter,hsolverT,UPDATELG,pSimPar->getOutputPathName(),pSimPar->useRestart(),
						pSimPar->getStepOutputFile(),pSimPar->getCumulativeSimulationTime(),pSimPar->getCPU_time());
			}
		}
		return 0;
	}
}
