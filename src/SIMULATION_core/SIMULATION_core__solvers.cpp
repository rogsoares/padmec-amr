/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"

namespace PRS{

int SIMULATION_core::solver(){
	switch ( simFlag ){
	case STEADY_STATE:
		this->steadyState();
		break;
	case TRANSIENT:
		this->transient();
		break;
	default:
		throw Exception(__LINE__,__FILE__, Exception::INIT_ERROR );
	}
}

int SIMULATION_core::steadyState(){
	PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
	bool adapt;
	double tol1 = pSimPar->getToleranceForAllElements();
	double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();

	pPPData->setSimulationState(true);
	do{
		pElliptic_eq->solver();
		if ( pSimPar->userRequiresAdaptation() ){
			adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),1);
			if ( adapt ){
				pMeshAdapt->run(pErrorAnalysis,theMesh);
				pInterpolateData(pIData);
				EBFV1_preprocessor(theMesh,pGCData);
				#ifdef __ADAPTATION_DEBUG__
				validate_EBFV1(pGCData,theMesh,pSimPar->setOfDomains);	// this function must be removed to preprocessor.cpp
				#endif
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
				updatePointersData();
			}
		}
	}while(adapt);
	PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
	return 0;
}

int SIMULATION_core::transient(){
	/*
	 * Open log file monitor
	 */
	LogFiles(OPENLG,0,0,0,0,pSimPar->getOutputPathName(),
			pSimPar->useRestart(),
			pSimPar->getTStepNumber(),			/* if restart=false, returns 0*/
			pSimPar->getCPU_time());			/* if restart=false, returns .0*/
	//int count = 0;
	double timeStep;
//	double time_step_summation = .0;
	while ( !pSimPar->finishSimulation() ){
		double t1 = pElliptic_eq->solver();															// calculate pressure field and gradients
		double t2 = pHyperbolic_eq->solver(timeStep);												// calculate saturation field and gradients
		// Output data (Log)
		LogFiles(UPDATELG,t1,t2,timeStep,pSimPar->getAccumulatedSimulationTime());
	}
	LogFiles(CLOSELG);
	return 0;
}
}
