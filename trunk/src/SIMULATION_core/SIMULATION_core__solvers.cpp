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
		return 0;
	}

	/*                       ATENCAO
	 * A chamada do solver de presso faz usa da malha apontada por theMesh, logo esta  quem deve ser adaptada. Isto  feito por meio
	 * da malha pIData->m1 e NUNCA por pIData->m2. Aps a adaptacao, a interpolacao  de pIData->m2 para pIData->m1, o preprocessamento
	 *  sobre pIData->m1. Quando o ciclo de adaptacao se encerrar, o solver, que faz usa de theMesh, ira calcular o novo campo de pressoes
	 * sobre pIData->m1, pois este ponteiro aponta para o mesmo endereco de memoria que theMesh.
	 *
	 * Rogrio, 15/05/2013
	 */
	int SIMULATION_core::steadyState(){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
#ifndef NOADAPTATION
		double tol1 = pSimPar->getToleranceForAllElements();
		double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();
		int numFields = 1;
		pIData->m1 = theMesh;
		PADMEC_mesh *pm = new PADMEC_mesh;
#endif
		bool adapt = false;
		pPPData->setSimulationState(false);
		int count = 0;

		do{
			pElliptic_eq->solver(theMesh);
			pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
#ifndef NOADAPTATION
			adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),numFields);
			pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
			if ( adapt ){
				makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSaturation_Old);
				pMeshAdapt->rodar(pErrorAnalysis,pIData->m1);
				deleteMesh(pIData->m1);pIData->m1 = 0;
				deleteMesh(theMesh); theMesh = 0;
				pIData->m1 = MS_newMesh(0);
				readmesh(pIData->m1,"Final_01.msh");
				theMesh = pIData->m1;
				makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);
				PADMEC_GAMBIARRA(pIData->m1);
				cout<< "Interpolador"<<endl;
				interpolation(pIData,pSimPar->getInterpolationMethod());
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
				EBFV1_preprocessor(pIData->m1,pGCData);
	#ifdef __ADAPTATION_DEBUG__
				validate_EBFV1(pGCData,pIData->m1,pSimPar->setOfDomains);
				if (!pSimPar->setOfDomains.size()){
					throw Exception(__LINE__,__FILE__,"Num domains NULL!\n");
				}
	#endif
				updatePointersData(theMesh);
				// After interpolation of pressure field calculate (not interpolate) the gradient pressure for the new mesh
				deleteMesh(pm);
				deleteMesh(pIData->m2); pIData->m2 = 0;
				pIData->m2 = MS_newMesh(0);
			}
			count++;
#endif
		}while (adapt);
		PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
		cout<< "Loops :"<< count <<endl;
		return 0;
	}
	
	int SIMULATION_core::transient(){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
#ifndef NOADAPTATION
		double adaptStep = 0;				// performs mesh adaptation every 20 timeSteps
		bool adapt;
		double tol1 = pSimPar->getToleranceForAllElements();
		double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();
		int numFields = 2;
		pIData->m1 = theMesh;
		PADMEC_mesh *pm = new PADMEC_mesh;
#endif
		LogFiles(OPENLG,0,0,0,0,pSimPar->getOutputPathName(),pSimPar->useRestart(),pSimPar->getTStepNumber(),pSimPar->getCPU_time());
				double timeStep;
				double time_step_summation = .0;
		while ( !pSimPar->finishSimulation() ){
			pElliptic_eq->solver(theMesh);
			pHyperbolic_eq->solver(theMesh,timeStep);
			pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
#ifndef NOADAPTATION
			if ( pSimPar->userRequiresAdaptation() ){
				adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),numFields);
				if (adapt){
					makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSaturation_Old);
					pMeshAdapt->rodar(pErrorAnalysis,pIData->m1);
					deleteMesh(pIData->m1);pIData->m1 = 0;
					deleteMesh(theMesh); theMesh = 0;
					pIData->m1 = MS_newMesh(0);
					readmesh(pIData->m1,"Final_01.msh");
					theMesh = pIData->m1;
					makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);
					PADMEC_GAMBIARRA(pIData->m1);
					cout<< "Interpolador"<<endl;
					interpolation(pIData,pSimPar->getInterpolationMethod());
					EBFV1_preprocessor(pIData->m1,pGCData);

					#ifdef __ADAPTATION_DEBUG__
						validate_EBFV1(pGCData,pIData->m1,pSimPar->setOfDomains);
						if (!pSimPar->setOfDomains.size()){
							throw Exception(__LINE__,__FILE__,"Num domains NULL!\n");
						}
					#endif

					updatePointersData(theMesh);
					deleteMesh(pm);
					deleteMesh(pIData->m2); pIData->m2 = 0;
					pIData->m2 = MS_newMesh(0);
				}
			}
#endif
		}
		return 0;
	}
}
