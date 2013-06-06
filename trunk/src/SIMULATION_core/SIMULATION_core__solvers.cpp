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
}

/*                       ATENCAO
 * A chamada do solver de press‹o faz usa da malha apontada por theMesh, logo esta Ž quem deve ser adaptada. Isto Ž feito por meio
 * da malha pIData->m1 e NUNCA por pIData->m2. Ap—s a adaptacao, a interpolacao Ž de pIData->m2 para pIData->m1, o preprocessamento
 * Ž sobre pIData->m1. Quando o ciclo de adaptacao se encerrar, o solver, que faz usa de theMesh, ira calcular o novo campo de pressoes
 * sobre pIData->m1, pois este ponteiro aponta para o mesmo endereco de memoria que theMesh.
 *
 * RogŽrio, 15/05/2013
 */
int SIMULATION_core::steadyState(){
	PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
	bool adapt;
	double tol1 = pSimPar->getToleranceForAllElements();
	double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();

	pPPData->setSimulationState(true);
	pIData->m1 = theMesh;
	PADMEC_mesh *pm = new PADMEC_mesh;
	int contador =0;
	do{
		contador ++;
		pElliptic_eq->solver(theMesh);
		pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
		if ( pSimPar->userRequiresAdaptation() ){
			adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),1);
			if ( adapt ){
				makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSaturation);
				//if (contador<2){
				std::list<pEntity> elementList;
				std::set<pEntity> elementSet;
				pErrorAnalysis->getRefUnrefElementsList(theMesh,elementList,elementSet);
				int i = 1;
				double Hold; 
				pMeshAdapt->rodar(pErrorAnalysis,pIData->m1);
				
				//pSimPar->printOutVTK(pIData->m1,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);//STOP();
				
				
				makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);
				//pSimPar->printOutVTK(pm,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
				//pSimPar->printOutVTK(pIData->m2,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);//STOP();
				// GAMBIARRA PARA ENCONTRAR FLAGS DAS ENTIDADES DE CONTORNO: NÓS COM POÇO DE INJEÇÃO, 
				// ARESTAS DE CONTORNO ENTRE DOMINIOS E ELEMENTOS DOS DOMINIOS
				// OBS.: SAULO, ISSO DEVE SER HERDADO JA NA ADAPÇÃO
				// return (ent->getClassification())?GEN_id( ent->getClassification() ):0;
				PADMEC_GAMBIARRA(pIData->m1);

				cout<< "Interpolador"<<endl;
				pInterpolateData(pIData);
				pSimPar->printOutVTK(pIData->m1,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);//STOP();

				cout<< "EBFV1_preprocessor"<<endl;
				EBFV1_preprocessor(pIData->m1,pGCData);
				#ifdef __ADAPTATION_DEBUG__
					cout<< "EBFV1"<<endl;
					validate_EBFV1(pGCData,pIData->m1,pSimPar->setOfDomains);	// this function must be removed to preprocessor.cpp
					if (!pSimPar->setOfDomains.size()){
						throw Exception(__LINE__,__FILE__,"Num domains NULL!\n");
					}
				#endif

				cout<< "UPDATEPOINTERS"<<endl;
				theMesh = pIData->m1;
				updatePointersData(theMesh);
				deleteMesh(pIData->m2);
				pIData->m2 = MS_newMesh(0);
				
				pIData->m2->modifyState(3,2,1);
				pIData->m2->modifyState(3,1,1);
				pIData->m2->modifyState(3,0);

				pIData->m2->modifyState(2,1);
				pIData->m2->modifyState(2,3);
				pIData->m2->modifyState(2,0);
			    
				pIData->m2->modifyState(1,3);
				pIData->m2->modifyState(1,2);
				pIData->m2->modifyState(1,0);
			    
				pIData->m2->modifyState(0,2);
				pIData->m2->modifyState(0,1);
				pIData->m2->modifyState(0,3);
				
				deleteMesh(pm);
				if (contador ==3){
				//STOP();
				}
			}
		}
	}while(adapt);
	PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
	cout<< "Loops :"<<contador <<endl;
	return 0;
}

int SIMULATION_core::transient(){
	/*
	 * Open log file monitor
	 */
//	LogFiles(OPENLG,0,0,0,0,pSimPar->getOutputPathName(),
//			pSimPar->useRestart(),
//			pSimPar->getTStepNumber(),			/* if restart=false, returns 0*/
//			pSimPar->getCPU_time());			/* if restart=false, returns .0*/
//	//int count = 0;
//	double timeStep;
//	//	double time_step_summation = .0;
//	while ( !pSimPar->finishSimulation() ){
//		double t1 = pElliptic_eq->solver();															// calculate pressure field and gradients
//		double t2 = pHyperbolic_eq->solver(timeStep);												// calculate saturation field and gradients
//		// Output data (Log)
//		LogFiles(UPDATELG,t1,t2,timeStep,pSimPar->getAccumulatedSimulationTime());
//	}
//	LogFiles(CLOSELG);
	return 0;
}
}
