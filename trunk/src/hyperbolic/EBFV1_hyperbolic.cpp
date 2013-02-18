/*
 * EBFV1-hyperbolic.cpp
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

	EBFV1_hyperbolic::EBFV1_hyperbolic(){
	}

	EBFV1_hyperbolic::EBFV1_hyperbolic(pMesh mesh, PhysicPropData *ppd,
			SimulatorParameters *sp, GeomData *gcd,
			MeshData *md, OilProductionManagement *popm, ErrorAnalysis *pEA){

		pGCData = gcd;
		pOPManager = popm;
		/*
		 * Experimental:
		 */
		pStruct = new PointerStruct;
		pStruct->pPPData = ppd;
		//pGCData = gcd;
		pStruct->pSimPar = sp;
		pStruct->theMesh = mesh;
		pStruct->pErrorAnalysis = pEA;
		pHOApproximation = new HighOrderApproximation(pStruct);
		pMData = md;
		//rowToImport = 0;
	}

	EBFV1_hyperbolic::~EBFV1_hyperbolic(){
		delete pStruct;
		delete pHOApproximation;
//		if (rowToImport){
//			delete[] rowToImport; rowToImport = 0;
//			MatDestroy(joinNodes);
//			MatDestroy(updateValues);
//		}
	}

	double EBFV1_hyperbolic::solver(double &timeStep){
	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "Hyperbolic solver...\n";
	#endif
		double hyp_time = .0;

		// set to 0 all fluxes across control volumes from previous time-step
		resetNodalNonviscTerms();

		/*
		 * For each domain, calculate:
		 * 1 - velocity field
		 * 2 - saturation gradient and slope limiters (if required by user)
		 * 3 - fluxes across surfaces from control volumes
		 */

		// loop over domains
		// calculate saturation gradient if adaptation or high order approximation were required
		if (  pStruct->pSimPar->userRequiresAdaptation() || pStruct->pSimPar->useHOApproximation())
			hyp_time += calculateSaturationGradient();

		int dom_counter = 0;
		for (SIter_const dom=pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end();dom++){
			hyp_time += calculateVelocityField(*dom,dom_counter);

			/*
			 * Some physical/numeric properties must be calculated before get high
			 * order approximations for saturation.
			 */
			if ( pStruct->pSimPar->useHOApproximation() ){
				NodeSlopeLimiter* pNodeSL = pHOApproximation->getNodeSL_Ptr();
				hyp_time += pNodeSL->defineSlopeLimiters();
			}
			hyp_time += calculateIntegralAdvectiveTerm(*dom);
			dom_counter++;
		}
		//pMData->unifyScalarsOnMeshNodes(PhysicPropData::getNonViscTerm,PhysicPropData::setNonViscTerm,pGCData,0);

		// take minimum time-step calculated by domain
		 timeStep = getTimeStep();

	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "<_SEEKFORBUGS_>  timeStep = : " << timeStep << endl;
	#endif

		// correct time-step value to print out the desired simulation moment
		 pStruct->pSimPar->correctTimeStep(timeStep);

	#ifdef _SEEKFORBUGS_
		if (timeStep==.0) throw Exception(__LINE__,__FILE__,"Time step NULL!");
	#endif
		// AccSimTime = AccSimTime + timeStep
		 pStruct->pSimPar->setAccumulatedSimulationTime(timeStep);

	#ifdef _SEEKFORBUGS_
		if (timeStep==.0) throw Exception(__LINE__,__FILE__,"Time step NULL!");
	#endif

		// Calculate saturation field: Sw(n+1)
		 hyp_time += calculateExplicitAdvanceInTime(timeStep);

		// Output data (VTK) VAZAMANETO NA SAIDA DO VTK
		pStruct->pSimPar->printOutVTK(pStruct->theMesh,pStruct->pPPData,pStruct->pErrorAnalysis,pStruct->pSimPar,exportSolutionToVTK);

		/*
		 * The following lines below will be condensed to a function member call
		 * and it will belong to EBFV1_hyperbolic
		 */
		if (pStruct->pSimPar->rankHasProductionWell()){
			pOPManager->printOilProduction(timeStep,
					pStruct->pSimPar->getAccumulatedSimulationTime(),
					pStruct->pSimPar->getSimTime(),
					getRecoveredOilValue());
		}

		if (!P_pid()) std::cout << "done.\n\n";
		return hyp_time;
	}

	// it also works for parallel simulation
	double EBFV1_hyperbolic::getTimeStep(){
		double dt = *min_element(timeStepByDomain.begin(),timeStepByDomain.end());
		timeStepByDomain.clear();
		return dt;
	}
}
