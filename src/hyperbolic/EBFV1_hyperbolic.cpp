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

	EBFV1_hyperbolic::EBFV1_hyperbolic(pMesh mesh, PhysicPropData *ppd,SimulatorParameters *sp, GeomData *gcd,MeshData *md, OilProductionManagement *popm, ErrorAnalysis *pEA){

		pGCData = gcd;
		pOPManager = popm;
		pStruct = new PointerStruct;
		pStruct->pPPData = ppd;
		pStruct->pSimPar = sp;
		pStruct->pErrorAnalysis = pEA;
		pHOApproximation = new HighOrderApproximation(pStruct);
		pMData = md;
		_cumulativeOil = .0;
	}

	EBFV1_hyperbolic::~EBFV1_hyperbolic(){
		delete pStruct;
		delete pHOApproximation;
	}

	double EBFV1_hyperbolic::solver(pMesh theMesh, double &timeStep){
	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "Hyperbolic solver...\n";
		if (!M_numEdges(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of edges: 0!\n");
		}
	#endif
		double hyp_time = .0;
		static int timestep_counter = 0;			// counts number of time steps every new VTK
		resetNodalNonviscTerms(theMesh);			// set to 0 all fluxes across control volumes from previous time-step

		/*
		 * For each domain, calculate:
		 * 1 - velocity field
		 * 2 - saturation gradient and slope limiters (if required by user)
		 * 3 - fluxes across surfaces from control volumes
		 */

		// calculate saturation gradient if adaptation or high order approximation were required
		if (  pStruct->pSimPar->userRequiresAdaptation() || pStruct->pSimPar->useHOApproximation()){
			hyp_time += calculateSaturationGradient(theMesh);
		}

		// initialize time step with a very high number
		timeStep = 1.0e10;
		int dom_counter = 0;
		for (SIter_const dom=pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end();dom++){
			hyp_time += calculateVelocityField(theMesh,*dom,dom_counter);
			if ( pStruct->pSimPar->useHOApproximation() ){
				NodeSlopeLimiter* pNodeSL = pHOApproximation->getNodeSL_Ptr();
				hyp_time += pNodeSL->defineSlopeLimiters();
			}
			hyp_time += calculateIntegralAdvectiveTerm(theMesh,*dom,dom_counter,timeStep);
			dom_counter++;
		}

	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "<_SEEKFORBUGS_>  timeStep = : " << timeStep << endl;
	#endif

		// correct time-step value to print out the desired simulation moment
		 pStruct->pSimPar->correctTimeStep(timeStep);

	#ifdef _SEEKFORBUGS_
		if (timeStep==.0) throw Exception(__LINE__,__FILE__,"Time step NULL!");
	#endif
		// AccSimTime = AccSimTime + timeStep
		 pStruct->pSimPar->setCumulativeSimulationTime(timeStep);

	#ifdef _SEEKFORBUGS_
		if (timeStep==.0) throw Exception(__LINE__,__FILE__,"Time step NULL!");
	#endif

		// Calculate saturation field: Sw(n+1)
		hyp_time += calculateExplicitAdvanceInTime(theMesh,timeStep);

		timestep_counter++;

		// oil production output
		if (pStruct->pSimPar->rankHasProductionWell() && pStruct->pSimPar->timeToPrintVTK()){
			pOPManager->printOilProduction(timeStep,
					                       pStruct->pSimPar->getCumulativeSimulationTime(),
					                       pStruct->pSimPar->getSimTime(),
					                       getRecoveredOil(),
					                       getCumulativeOil(),
					                       timestep_counter);
			timestep_counter = 0;
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
