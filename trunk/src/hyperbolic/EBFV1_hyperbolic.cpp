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
//		pStruct = new PointerStruct;
		pPPData = ppd;
		pSimPar = sp;
//		pStruct->pErrorAnalysis = pEA;
//		pStruct->pSimPar = sp;
		pMData = md;
		_cumulativeOil = .0;
	}

	EBFV1_hyperbolic::~EBFV1_hyperbolic(){
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
		pPPData->resetNonvisc();
		alpha_max = .0;
		int dim = theMesh->getDim();
		// calculate saturation gradient if adaptation or high order approximation were required
		if (  pSimPar->userRequiresAdaptation() || pSimPar->useHOApproximation()){
			calculateSaturationGradient(theMesh);
		}

		// initialize time step with a very high number
		timeStep = 1.0e+10;
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calculateVelocityField(dom,dim);
			calculateIntegralAdvectiveTerm(dom,timeStep);
		}
		pSimPar->correctTimeStep(timeStep);					// correct time-step value to print out the desired simulation moment
		pSimPar->setCumulativeSimulationTime(timeStep); 	// AccSimTime = AccSimTime + timeStep
		calculateExplicitAdvanceInTime(timeStep);	// Calculate saturation field: Sw(n+1)

	#ifdef _SEEKFORBUGS_
		if (!P_pid()) std::cout << "<_SEEKFORBUGS_>  timeStep = : " << timeStep << endl;
		if (timeStep==.0) throw Exception(__LINE__,__FILE__,"Time step NULL!");
	#endif

		timestep_counter++;

		// oil production output
		if (pSimPar->rankHasProductionWell() && pSimPar->timeToPrintVTK()){
			pOPManager->printOilProduction(timeStep,pSimPar->getCumulativeSimulationTime(),pSimPar->getSimTime(),getRecoveredOil(),getCumulativeOil(),timestep_counter);
			timestep_counter = 0;
		}

		if (!P_pid()) std::cerr << "done.\n\n";
		return hyp_time;
	}
}
