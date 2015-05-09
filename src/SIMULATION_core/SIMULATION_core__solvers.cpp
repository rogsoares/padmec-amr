/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"
//#include <sstream>

namespace PRS{

	int SIMULATION_core::solver(){
		double timeStep;
		if (simFlag==STEADY_STATE){
			bool adapt = false;	// ............ If adaptation is not required, get out while loop.
			int count = 0;
			do{
				pElliptic_eq->solver(theMesh);
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);
				#ifndef NOADAPTATION
					adapt = adaptation();
				#endif
				count++;
			}while (adapt);
			if (adapt){
				cout<< "Loops :"<< count <<endl;
			}
		}
		else if (simFlag==TRANSIENT){
			while ( !pSimPar->finishSimulation() ){
				pElliptic_eq->solver(theMesh);
				pHyperbolic_eq->solver(theMesh,timeStep);

				std::list<int> elemList;
				std::map<int,double> nodeMap;
				//calculate_ErrorAnalysis(pErrorAnalysis,pSimPar,pGCData,PhysicPropData::getGradient,elemList,nodeMap);
				//pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);

				//STOP();

				//#ifndef NOADAPTATION
				//	adaptation();
				//#endif
			}
		}
		return 0;
	}
}
