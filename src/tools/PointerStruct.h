/*
 * PointerStruct.h
 *
 *  Created on: 11/07/2009
 *      Author: rogerio
 */

#ifndef POINTERSTRUC_H_
#define POINTERSTRUC_H_

#include "PhysicPropData.h"
#include "ErrorAnalysis.h"

/*
 * PointerStruct holds all pointers necessary for the simulator.
 * Each pointer has a specific task providing parameters from input file like
 * numereic, physical and pre-processor data. Almost all simulator class need
 * them.
 */

namespace PRS{

	struct PointerStruct{

		PhysicPropData *pPPData;
		SimulatorParameters *pSimPar;
		GeomData *pGCData;
		pMesh theMesh;
		ErrorAnalysis *pErrorAnalysis;
	};

}

#endif /* POINTERSTRUC_H_ */
