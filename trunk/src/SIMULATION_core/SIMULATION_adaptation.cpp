/*
 * SIMULATION_adaptation.cpp
 *
 *  Created on: Jan 31, 2014
 *      Author: rogerio
 */

#include "SIMULATION_core.h"
//#ifndef NOADAPTATION

namespace PRS{
	bool SIMULATION_core::adaptation(){
		bool adapt = false;
		pSimPar->set_adapt_occur(false);

		if ( pSimPar->userRequiresAdaptation() ){
			// performs an error estimation on saturation and/or pressure solution and verifies if tolerances are obeyed.
			// If error is greater than tolerance then mesh will be adapted to improve solution quality

			std::list<int> elemList;
			std::map<int,double> nodeMap;
			bool adapt = calculate_ErrorAnalysis(pErrorAnalysis,pSimPar,pGCData,PhysicPropData::getGradient,elemList,nodeMap);
			pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);


			if (adapt){
				// let other simulation code parts aware of the occurrence of the adaptation process
				pSimPar->set_adapt_occur(true);

				// retrieve cumulative simulation time
				pSimPar->retrieveSimulationTime();

				pIData->m1 = theMesh;				// auxiliary pointer
				pIData->m2 = MS_newMesh(0);			// initialize auxiliary pointer
				PADMEC_mesh *pm = new PADMEC_mesh;	// temporary pointer

				// preserves current mesh and create a new one adapted
				makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSw_old);

				// mesh adaptation (adapted mesh saved in file. Bad, bad, bad!)
				//pMeshAdapt->run(pIData->m1,elemToAdapt,backgroundNodes);

				// delete any trace of FMDB data structure
				deleteMesh(theMesh); theMesh = 0;

				//create a new FMDB data structure with data in file
				pIData->m1 = MS_newMesh(0);
				readmesh(pIData->m1,"simulation-parameters/FS-2D-homogeneo/pp-data-files/fs2dhmg__3.msh");
				theMesh = pIData->m1;

				// take mesh (coords and connectivities) from temporary matrix to m2 (avoid conflicts with FMDB when deleting m1 and theMesh)
				makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);

				// before interpolate data from the old to the new mesh, vectors where pressure and saturation are stored
				// must be allocated for the new mesh (it has just to be created)
				pPPData->allocateTemporaryData(M_numVertices(pIData->m1));

				// interpolate data from m2 to m1
				//interpolation(pIData,pSimPar->getInterpolationMethod());

				// transfer Sw and pressure from tmp to main struct
				pPPData->transferTmpData();

				//pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,pGCData,exportSolutionToVTK);

				// waste old data and get ready for new data
				pGCData->deallocatePointers(0);

				// mesh pre-processor
				EBFV1_preprocessor(pIData->m1,pGCData);

				pGCData->initilize(theMesh);

				// objects (pSimPar, pMData, pGCData, pPPData) must store new data
				updatePointersData(theMesh);

				// data transfer from FMDB to matrices
				pGCData->dataTransfer(theMesh);

				// temporary mesh not necessary any more (goodbye!)
				deleteMesh(pm);
				delete pm; pm = 0;
				deleteMesh(pIData->m2); delete pIData->m2; pIData->m2 = 0;
			}
		}
		return adapt;
	}
}

//#endif
