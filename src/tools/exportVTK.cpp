#include "exportVTK.h"
#include "CPU_Profiling.h"

void exportSolutionToVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, void *pData4, string filename){
	CPU_Profile::Start();
	// open file
	ofstream fid;
	open_file(fid,filename,__LINE__,__FILE__);

	PRS::PhysicPropData *pPPData  = (PRS::PhysicPropData*)pData1;
	ErrorAnalysis *pErrorAnalysis  = (ErrorAnalysis*)pData2;
	//PRS::SimulatorParameters *pSimPar = (PRS::SimulatorParameters*)pData3;
	PRS::GeomData *pGCData = (PRS::GeomData*)pData4;

	// print data to file
	fid << "# vtk DataFile Version 2.0\n";
	fid << "Two phases flow simulation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << M_numVertices(theMesh) << " float\n";

	int dim = theMesh->getDim();
	int numElements = (dim==2)?M_numFaces(theMesh):M_numRegions(theMesh);
	
	throw_exception(!numElements,"Number of elements NULL!",__LINE__,__FILE__);

	printVerticesCoordenates(fid,theMesh);
	printElementConnectivities(fid,theMesh,dim,numElements);
	printCellTypeList(fid,dim,numElements);

	// start print nodal values
	// LIST HERE ALL NODAL FIELDS
	fid << "\nPOINT_DATA "<< M_numVertices(theMesh) << endl;
	printPressure(fid,theMesh,pPPData);
	printSaturation(fid,theMesh,pPPData,pGCData);
	printPressureGradient(fid,pGCData,pPPData);

	// LIST HERE ALL ELEMENT FIELDS
	fid << "\nCELL_DATA "<< pGCData->getNumElements() << endl;
	printElementError(fid,pGCData,pErrorAnalysis);
	print_h_ratio(fid,pGCData,pErrorAnalysis);
	print_singular_regions(fid,pGCData,pErrorAnalysis);
	print_elements_to_remove(fid,pGCData,pEA);
	
	fid.close();
	CPU_Profile::End("VTK");
}

// print vertex coordenates and transfer to each one all computed values
// to be printed later
// = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = =
void printVerticesCoordenates(ofstream &fid, pMesh theMesh){
	//CPU_Profile::Start();
	pEntity e;
	int count = 0;
	VIter vit = M_vertexIter(theMesh);
	while( (e = VIter_next(vit)) ) {
		double coord[3] = {.0, .0, .0};
		V_coord(e,coord);
		for (int i=0; i<3; i++){
			fid << coord[i] << " ";
		}
		fid << endl;
		EN_attachDataInt(e,MD_lookupMeshDataId("mLN"),count++);
	}
	VIter_delete(vit);
	//CPU_Profile::End("VTK__coord");
}

// print elements connectivities
void printElementConnectivities(ofstream &fid, pMesh theMesh, int dim, int numElements){
	////CPU_Profile::Start();
	fid << "\nCELLS " << numElements << " " << (dim+2)*numElements << endl;
	pEntity elem;
	if (dim==2){
		FIter fit = M_faceIter(theMesh);
		while( (elem = FIter_next(fit)) )
			if ( !theMesh->getRefinementDepth(elem) )
				printElementConnectivities(fid,elem,dim);
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter(theMesh);
		while ( (elem = RIter_next(rit)) )
			if ( !theMesh->getRefinementDepth(elem) )
				printElementConnectivities(fid,elem,dim);
		RIter_delete(rit);
	}
	////CPU_Profile::End("VTK__connect01");
}

void printElementConnectivities(ofstream &fid, pEntity elem, int dim){
	////CPU_Profile::Start();
	int mappedLNodes;
	fid << dim+1 << " ";
	for(int i=0; i<dim+1; i++){
		EN_getDataInt((pVertex)elem->get(0,i),MD_lookupMeshDataId("mLN"),&mappedLNodes);
		fid << mappedLNodes << " ";
	}
	fid << endl;
	////CPU_Profile::End("VTK__connect02");
}

void printCellTypeList(ofstream &fid, int dim, int numElements){
	////CPU_Profile::Start();
	fid << "\nCELL_TYPES " << numElements << endl;
	int type = (dim==2)?5:10;
	for(int i=0; i<numElements; i++){
		fid << type << endl;
	}
	////CPU_Profile::End("VTK__typelist");
}

void printPressure(ofstream &fid, pMesh theMesh, PRS::PhysicPropData *pPPData){
	//CPU_Profile::Start();
	fid << "SCALARS Pressure float 1\n";
	fid << "LOOKUP_TABLE default\n";
	pEntity node;
	double p;
	int idx=0;
	VIter vit = M_vertexIter(theMesh);
	while( (node = VIter_next(vit)) ){
		pPPData->getPressure(idx,p);
		fid << p << endl;
		idx++;
	}
	VIter_delete(vit);
	//CPU_Profile::End("VTK__pressure");
}

void printSaturation(ofstream &fid, pMesh theMesh, PRS::PhysicPropData* pPPData, PRS::GeomData* pGCData){
	////CPU_Profile::Start();
	fid << "SCALARS Saturation float 1\n";
	fid << "LOOKUP_TABLE default\n";
	double Sw;
	int nnodes = M_numVertices(theMesh);
	for(int i=0; i<nnodes; i++){
		pPPData->getSaturation(i,Sw);
		fid << setprecision(8) << fixed << Sw << endl;
	}
	////CPU_Profile::End("VTK__saturation");
}

void printElementError(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	throw_exception(!pEA,"ErrorAnalysis* pEA = NULL!",__LINE__,__FILE__);

	fid << "SCALARS Element_Error float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->getElementError(k++) << endl;
		}
	}
}

void printPressureGradient(ofstream& fid, PRS::GeomData* pGCData, PRS::PhysicPropData *pPPData){
	fid << "VECTORS p_grad float\n";
	double* p_grad = NULL;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int node=0; node<pGCData->getNumNodesPerDomain(dom); node++){
			pPPData->get_pw_Grad(dom,node,p_grad);
			fid << p_grad[0] << " " << p_grad[1] << " " << p_grad[2] << endl;
		}
	}
}

void print_h_ratio(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	throw_exception(!pEA,"ErrorAnalysis* pEA = NULL!",__LINE__,__FILE__);

	fid << "SCALARS h_ratio float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pGCData->getElem_HR(k++) << endl;
		}
	}
}

void print_singular_regions(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	throw_exception(!pEA,"ErrorAnalysis* pEA = NULL!",__LINE__,__FILE__);

	fid << "SCALARS singular float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->isSingular(k++) << endl;
		}
	}
}

void print_elements_to_remove(ofstream &fid, GeomData* pGCData, ErrorAnalysis *pEA){
	throw_exception(!pEA,"ErrorAnalysis* pEA = NULL!",__LINE__,__FILE__);

	fid << "SCALARS singular float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	int k = 0;
	for (int dom=0; dom<pGCData->getNumDomains(); dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			fid << pEA->isToRemove(k++) << endl;
		}
	}
}
