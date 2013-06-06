#include "exportVTK.h"

void exportSolutionToVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, string filename){
	// open file
	ofstream fid;
	fid.open(filename.c_str());

	// check if it's OK, otherwise terminate program
	char msg[512]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n"
			"\tCheck if directory was typed correctly.\n",filename.c_str());

	if ( !fid.is_open() ){
		throw Exception(__LINE__,__FILE__,msg);
	}

	PRS::PhysicPropData *pPPData  = (PRS::PhysicPropData*)pData1;
	ErrorAnalysis *pErrorAnalysis  = (ErrorAnalysis*)pData2;
	PRS::SimulatorParameters *pSimPar = (PRS::SimulatorParameters*)pData3;
	// print data to file
	fid << "# vtk DataFile Version 2.0\n";
	fid << "Two phases flow simulation\n";
	fid << "ASCII\nDATASET UNSTRUCTURED_GRID\n";
	fid << "POINTS " << M_numVertices(theMesh) << " float\n";

	int dim = theMesh->getDim();
	pErrorAnalysis->countElements(theMesh,false);
	int numElements = pErrorAnalysis->getNumElements();
	
	if (!numElements){
	  throw Exception(__LINE__,__FILE__,"Number of elements NULL!");
	}

	printVerticesCoordenates(fid,theMesh);
	printElementConnectivities(fid,theMesh,dim,numElements);
	printCellTypeList(fid,dim,numElements);
	printPressure(fid,theMesh,pPPData);
	printSaturation(fid,theMesh,pPPData);

	if ( pSimPar->userRequiresAdaptation() ){
//		print_Swgrad(fid,theMesh,pSimPar,pPPData);
//		print_Sw_GradientNorm2(fid,theMesh,pErrorAnalysis,pSimPar,pPPData);
		printDegreeOfRefinement(fid,theMesh,pErrorAnalysis);
//		printElementError(fid,theMesh,pErrorAnalysis);
//		printCharac_Lenth(fid,theMesh,pErrorAnalysis);
//		print_pw_GradientNorm(fid,theMesh,pErrorAnalysis,pSimPar,pPPData);
//		print_Sw_GradientNorm(fid,theMesh,pErrorAnalysis,pSimPar,pPPData);
	}

	fid.close();
}

// print vertex coordenates and transfer to each one all computed values
// to be printed later
// = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = =
void printVerticesCoordenates(ofstream &fid, pMesh theMesh){
	pEntity e;
	int count = 0;
	VIter vit = M_vertexIter(theMesh);
	while( (e = VIter_next(vit)) ) {
		double coord[3] = {.0, .0, .0};
		V_coord(e,coord);
		for (int i=0; i<3; i++) fid << coord[i] << " ";
		fid << endl;
		EN_attachDataInt(e,MD_lookupMeshDataId("mLN"),count++);
	}
	VIter_delete(vit);
}

// print elements connectivities
void printElementConnectivities(ofstream &fid, pMesh theMesh, int dim, int numElements){
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
}

void printElementConnectivities(ofstream &fid, pEntity elem, int dim){
	int mappedLNodes;
	fid << dim+1 << " ";
	for(int i=0; i<dim+1; i++){
		EN_getDataInt((pVertex)elem->get(0,i),MD_lookupMeshDataId("mLN"),&mappedLNodes);
		fid << mappedLNodes << " ";
	}
	fid << endl;
}

void printCellTypeList(ofstream &fid, int dim, int numElements){
	fid << "\nCELL_TYPES " << numElements << endl;
	int type = (dim==2)?5:10;
	for(int i=0; i<numElements; i++)
		fid << type << endl;
}

void printPressure(ofstream &fid, pMesh theMesh, PRS::PhysicPropData *pPPData){
	// print pressure
	// = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = =
	fid << "\nPOINT_DATA "<< M_numVertices(theMesh) << endl;
	fid << "SCALARS Pressure float 1\n";
	fid << "LOOKUP_TABLE default\n";
	pEntity node;
	double val;
	VIter vit = M_vertexIter(theMesh);
	while( (node = VIter_next(vit)) ){
		val = pPPData->getPressure(node);
		fid << val << endl;
	}
	VIter_delete(vit);
}

void printSaturation(ofstream &fid, pMesh theMesh, PRS::PhysicPropData *pPPData){
	// print saturation
	// = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = =
	fid << "SCALARS Saturation float 1\n";
	fid << "LOOKUP_TABLE default\n";
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while( (node = VIter_next(vit)) ){
		double val = .0;
		val = pPPData->getSaturation(node);
		fid << val << endl;
	}
	VIter_delete(vit);
}

void print_Sw_GradientNorm2(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData){
	// print saturation
	// = = = = = = = = = = = = = = = =  = = = = = = = = = = = = = = = = = =
	double grad[3];
	int dom_counter = 0; // one domain for while (soon multiple domains)
	int row, dim = theMesh->getDim();
	char tag[4]; sprintf(tag,"%d",dom_counter);

	fid << "SCALARS Satgrad float 1\n";
	fid << "LOOKUP_TABLE default\n";
	pEntity node;
	VIter vit = M_vertexIter(theMesh);
	while( (node = VIter_next(vit)) ){
		pSimPar->getLocalNodeIDNumbering(node,tag,row);
		pPPData->get_Sw_Grad(dom_counter,row,grad);
		double norm_0 = .0;
		for (int j=0;j<dim;j++){
			norm_0 += grad[j]*grad[j];
		}
		fid << sqrt(norm_0) << endl;
	}
	VIter_delete(vit);
}

void printDegreeOfRefinement(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis){
	fid << "\nCELL_DATA " << pErrorAnalysis->getNumElements() << endl;
	fid << "SCALARS LevelOfRefinement int 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	FIter fit = M_faceIter(theMesh);
	while(pFace face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){
			fid << pErrorAnalysis->getLevelOfRefinement(face) <<endl;
		}
	}
	FIter_delete(fit);
}

void printElementError(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis){
	fid << "SCALARS Element_Error float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	FIter fit = M_faceIter(theMesh);
	while(pFace face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){
			fid << pErrorAnalysis->getElementError(face) << endl;
		}
	}
	FIter_delete(fit);
}

void printCharac_Lenth(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis){
	fid << "SCALARS Charac_Lenth float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;
	FIter fit = M_faceIter(theMesh);
	while(pFace face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){
			fid << pErrorAnalysis->getElement_CDL(face) << endl;
		}
	}
	FIter_delete(fit);
}

void print_pw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData){
	fid << "SCALARS pwGrad_Norm float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;

	double grad[3];
	int dom_counter = 0; // one domain for while (soon multiple domains)
	int row, dim = theMesh->getDim();
	char tag[4]; sprintf(tag,"%d",dom_counter);

	FIter fit = M_faceIter(theMesh);
	while(pFace face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){
			double norm_1 = .0;
			for (int i=0;i<3;i++){
				pSimPar->getLocalNodeIDNumbering(face->get(0,i),tag,row);
				pPPData->get_pw_Grad(dom_counter,row,grad);
				double norm_0 = .0;
				for (int j=0;j<dim;j++){
					norm_0 += grad[j]*grad[j];
				}
				norm_1 += sqrt(norm_0);
			}
			norm_1 /= 3.0;
			fid << norm_1 << endl;
		}
	}
	FIter_delete(fit);
}

void print_Sw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData){
	fid << "SCALARS SwGrad_Norm float 1 " << endl;
	fid << "LOOKUP_TABLE default " << endl;

	double grad[3];
	int dom_counter = 0; // one domain for while (soon multiple domains)
	int row, dim = theMesh->getDim();
	char tag[4]; sprintf(tag,"%d",dom_counter);

	FIter fit = M_faceIter(theMesh);
	while(pFace face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){
			double norm_1 = .0;
			for (int i=0;i<3;i++){
				pSimPar->getLocalNodeIDNumbering(face->get(0,i),tag,row);
				pPPData->get_Sw_Grad(dom_counter,row,grad);
				double norm_0 = .0;
				for (int j=0;j<dim;j++){
					norm_0 += grad[j]*grad[j];
				}
				norm_1 += sqrt(norm_0);
			}
			norm_1 /= 3.0;
			fid << norm_1 << endl;
		}
	}
	FIter_delete(fit);
}

void print_Swgrad(ofstream &fid, pMesh theMesh, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData){
	fid << "VECTORS p_grad float\n";
	pEntity node;
	double grad[3];
	int dom_counter = 0; // one domain for while (soon multiple domains)
	int row;
	//int dim = theMesh->getDim();
	char tag[4]; sprintf(tag,"%d",dom_counter);

	VIter vit = M_vertexIter(theMesh);
	while( (node = VIter_next(vit)) ){
		pSimPar->getLocalNodeIDNumbering(node,tag,row);
		pPPData->get_pw_Grad(dom_counter,row,grad);
		fid << grad[0] << " " << grad[1] << " .0\n";
	}
	VIter_delete(vit);
}
