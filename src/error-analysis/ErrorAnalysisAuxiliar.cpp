/*
 * ErrorAnalysis_auxiliar.cpp
 *
 *  Created on: 29/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

void getMaxMinData(pMesh theMesh, pEntity elem, ErrorAnalysis* pEA, int &maxDepth, int &maxRefLevel, int &minRefLevel);

double ErrorAnalysis::calculate_ErrorSum(pMesh theMesh){
	double error_sum = .0;
	if (theMesh->getDim()==2){
		pEntity face;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				error_sum -= pow(this->getElementError(face),2);
			}
		}
		FIter_delete(fit);
	}
	else{
		pEntity tetra;
		RIter rit = M_regionIter (theMesh);
		while ( (tetra = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(tetra)){
				error_sum += this->getElementError(tetra);
			}
		}
		RIter_delete(rit);
	}
	return error_sum;
}

void ErrorAnalysis::initializeParameters(pMesh theMesh){
	if (theMesh->getDim()==2){
		pEntity face;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			//if(!theMesh->getRefinementDepth(face)){
			setElementError(face,.0);
			setLevelOfRefinement(face,0);	// force element level to assume a real value and not a fictitious one like -1000.
			//}
		}
		FIter_delete(fit);
	}
	else{
		pEntity tetra;
		RIter rit = M_regionIter (theMesh);
		while ( (tetra = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(tetra)){
				setElementError(tetra,.0);
				setLevelOfRefinement(tetra,0);	// force element level to assume a real value and not a fictitious one like -1000.
			}
		}
		RIter_delete(rit);
	}
}

void ErrorAnalysis::countElements(pMesh theMesh, bool excludingSingularities=false){
	int numElem = 0;
	int numElem_exSing = 0;
	if (theMesh->getDim()==2){
		pEntity face;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				if (excludingSingularities && !this->isSingular(face)){
					numElem_exSing++;
				}
				else{
					numElem++;
				}
			}
		}
		FIter_delete(fit);
	}
	else{
		pEntity tetra;
		RIter rit = M_regionIter (theMesh);
		while ( (tetra = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(tetra)){
				if (excludingSingularities && !this->isSingular(tetra)){
					numElem_exSing++;
				}
				else{
					numElem++;
				}
			}
			RIter_delete(rit);
		}
	}
	if (excludingSingularities){
		setNumElements_excludingSingularities(numElem_exSing);
	}
	else{
		setNumElements(numElem);
	}
}

void ErrorAnalysis::calculate_MaxMinErrorData(pMesh theMesh){
	int maxDepth = 0;
	int maxRefLevel = -10;
	int minRefLevel = 10;
	pEntity elem;
	if (theMesh->getDim()==2){
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			getMaxMinData(theMesh,elem,this,maxDepth,maxRefLevel,minRefLevel);
		}
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter (theMesh);
		while ( (elem = RIter_next(rit)) ){
			getMaxMinData(theMesh,elem,this,maxDepth,maxRefLevel,minRefLevel);
		}
		RIter_delete(rit);
	}
	this->setMaxDepth(maxDepth);
	this->setMaxRefinementFlag(maxRefLevel);
	this->setMinRefinementFlag(minRefLevel);
}

void ErrorAnalysis::monitoring(pMesh theMesh){
	static int i = 0;
	static bool openfile = true;
	if (openfile){
		fid.open("Error_Analysis_Monitor.txt");
		
		fid << "                            ERROR ANALISYS MONITOR\n"
		"================================================================================\n"
		"Legend:\n"
		"i -------- Number of steps until adaptation is not necessary anymore\n"
		"ge ------- GlobalError\n"
		"sgn ------ Smoothe Gradient Norm\n"
		"sgns  ---- Smoothe Gradient Norm Excluding Singularity\n"
		"me ------- Mean Error\n"
		"mes ------ Mean Error Excluding Singulariry\n"
		"ne ------- Number of Elements\n"
		"================================================================================\n\n"
		"--------------------------------------------------------------------------------\n"
		"i       ge           sgn          sgns         me           mes          ne\n"
		"--------------------------------------------------------------------------------\n";
		//		fid << "GlobalError   SmootheGradNorm     SmootheGradNorm_singularity    MeanError  MeanError_singulariry  numElements\n\n";
		openfile = false;
	}
	fid << setprecision(5) << scientific;
	fid << ++i << "   " << getGlobalError() << "  " << getSmoothedGradNorm() << "  " <<  getSmoothedGradNorm_singular() <<
	"  " << getAverageError() << "  " <<  getAverageError_excludingSingularities() << "     " << this->getNumElements() << endl;
	
}

void getMaxMinData(pMesh theMesh, pEntity elem, ErrorAnalysis* pEA, int &maxDepth, int &maxRefLevel, int &minRefLevel){
	maxDepth = std::max(maxDepth,theMesh->getRefinementDepth(elem));
	maxRefLevel = std::max(maxRefLevel,pEA->getLevelOfRefinement(elem));
	minRefLevel = std::min(minRefLevel,pEA->getLevelOfRefinement(elem));
}

bool ErrorAnalysis::checkMaximumNumberOfSubdivision(pMesh theMesh, const int &maxNumberOfSubdivision){
	pEntity elem;
	int level = -1000;		/// initialize variable to get maximum number of refinement
	if (theMesh->getDim()==2){
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(elem)){
				level = std::max(theMesh->getRefinementLevel(elem),level);
			}
		}
		FIter_delete(fit);
	}
	else{
		RIter rit = M_regionIter (theMesh);
		while ( (elem = RIter_next(rit)) ){
			if(!theMesh->getRefinementDepth(elem)){
				level = std::max(theMesh->getRefinementLevel(elem),level);
			}
		}
		RIter_delete(rit);
	}
	return (level==maxNumberOfSubdivision);
}

void ErrorAnalysis::getRefUnrefElementsList(pMesh theMesh, std::list<pEntity> &elementList, std::set<pEntity>& nodesBGMesh){
	// avoid messing things
	if (!elementList.size()){
		elementList.clear();
	}
	
	pEntity elem;
	std::list<double> heightList;
	double param1 = _pSimPar->Remeshing_param1();
	double param2 = _pSimPar->Remeshing_param2();
	printf("param1: %.5f\tparam2: %.5f\n",param1,param2);
	double avgError, h_old, h_new, error;
	if (!theMesh){
		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
	}
	int dim = theMesh->getDim();
	if (dim==2){
		avgError = getAverageError();
		FIter fit = M_faceIter (theMesh);
		while ( (elem = FIter_next(fit)) ){
			error = getElementError(elem);
			h_old = getElement_CDL(elem);
			h_new = h_old*(avgError/error);
			if ( h_new < param1*h_old || h_new > param2*h_old ){
				elementList.push_back(elem);
				heightList.push_back(h_new);
			}
		}
		FIter_delete(fit);
	}
	else if (dim==3){
		throw Exception(__LINE__,__FILE__,"Under construction!");
		//		RIter rit = M_regionIter (theMesh);
		//		while ( (elem = RIter_next(rit)) ){
		//			if(!theMesh->getRefinementDepth(elem)  && this->getLevelOfRefinement(elem)){
		//				elementList.push_back(elem);
		//			}
		//		}
		//		RIter_delete(rit);
	}
	else{
		throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
	}
	
	// calculate average element height per node. Each element node have the average new height of all elements choose to be (un)refined
	int i;
	double height;
	pEntity v;
	std::list<pEntity>::iterator iter1 = elementList.begin();
	std::list<double>::iterator iter2 = heightList.begin();
	for ( ;iter1 != elementList.end(); iter1++, iter2++){
		for (i=0; i<dim; i++){
			v = (pEntity)(*iter1)->get(0,i);
			nodesBGMesh.insert(v);
			height = .0;
			EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
			height += *iter2;
			EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),height);
		}
	}
	
	std::set<pEntity>::iterator iter3 = nodesBGMesh.begin();
	for ( ;iter3 != nodesBGMesh.end(); iter3++){
		v = *iter3;
		height = .0;
		EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
		//printf("height: %.5f\t",height);
		int n = V_numFaces(v);
		//printf("V_numFace: %d\t",n);
		height /= n;
		//printf("height_avg: %.5f\n",height);
		EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),height);
	}
	cout << "nodesBGMesh.size() = " << nodesBGMesh.size() << endl;
}


//void ErrorAnalysis::getRefUnrefElementsList(pMesh theMesh, std::list<pEntity> &elementList, std::set<pEntity>& nodesBGMesh){
//	// avoid messing things
//	if (elementList.size()){
//		elementList.clear();
//	}
//	if (nodesBGMesh.size()) {
//		nodesBGMesh.clear();
//	}
//	
//	pEntity elem;
//	std::list<double> heightList;
//	double param1 = _pSimPar->Remeshing_param1();
//	double param2 = _pSimPar->Remeshing_param2();
//	
//	
//	theMesh->modifyState(3,2,1);
//	theMesh->modifyState(3,1,1);
//	theMesh->modifyState(3,0);
//	
//	theMesh->modifyState(2,1);
//	theMesh->modifyState(2,3);
//	theMesh->modifyState(2,0);
//	
//	theMesh->modifyState(1,3);
//	theMesh->modifyState(1,2);
//	theMesh->modifyState(1,0);
//	
//	theMesh->modifyState(0,2);
//	theMesh->modifyState(0,1);
//	theMesh->modifyState(0,3);
//	
//	double avgError, h_old, h_new, error;
//	int dim = theMesh->getDim();
//	if (dim==2){
//		avgError = getAverageError();
//		FIter fit = M_faceIter (theMesh);
//		while ( (elem = FIter_next(fit)) ){
//		pEdge e1 = F_edge(elem = FIter_next(fit), 0);
//		pEdge e2 = F_edge(elem, 1);
//		pEdge e3 = F_edge(elem, 2);
//		while(( E_numFaces(e1) <= 1 ) || ( E_numFaces(e2) <= 1 ) || ( E_numFaces(e3) <= 1 )){
//			for(int k = 0; k < rand()%200; k++) {
//			    cout << EN_id(elem);
//			    elem = FIter_next(fit);
//			    e1 = F_edge(elem,0);
//			    e2 = F_edge(elem,1);
//			    e3 = F_edge(elem,2);
//			}
//		}
//		error = getElementError(elem);
//		h_old = getElement_CDL(elem);
//		h_new = h_old*(avgError/error);
//		cout << "hold: "<<h_old <<" error: "<< error << endl;
//		
//		if ( h_new < param1*h_old || h_new > param2*h_old ){
//			elementList.push_back(elem);
//			heightList.push_back(h_new);
//		}
//		
//		elementList.push_back(elem);
//		heightList.push_back(0.5f);
//		
//		cout << "elemento base inserido" << endl;
//		
//		
//		set<pFace> faces;
//		//for{
//			//pEdge 
//			e1 = F_edge(elem,0);
//			e2 = F_edge(elem,1);
//			e3 = F_edge(elem,2);
//			pFace f1 ;
//			//faces adjacentes
//			
//			f1 = (E_face(e1, 1) != elem) ? E_face(e1, 1) : E_face(e1, 0);
//			cout << "F1: " << EN_id(f1) << endl;
//			elementList.push_back(f1);
//			heightList.push_back(0.5f);
//			
//			pFace f2;
//			if(E_numFaces(e2) > 1) {
//				f2 = (E_face(e2, 1) != elem) ? E_face(e2, 1) : E_face(e2, 0);
//				cout << "F2: " << EN_id(f2) << endl;
//			    elementList.push_back(f2);
//				heightList.push_back(0.5f);
//			}
//			pFace f3;
//			if(E_numFaces(e3) > 1) {
//				f3 = (E_face(e3, 1) != elem) ? E_face(e3, 1) : E_face(e3, 0);
//			    cout << "F3: " << EN_id(f3) << endl;
//				elementList.push_back(f3);
//				heightList.push_back(0.5f);
//			}
//		//}
//		
//		
//		
//		//}
//		FIter_delete(fit);
//	}
//	else if (dim==3){
//		throw Exception(__LINE__,__FILE__,"Under construction!");
//		//		RIter rit = M_regionIter (theMesh);
//		//		while ( (elem = RIter_next(rit)) ){
//		//			if(!theMesh->getRefinementDepth(elem)  && this->getLevelOfRefinement(elem)){
//		//				elementList.push_back(elem);
//		//			}
//		//		}
//		//		RIter_delete(rit);
//	}
//	else{
//		throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
//	}
//	// calculate average element height per node. Each element node have the average new height of all elements choose to be (un)refined
//	int i;
//	double height;
//	double h2;
//	pEntity v;
//	std::list<pEntity>::iterator iter1 = elementList.begin();
//	std::list<double>::iterator iter2 = heightList.begin();
//	for ( ;(iter1 != elementList.end()) && (iter2 != heightList.end()) ; iter1++, iter2++){
//		for (i=0; i<dim; i++){
//		    cout << "Chegou aqui1" << endl;
//			v = (pEntity)(*iter1)->get(0,i);
//			cout << "Chegou aqui2 v: " << v << endl;
//			nodesBGMesh.insert(v);
//			cout << "Chegou aqui3" << endl;
//			height = .0;
//			EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
//			cout << "Chegou aqui4 height:" << height << "iter2: " << *iter2 << endl;
//			height += *iter2;
//			h2 = height;
//			cout << "Chegou aqui4.5 ID:" << EN_id(v) << endl;
//			EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),h2);
//			cout << "Chegou aqui5 height:" << height << endl;
//		}
//	}
//	cout << "Chegou aqui" << endl;
//	std::set<pEntity>::iterator iter3 = nodesBGMesh.begin();
//	for ( ;iter3 != nodesBGMesh.end(); iter3++){
//		v = *iter3;
//		height = .0;
//		EN_getDataDbl(v,MD_lookupMeshDataId( "elem_height" ),&height);
//		printf("height: %.5f\t",height);
//		int n = V_numFaces(v);
//		printf("V_numFace: %d\t",n);
//		height /= n;
//		printf("height_avg: %.5f\n",height);
//		EN_attachDataDbl(v,MD_lookupMeshDataId( "elem_height" ),height);
//	}
//	cout << "nodesBGMesh.size() = " << nodesBGMesh.size() << endl;
//	STOP();
//}
