#include "ErrorEstimator.h"
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>
#include <cmath>


ErrorEstimator::ErrorEstimator (pMesh theMesh){}

ErrorEstimator::~ErrorEstimator (){
	VCL.clear();
	AllNodes.clear();
	ListedFaces.clear();
}

void ErrorEstimator::LeitorDeLista (ErrorAnalysis *pErrorAnalysis, pMesh theMesh){

	std::list<pEntity> elementList;
	std::set<pEntity> elementSet;
	pErrorAnalysis->getRefUnrefElementsList(theMesh,elementList,elementSet);

	for (std::list<pEntity>::iterator iter = elementList.begin(); iter!=elementList.end(); iter++){
		ListedFaces.insert(*iter);

	}
}




void ErrorEstimator::VertexList (pMesh theMesh){ // Pega os elementos do set de faces a remover e levanta todos os vertex

	set<pFace>::iterator itfaces;
	for (itfaces=ListedFaces.begin() ; itfaces != ListedFaces.end(); itfaces++ ){
		pVertex vt;
		for(int y=0;y<=2;y++){
			vt=F_vertex(*itfaces,y);
			AllNodes.insert(vt);
		}
	}

	cout << "A malha possui " << M_numVertices(theMesh) << " nodes ao total" << endl;
	cout << "AllNodes.size(): " << AllNodes.size() << endl;
}












