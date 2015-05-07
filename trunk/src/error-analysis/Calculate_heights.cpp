#include "ErrorAnalysis.h"

void ErrorAnalysis::calculate_h_new(GeomData* pGCData, double average_error, FIELD field){
	for (int k=0; k<pGCData->getNumElements(); k++){
		// if element is singular, h_new is already known: h_new = h_mim
		if (!isSingular(k)){
			double error = getElementError(k);
			double h_old = pGCData->getElem_CDL(k);
			double h_new = h_old;
			if ( fabs(error) > 1e-8){
				h_new = h_old*(average_error/error);
			}
			set_h_new(h_new,k,field);
//			if ( (double)(h_new/h_old) < 1){
//				cout << setprecision(5) << h_old << " " << h_new << "\tratio: " << h_new/h_old << endl;
//			}
		}
	}
}

void ErrorAnalysis::identify_singular_regions(GeomData* pGCData, FIELD field){
	int counter = 0;
	const double h_min = get_h_min();
	for (int k=0; k<pGCData->getNumElements(); k++){
		double h_new = get_h_new(k,field);
		if (h_new < h_min){
			setElementAsSingular(k);
			set_h_new(h_new,k,field);			// put a limit to h_new size. It cannot be less than h_min
			counter++;
		}
	}
	cout << "Number of singular elements: " << counter << endl;
}

void ErrorAnalysis::calculate_h_ratio(GeomData* pGCData){
	for (int k=0; k<pGCData->getNumElements(); k++){
		// take the smallest h_new from both pressure and saturation fields
		double h_new = std::min( get_h_new(k,0), get_h_new(k,1) );
		double h_old = pGCData->getElem_CDL(k);
		pGCData->setElem_HR(k,h_new/h_old);
	}
}

void ErrorAnalysis::getElementsForAdaptation(double param1, double param2, GeomData* pGCData, std::list<int>& elemList){
	int k = 0;
	int dim = pGCData->getMeshDim();
	int pos = dim+1;
	int ndom = pGCData->getNumDomains();
	const int* indices = NULL;
	for (int dom=0; dom<ndom; dom++){
		for (int row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			double h_ratio = pGCData->getElem_HR(k);

			if ( (h_ratio > param2) || (h_ratio < param1) ){
				pGCData->getElement(dom,row,indices);
				this->pElmToRemove[k] = true;
				for (int i=0;i<dim+1;i++){
					elemList.push_back(indices[i+pos]);
				}
			}
			k++;
		}
	}
	cout << "Number of elements to be removed: " << elemList.size() << endl;
}

void ErrorAnalysis::getNodesForAdaptation(GeomData* pGCData, std::map<int,double>& nodeMap){

	int i,dom,row,ID;
	int numNodes = pGCData->getNumNodes();
	int k = 0;
	int dim = pGCData->getMeshDim();
	int pos1 = dim+1;
	int pos2 = 2*dim-1;
	int ndom = pGCData->getNumDomains();
	const int* indices = NULL;
	double h1, h2;

	// initialize
	for(i=0; i<numNodes; i++){
		ID = pGCData->getNodeID(i);
		nodeMap[ID] = 0;
	}

	// summation
	for (dom=0; dom<ndom; dom++){
		for (row=0; row<pGCData->getNumElemPerDomain(dom); row++){
			pGCData->getElement(dom,row,indices);

			h1 = std::min( get_h_new(k,0), get_h_new(k,1) );
			for (i=pos1; i<pos2; i++){
				h2 = nodeMap[ indices[i] ];
				nodeMap[ indices[i] ] = h1 + h2;
			}
			k++;
		}
	}

	// weighting heights
	for(i=0; i<numNodes; i++){
		nodeMap[pGCData->getNodeID(i)] /= pGCData->getNumFacesSharingVertex(i);
	}
}














