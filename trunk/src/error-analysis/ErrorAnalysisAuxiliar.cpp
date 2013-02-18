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
