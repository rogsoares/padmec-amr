/*
 * CalculateSmoothGradientNorm.cpp
 *
 *  Created on: 24/04/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"


/*
 * Calculate smooth gradient Norm (SGN):
 * 	SNG = sum((|gradSw_I|+|gradSw_J|+|gradSw_K|)/3.0)^2
 * 	where sum is for each triangle element with nodes I, J, K.
 */
void ErrorAnalysis_2D::calculate_SmoothedGradientNorm(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad pGetGradient, FIELD field){
	double gradDiff, grad[3], SGN;
	int row,i;
	pEntity face;
	int dom_counter = 0;
	gradDiff = .0;
	SIter_const dom=pSimPar->setDomain_begin();
	for (; dom!=pSimPar->setDomain_end();dom++){		// loop over domains to calculate error element
		char tag[4]; sprintf(tag,"%d",dom_counter);
		SGN = .0;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face)){
				for (i=0; i<3; i++){
					mVertex* vertex = (mVertex*)face->get(0,i);		// get face vertex
					pSimPar->getLocalNodeIDNumbering(vertex,tag,row);
					pGetGradient((pVertex)vertex,dom_counter,row,grad);
					gradDiff += dot(grad,grad,2);
				}
				gradDiff /= 3.0;
				SGN += gradDiff;
			}
		}
		FIter_delete(fit);
		// Allow only one loop over elements to calculate error based on saturation field. Saturation gradient is not per domain like pressure.
		if (field == SATURATION){
			break;
		}
		dom_counter++;
	}
	setSmoothedGradNorm(sqrt(SGN));
}

void ErrorAnalysis_2D::calculate_SmoothedGradientNorm_Singularity(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad pGetGradient, FIELD field){
	double grad[3], gradDiff, SGN;
	int i, row, dom_counter = 0;
	pEntity face;
	gradDiff = .0;
	SIter_const dom=pSimPar->setDomain_begin();
	for (; dom!=pSimPar->setDomain_end();dom++){		// loop over domains to calculate error element
		char tag[4]; sprintf(tag,"%d",dom_counter);
		SGN = .0;
		FIter fit = M_faceIter (theMesh);
		while ( (face = FIter_next(fit)) ){
			if(!theMesh->getRefinementDepth(face) && !this->isSingular(face)){
				for (i=0;i<3;i++){
					mVertex* vertex = (mVertex*)face->get(0,i);		// get face vertex
					pSimPar->getLocalNodeIDNumbering(vertex,tag,row);
					pGetGradient((pVertex)vertex,dom_counter,row,grad);
					gradDiff += dot(grad,grad,2);
				}
				gradDiff /= 3.0;
				SGN += gradDiff;
			}
		}
		FIter_delete(fit);
		// Allow only one loop over elements to calculate error based on saturation field. Saturation gradient is not per domain like pressure.
		if (field == SATURATION){
			break;
		}
		dom_counter++;
	}
	setSmoothedGradNorm_Singularity(sqrt(SGN));
}
