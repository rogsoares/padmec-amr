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
void ErrorAnalysis_2D::calculate_SmoothedGradientNorm(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad pGetGradient){
	double grad[3] = {.0,.0,.0};
	double SGN = .0;
	int row, dom_counter = 0;
	char tag[4]; sprintf(tag,"%d",dom_counter);
	pEntity face;
	double gradDiff = .0;
	FIter fit = M_faceIter (theMesh);
	while ( (face = FIter_next(fit)) ){
		if(!theMesh->getRefinementDepth(face)){
			for (int i=0;i<3;i++){
				mVertex* vertex = (mVertex*)face->get(0,i);		// get face vertex
				pSimPar->getLocalNodeIDNumbering(vertex,tag,row);
				pGetGradient(dom_counter,row,grad);
				//gradDiff += fabs(grad[0] - grad[1]);			// sum gradient difference (x,y)
				gradDiff += dot(grad,grad,2);
			}
			gradDiff /= 3.0;
			SGN += gradDiff;
		}
	}
	FIter_delete(fit);
	this->setSmoothedGradNorm(sqrt(SGN));
}

void ErrorAnalysis_2D::calculate_SmoothedGradientNorm_Singularity(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad pGetGradient){
	double grad[3];
	double SGN = .0;
	int row, dom_counter = 0;
	char tag[4]; sprintf(tag,"%d",dom_counter);
	pEntity face;
	FIter fit = M_faceIter (theMesh);
	while ( (face = FIter_next(fit)) ){
		if(!theMesh->getRefinementDepth(face) && !this->isSingular(face)){
			double gradDiff = .0;
			for (int i=0;i<3;i++){
				mVertex* vertex = (mVertex*)face->get(0,i);		// get face vertex
				pSimPar->getLocalNodeIDNumbering(vertex,tag,row);
				pGetGradient(dom_counter,row,grad);
				//gradDiff += fabs(grad[0] - grad[1]);			// sum gradient difference (x,y)
				gradDiff += dot(grad,grad,2);
			}
			gradDiff /= 3.0;
			SGN += gradDiff;
		}
	}
	FIter_delete(fit);
	this->setSmoothedGradNorm_Singularity( sqrt(SGN) );
}
