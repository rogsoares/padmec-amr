/*
 * CalculateElementsError.cpp
 *
 *  Created on: 24/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

void ErrorAnalysis_2D::calculate_ElementsError(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad pGetGradient){

	int dim = theMesh->getDim();
	double v0_coords[3], v1_coords[3], v2_coords[3];	// face's nodes coords array
	double vector_0[3], vector_1[3], vector_2[3];		// face's edge vectors
	double grad_0[3], grad_1[3], grad_2[3];				// face's node gradients

	// tag to get local node ID per domain
	int dom_counter = 0; // one domain for while (soon multiple domains)
	int row_I, row_J, row_K;
	char tag[4]; sprintf(tag,"%d",dom_counter);
	
	#ifdef __ERROR_ANALYSIS_DEBUG__
	double grad[3] = {.0,.0,.0};
	#endif

	// loop over all leave elements (without children)
	FIter fit = M_faceIter (theMesh);
	while (pEntity face = FIter_next(fit)){
		if ( !theMesh->getRefinementDepth(face) ){

			// get face's vertices
			pVertex vertices[3] = {face->get(0,0), face->get(0,1), face->get(0,2)};

			// get face's vertices coordinates
			V_coord(vertices[0],v0_coords);
			V_coord(vertices[1],v1_coords);
			V_coord(vertices[2],v2_coords);

			// make vector
			makeVector(v0_coords,v1_coords,vector_0);
			makeVector(v1_coords,v2_coords,vector_1);
			makeVector(v2_coords,v0_coords,vector_2);

			// get all three face's node gradients for 'i' field
			pSimPar->getLocalNodeIDNumbering(vertices[0],tag,row_I);
			pSimPar->getLocalNodeIDNumbering(vertices[1],tag,row_J);
			pSimPar->getLocalNodeIDNumbering(vertices[2],tag,row_K);
			pGetGradient(dom_counter,row_I,grad_0);
			pGetGradient(dom_counter,row_J,grad_1);
			pGetGradient(dom_counter,row_K,grad_2);
			
			#ifdef __ERROR_ANALYSIS_DEBUG__
			for (int i=0;i<dim;i++){
				//cout << grad_0[0] << " " << grad_0[2] << " " << grad_1[0] << " " << grad_1[2] << " " << grad_2[0] << " " << grad_2[2] << endl;
				grad[i] += grad_0[i] + grad_1[i] + grad_2[i];
			}
			#endif
			
			// get delta gradient
			double delta_grad0[3] = {grad_1[0]-grad_0[0], grad_1[1]-grad_0[1], grad_1[2]-grad_0[2]};		// edge: 0 - 1
			double delta_grad1[3] = {grad_2[0]-grad_1[0], grad_2[1]-grad_1[1], grad_2[2]-grad_1[2]};		// edge: 0 - 2
			double delta_grad2[3] = {grad_0[0]-grad_2[0], grad_0[1]-grad_2[1], grad_0[2]-grad_2[2]};		// edge: 1 - 2
			double error = (fabs(dot(delta_grad0,vector_0,2)) + fabs(dot(delta_grad1,vector_1,2)) + fabs(dot(delta_grad2,vector_2,2)))/3.0;
			this->setElementError(face, sqrt(error));
		}
	}
	FIter_delete(fit);
	#ifdef __ERROR_ANALYSIS_DEBUG__
	if (fabs(grad[0]+grad[1]+grad[2]) < 1e-8){
		throw Exception(__LINE__,__FILE__,"Gradients are NULL!");
		
	}
	///STOP();
	#endif
}
