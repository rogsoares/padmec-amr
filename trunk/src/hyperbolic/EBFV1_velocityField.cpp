/*
 * velocityField.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

// computes total velocity on edges for domain dom
double EBFV1_hyperbolic::calculateVelocityField(pMesh theMesh, int dom, int dom_counter){
	cout << "calculateVelocityField\n";

	double startt = MPI_Wtime();
	int i,j;
	int dim = pGCData->getMeshDim();			// mesh dimension
	dblarray edIJ(dim,.0),edIJ2(dim,.0);//, v(dim,.0);
	double pw_grad_I[3], pw_grad_J[3], pw_grad_IJ[3], vel[3];
	pEntity edge;
	double sign = 1.0;
	int row_I, row_J;
	char tag[4]; sprintf(tag,"%d",dom_counter);

	// check if pressure gradient exist for all mesh nodes
#ifdef _SEEKFORBUGS_
//	pEntity node;
//	VIter vit = M_vertexIter(theMesh);
//	while ( (node = VIter_next(vit)) ){
//		pStruct->pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
//		pStruct->pPPData->get_pw_Grad(node,dom_counter,row_I,pw_grad_I);
//		if ( fabs(pw_grad_I[0])<1e-8 && fabs(pw_grad_I[1])<1e-8 ){
//			char msg[256]; sprintf(msg,"Gradient has all its componets null for node %d\n",EN_id(node));
//			throw Exception(__LINE__,__FILE__,msg);
//		}
//	}
//	VIter_delete(vit);
#endif

	int row = 0;
	// loop over all edges from domain 'dom'
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		sign = 1.0;
		if ( pGCData->edgeBelongToDomain(edge,dom) ){
			/// get nodes I and J
			pEntity I = (pVertex)edge->get(0,0);
			pEntity J = (pVertex)edge->get(0,1);

			// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
			if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
				std::swap(I,J);
				sign = -sign;
			}

			pGCData->getEdgeVec_Unitary(edge,edIJ);
			for (i=0; i<dim; i++){
				edIJ[i] *= sign;
			}
			for (i=0; i<dim; i++){
				edIJ[i] = -edIJ[i];
			}

			pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
			pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
			pStruct->pPPData->get_pw_Grad(I,dom_counter,row_I,pw_grad_I);
			pStruct->pPPData->get_pw_Grad(J,dom_counter,row_J,pw_grad_J);
			for (i=0; i<dim; i++){
				pw_grad_IJ[i] = .5*(pw_grad_I[i] + pw_grad_J[i]);
			}

			// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
			// what must be calculated:
			//				VIJN = -PermIJ*MobIJ*(gradIJ - (gradIJ*edgVec)*edgVec)
			//				VIJPDF = -PermIJ*MobIJ*((p_J - p_I)/length)*edgVec
			// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

			// velocity normal component
			double VIJN[3] = {.0,.0,.0};
			double a = .0;
			for (i=0; i<dim; i++){
				a += edIJ[i]*pw_grad_IJ[i];
			}
			for (i=0; i<dim; i++){
				VIJN[i] = pw_grad_IJ[i] - a*edIJ[i];
			}

			// velocity parallel component - finite difference approximation
			a = (pStruct->pPPData->getPressure(J) - pStruct->pPPData->getPressure(I))/pGCData->getEdgeLength(edge);

			double VIJFD[3] = {edIJ[0]*a, edIJ[1]*a, edIJ[2]*a};

			// total velocity: vIJ
			double aux[3] = {VIJFD[0]+VIJN[0], VIJFD[1]+VIJN[1], VIJFD[2]+VIJN[2]};

			// make v_old = v_new (for t=0, v_new = 0)
//			pStruct->pPPData->getVelocity_new(edge,dom,v);
//			pStruct->pPPData->setVelocity_old(edge,dom,v);

			pStruct->pPPData->getVelocity_new(dom_counter,row,vel);
			pStruct->pPPData->setVelocity_old(dom_counter,row,vel);

			// update v_new
			//dblarray vel(dim,.0);

			// get absolute permeability tensor: dim x dim
			const double *K = pStruct->pSimPar->getPermeability(dom);

			// get mobility for nodes I and J
			const double MobI = pStruct->pPPData->getTotalMobility(I);
			const double MobJ = pStruct->pPPData->getTotalMobility(J);
			const double MobIJ = .5*(MobI + MobJ);

			int pos = 0;
			for (i=0; i<dim; i++){
				vel[i] = 0;
			}
			for (i=0; i<dim; i++){
				for (j=0; j<dim; j++){
					vel[i] += K[dim*i+j]*aux[j];
				}
			}

			for (i=0; i<dim; i++){
				vel[i] *= -MobIJ;
			}

			// total velocity: VTIJ = VIJFD + VIJN
			//pStruct->pPPData->setVelocity_new(edge,dom,vel);
			pStruct->pPPData->setVelocity_new(dom_counter,row,vel);
			row++;
		}
	}
	EIter_delete(eit);
	double endt = MPI_Wtime();
	return endt-startt;
}
}
