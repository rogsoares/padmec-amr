/*
 * velocityField.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{

// computes total velocity on edges for domain dom
double EBFV1_hyperbolic::calculateVelocityField(int dom, int dom_counter){

#ifdef _SEEKFORBUGS_
	bool check = false;
#endif

	double startt = MPI_Wtime();
	int i,j;
	int dim = pGCData->getMeshDim();			// mesh dimension
	dblarray edIJ(dim,.0),edIJ2(dim,.0), v(dim,.0);
	//dblarray pw_grad_I(3), pw_grad_J(3), pw_grad_IJ(3);
	double pw_grad_I[3], pw_grad_J[3], pw_grad_IJ[3];
	pEntity edge;
	double sign = 1.0;

//	ofstream fid;
//	fid.open("SimPar.txt");

	int row_I, row_J;
	char tag[4]; sprintf(tag,"%d",dom_counter);
	// loop over all edges from domain 'dom'
	EIter eit = M_edgeIter(pStruct->theMesh);
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
			for (i=0; i<dim; i++) edIJ[i] *= sign;
			for (i=0; i<dim; i++) edIJ[i] = -edIJ[i];
			double length = pGCData->getEdgeLength(edge);

//			pStruct->pPPData->get_pw_Grad(I,dom,pw_grad_I);
//			pStruct->pPPData->get_pw_Grad(J,dom,pw_grad_J);

			pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
			pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
			pStruct->pPPData->get_pw_Grad(dom_counter,row_I,pw_grad_I);
			pStruct->pPPData->get_pw_Grad(dom_counter,row_J,pw_grad_J);
			for (i=0; i<dim; i++) pw_grad_IJ[i] = .5*(pw_grad_I[i] + pw_grad_J[i]);

			// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
			// what must be calculated:
			//				VIJN = -PermIJ*MobIJ*(gradIJ - (gradIJ*edgVec)*edgVec)
			//				VIJPDF = -PermIJ*MobIJ*((p_J - p_I)/length)*edgVec
			// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

			// velocity normal component
			double VIJN[3] = {.0,.0,.0};
			//double a = std::inner_product(edIJ.begin(),edIJ.end(),pw_grad_IJ.begin(),.0);
//			double a = edIJ[0]*pw_grad_IJ[0] + edIJ[1]*pw_grad_IJ[1] + edIJ[2]*pw_grad_IJ[2];
			double a = .0;
			for (i=0; i<dim; i++) a += edIJ[i]*pw_grad_IJ[i];
			for (i=0; i<dim; i++) VIJN[i] = pw_grad_IJ[i] - a*edIJ[i];

			// velocity parallel component - finite difference approximation
			a = (pStruct->pPPData->getPressure(J) - pStruct->pPPData->getPressure(I))/length;

			double VIJFD[3] = {edIJ[0]*a, edIJ[1]*a, edIJ[2]*a};

			// total velocity: vIJ
			double aux[3] = {VIJFD[0]+VIJN[0], VIJFD[1]+VIJN[1], VIJFD[2]+VIJN[2]};

			// make v_old = v_new (for t=0, v_new = 0)
			pStruct->pPPData->getVelocity_new(edge,dom,v);
			pStruct->pPPData->setVelocity_old(edge,dom,v);

			// update v_new
			dblarray vel(dim,.0);

			// get absolute permeability tensor: dim x dim
			const double *K = pStruct->pSimPar->getPermeability(dom);

			// get mobility for nodes I and J
			const double MobI = pStruct->pPPData->getTotalMobility(I);
			const double MobJ = pStruct->pPPData->getTotalMobility(J);
			const double MobIJ = .5*(MobI + MobJ);

			int pos = 0;
			if ( pStruct->pSimPar->is_K_Isotropic() ){
				for (i=0; i<dim; i++){
					vel[i] = K[pos]*aux[i];
					pos += dim+1;
				}
			}
			else{
				for (i=0; i<dim; i++)
					for (j=0; j<dim; j++) vel[i] += K[dim*i+j]*aux[j];
			}
			for (i=0; i<dim; i++) vel[i] *= -MobIJ;

#ifdef _SEEKFORBUGS_
			for (i=0; i<dim; i++)
				if ( fabs(vel[i]) > 0.0 ) check = true;
#endif

			// total velocity: VTIJ = VIJFD + VIJN
			pStruct->pPPData->setVelocity_new(edge,dom,vel);
			//char text[256];
//			sprintf(text,"vel[%d - %d]: %.8E\t%.8E\tgrad_mean: %.8E\t%.8E\n",EN_id(I),EN_id(J),vel[0],vel[1],pw_grad_IJ[0],pw_grad_IJ[1]);
//			sprintf(text,"vel[%d - %d]: VIJFD: %.8E %.8E\n",EN_id(I),EN_id(J),VIJFD[0],VIJFD[1]);
//			fid << text;
		}
	}
	EIter_delete(eit);

//	fid.close();
//	STOP();


#ifdef _SEEKFORBUGS_
	if (!check) throw Exception(__LINE__,__FILE__,"Velocity field null!\n");
#endif

	double endt = MPI_Wtime();
	return endt-startt;
}
}
