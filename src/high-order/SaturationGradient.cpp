/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{
void printSaturationGradient(pMesh theMesh, PointerStruct* pStruct, GeomData* pGCData );

double EBFV1_hyperbolic::calculateSaturationGradient(pMesh theMesh){
	double start = MPI_Wtime();
	int dim = pGCData->getMeshDim();
	resetSaturationGradient(theMesh);
	int dom_counter = 0;
	for (SIter_const dom= pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end();dom++){
		calc_Sw_grad_1(theMesh,*dom,dom_counter,dim);
		calc_Sw_grad_2(theMesh,*dom,dim);
		dom_counter++;
	}
	calc_Sw_grad_31(theMesh);

	// compute saturation gradient on partition boundary nodes
	//	pMData->unifyVectorsOnMeshNodes(pStruct->pPPData->get_Sw_Grad,pStruct->pPPData->set_Sw_Grad,pGCData,3300);
	calc_Sw_grad_4(theMesh,dim);
	double end = MPI_Wtime();
	return end-start;
}

void EBFV1_hyperbolic::calc_Sw_grad_1(pMesh theMesh, int dom, int dom_counter, int dim){
	pEntity edge, I, J;
	//dblarray Cij(dim,.0);
	//int dom_counter = 0;
	int row_I, row_J;
	double Sw_grad_I[3]={.0,.0,.0}, Sw_grad_J[3]={.0,.0,.0};
	char tag[4]; sprintf(tag,"%d",dom_counter);

	double Cij[3];
	int row = 0;
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		if ( pGCData->edgeBelongToDomain(edge,dom) ){
			// get Cij vector (normal to control volume surface)
			//pGCData->getCij(edge,dom,Cij);
			pGCData->getCij(dom_counter,row,Cij); row++;

			// get nodes I and J
			I = (pVertex)edge->get(0,0);
			J = (pVertex)edge->get(0,1);

			/// NAO TIRAR ESSE IF DAQUI!!!  ROGERIO: 23/08/2013
			if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
				std::swap(I,J);
			}

			// get saturation gradient for nodes I and J
//			pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
//			pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
			pStruct->pPPData->get_Sw_Grad(I,dom,row_I,Sw_grad_I);
			pStruct->pPPData->get_Sw_Grad(J,dom,row_J,Sw_grad_J);

			// get saturation for nodes I and J
			double Sw_I = pStruct->pPPData->getSaturation(I);
			double Sw_J = pStruct->pPPData->getSaturation(J);

			//double nrc = pStruct->pGCData->getNumRemoteCopies(edge) + 1.0;
			double val = 0.5*(Sw_I + Sw_J);//nrc;
			for (int i=0; i<dim; i++){
				Sw_grad_I[i] += val*Cij[i];
				Sw_grad_J[i] += -val*Cij[i];
			}

			// update saturation
			pStruct->pPPData->set_Sw_Grad(I,dom,row_I,Sw_grad_I);
			pStruct->pPPData->set_Sw_Grad(J,dom,row_J,Sw_grad_J);
		}
	}
	EIter_delete(eit);
}

void EBFV1_hyperbolic::calc_Sw_grad_2(pMesh theMesh, int dom, int dim){
	pEntity edge, I, J, K;
	dblarray Dij(dim,.0);
	double Sw_grad_I[3], Sw_grad_J[3], Sw_grad_K[3];
	int dom_counter = 0;
	int row_I, row_J, row_K;
	char tag[4]; sprintf(tag,"%d",dom_counter);

	// loop over domains
	if (dim==2){
		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			if ( pGCData->edgeBelongToDomain(edge,dom) ){
				if (pGCData->belongsToBoundary(edge)){
					Dij[0] = .0; Dij[1] = .0;
					pGCData->getDij(edge,dom,Dij);
					// get nodes I and J
					I = (pVertex)edge->get(0,0);
					J = (pVertex)edge->get(0,1);

					// get saturation for nodes I and J
					double Sw_I = pStruct->pPPData->getSaturation(I);
					double Sw_J = pStruct->pPPData->getSaturation(J);

					// get saturation gradient for nodes I and J
//					pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
//					pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
					pStruct->pPPData->get_Sw_Grad(I,dom,row_I,Sw_grad_I);
					pStruct->pPPData->get_Sw_Grad(J,dom,row_J,Sw_grad_J);

					for (int i=0; i<dim; i++){
						Sw_grad_I[i] += ((5.*Sw_I + Sw_J)/6.)*Dij[i];
						Sw_grad_J[i] += ((Sw_I + 5.*Sw_J)/6.)*Dij[i];
					}
					pStruct->pPPData->set_Sw_Grad(I,dom,row_I,Sw_grad_I);
					pStruct->pPPData->set_Sw_Grad(J,dom,row_J,Sw_grad_J);
				}
			}
		}
		EIter_delete(eit);
	}
	//	else{
	//		pEntity face;
	//		FIter fit = M_faceIter(theMesh);
	//		while ( (face = FIter_next(fit)) ){
	//			if (pGCData->getDij(face,dom,Dij)){
	//				// get nodes I, J and K
	//				pEntity I = (pVertex)face->get(0,0);
	//				pEntity J = (pVertex)face->get(0,1);
	//				pEntity K = (pVertex)face->get(0,2);
	//
	//				// get nodal pressure
	//				double Sw_I = pStruct->pPPData->getSaturation(I);
	//				double Sw_J = pStruct->pPPData->getSaturation(J);
	//				double Sw_K = pStruct->pPPData->getSaturation(K);
	//
	//				// get nodal gradient
	//				pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
	//				pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
	//				pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_K);
	//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_I,Sw_grad_I);
	//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_J,Sw_grad_J);
	//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_K,Sw_grad_K);
	//
	//				for (int i=0; i<3; i++){
	//					Sw_grad_I[i] += ((6.*Sw_I + Sw_J + Sw_K)/8.0)*Dij[i];
	//					Sw_grad_J[i] += ((Sw_I + 6.*Sw_J + Sw_K)/8.0)*Dij[i];
	//					Sw_grad_K[i] += ((Sw_I + Sw_J + 6.*Sw_K)/8.0)*Dij[i];
	//				}
	//
	//				// update gradient
	//				pStruct->pPPData->set_Sw_Grad(dom_counter,row_I,Sw_grad_I);
	//				pStruct->pPPData->set_Sw_Grad(dom_counter,row_J,Sw_grad_J);
	//				pStruct->pPPData->set_Sw_Grad(dom_counter,row_K,Sw_grad_K);
	//			}
	//		}
	//		FIter_delete(fit);
	//	}
}
//
//void EBFV1_hyperbolic::calc_Sw_grad_3(pMesh theMesh, int dom, int dim){
//	pEntity node;
//	double Sw_grad[3],tmp[3];
//	int dom_counter = 0;
//	int row_I;
//	// weighting by domain's volume
//	double vol;
//	char tag[4]; sprintf(tag,"%d",dom_counter);
//	VIter vit = M_vertexIter(theMesh);
//	while ( (node = VIter_next(vit)) ){
//		if ( pGCData->nodeBelongToDomain(node,dom) ){
//			vol = pGCData->getVolume(node,dom);
////			pStruct->pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
//			pStruct->pPPData->get_Sw_Grad(node,dom,row_I,Sw_grad);
//			for (int i=0; i<dim; i++) {
//				Sw_grad[i] /= vol;
//#ifdef _SEEKFORBUGS_
//				if (vol == .0){
//					throw Exception(__LINE__,__FILE__,"Null volume!");
//				}
//#endif
//			}
//			pStruct->pPPData->set_Sw_Grad(node,dom,row_I,Sw_grad);
//		}
//	}
//	VIter_delete(vit);
//}

void EBFV1_hyperbolic::calc_Sw_grad_31(pMesh theMesh){
	pEntity node;
	int dim = theMesh->getDim();
	double tmp[3] = {.0,.0,.0};
	int row_I;
	char tag[4]; sprintf(tag,"%d",0);
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		double Sw_grad[3] = {.0,.0,.0};
		double vol = .0;
		for (SIter_const dom = pStruct->pSimPar->setDomain_begin(); dom != pStruct->pSimPar->setDomain_end(); dom++){
			pStruct->pPPData->get_Sw_Grad(node,*dom,row_I,tmp);
			for (int i=0; i<dim; i++){
				Sw_grad[i] += tmp[i];
			}
			vol += pGCData->getVolume(node,*dom);
		}
		for (int i=0; i<dim; i++){
			Sw_grad[i] /= vol;
		}
		pStruct->pPPData->set_Sw_Grad(node,0,row_I,Sw_grad);
	}
	VIter_delete(vit);
}

void EBFV1_hyperbolic::calc_Sw_grad_4(pMesh theMesh, int dim){
	pEntity edge, I, J, K;
	double Sw_grad_I[3], Sw_grad_J[3], Sw_grad_K[3];
	int dom_counter = 0;
	int row_I, row_J, row_K;
	char tag[4]; sprintf(tag,"%d",dom_counter);

//	SIter_const dom = pStruct->pSimPar->setDomain_begin();
//	for (; dom != pStruct->pSimPar->setDomain_end(); dom++){
		if (dim==2){
			EIter eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				int flag = GEN_tag(edge->getClassification());
				// get only edges on domain's boundaries
				if ( pGCData->belongsToBoundary(edge) && E_numFaces(edge)==1){

					// get nodes I and J
					I = (pVertex)edge->get(0,0);
					J = (pVertex)edge->get(0,1);

					// Nodal Coordinates
					VPoint point_I = E_getVertexPoint(edge,0);
					VPoint point_J = E_getVertexPoint(edge,1);

					// Edge Cossines/Projections
					double delt[2] = {point_J(0)-point_I(0), point_J(1)-point_I(1)};

					// Edge Length
					double length = E_length(edge);

					// Paralell Edge Direction (Unitary vector in edge's direction)
					double versorp[2] = {delt[0]/length,delt[1]/length};

					// Flux Boundary Conditions
					pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
					pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
					pStruct->pPPData->get_Sw_Grad(I,0,row_I,Sw_grad_I);
					pStruct->pPPData->get_Sw_Grad(J,0,row_J,Sw_grad_J);

					// Parallel Component of The Saturation Gradient
					double innerprod1 = Sw_grad_I[0]*versorp[0] + Sw_grad_I[1]*versorp[1];
					double innerprod2 = Sw_grad_J[0]*versorp[0] + Sw_grad_J[1]*versorp[1];
					for (int i=0; i<2; i++){
						Sw_grad_I[i] = innerprod1*versorp[i];
						Sw_grad_J[i] = innerprod2*versorp[i];
					}

					// exclude injection wells
					int flags[2] = {GEN_tag(I->getClassification()), GEN_tag(J->getClassification())};
					if ( !pStruct->pSimPar->isInjectionWell(flags[0]) ){
						pStruct->pPPData->set_Sw_Grad(I,0,row_I,Sw_grad_I);
					}

					if ( !pStruct->pSimPar->isInjectionWell(flags[1]) ){
						pStruct->pPPData->set_Sw_Grad(J,0,row_J,Sw_grad_J);
					}
				}
			}
			EIter_delete(eit);
		//}
		//		else{
		//			int i;
		//			pEntity face;
		//			dblarray Dij(dim,.0);
		//			FIter fit = M_faceIter(theMesh);
		//			while ( face = FIter_next(fit) ){
		//				// get nodes I, J and K
		//				pEntity I = (pVertex)face->get(0,0);
		//				pEntity J = (pVertex)face->get(0,1);
		//				pEntity K = (pVertex)face->get(0,2);
		//
		//				// get saturation gradient for nodes I, J, K
		//				pStruct->pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
		//				pStruct->pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
		//				pStruct->pSimPar->getLocalNodeIDNumbering(K,tag,row_K);
		//
		//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_I,Sw_grad_I);
		//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_J,Sw_grad_J);
		//				pStruct->pPPData->get_Sw_Grad(dom_counter,row_K,Sw_grad_K);
		//
		//				// get vector Dij orthogonal to face IJK
		//				pGCData->getDij(face,*dom,Dij);
		//
		//				// project each nodal gradient on face
		//				double norma = sqrt(Dij[0]*Dij[0] + Dij[1]*Dij[1] + Dij[2]*Dij[2]);
		//				double n[3] = {Dij[0]/norma, Dij[1]/norma, Dij[2]/norma};
		//
		//				double scalar[3] = {.0,.0,.0};
		//				for (i=0; i<3; i++){
		//					scalar[0] += Sw_grad_I[i]*n[i];
		//					scalar[1] += Sw_grad_J[i]*n[i];
		//					scalar[2] += Sw_grad_K[i]*n[i];
		//				}
		//
		//				for (i=0; i<3; i++){
		//					Sw_grad_I[i] = Sw_grad_I[i] - scalar[0]*n[i];
		//					Sw_grad_J[i] = Sw_grad_J[i] - scalar[1]*n[i];
		//					Sw_grad_K[i] = Sw_grad_K[i] - scalar[2]*n[i];
		//				}
		//
		//				// exclude injection wells
		//				int flags[3] = {GEN_tag(I->getClassification()), GEN_tag(J->getClassification()), GEN_tag(K->getClassification())};
		//
		//				if ( !pStruct->pSimPar->isInjectionWell(flags[0]) ){
		//					pStruct->pPPData->set_Sw_Grad(dom_counter,row_I,Sw_grad_I);
		//				}
		//				if ( !pStruct->pSimPar->isInjectionWell(flags[1]) ){
		//					pStruct->pPPData->set_Sw_Grad(dom_counter,row_J,Sw_grad_J);
		//				}
		//				if ( !pStruct->pSimPar->isInjectionWell(flags[2]) ){
		//					pStruct->pPPData->set_Sw_Grad(dom_counter,row_I,Sw_grad_K);
		//				}
		//			}
		//			FIter_delete(fit);
		//		}
	}
}

void EBFV1_hyperbolic::resetSaturationGradient(pMesh theMesh){
	double Sw_grad_I[3] = {.0,.0,.0};
	int dom_counter = 0;
	int row_I;
	char tag[4]; sprintf(tag,"%d",dom_counter);

	SIter_const dom = pStruct->pSimPar->setDomain_begin();
	for (; dom != pStruct->pSimPar->setDomain_end(); dom++){
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			//pStruct->pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
			pStruct->pPPData->set_Sw_Grad(node,*dom,row_I,Sw_grad_I);
		}
		VIter_delete(vit);
	}
	VIter vit = M_vertexIter(theMesh);
	while (pEntity node = VIter_next(vit)){
		pStruct->pPPData->set_Sw_Grad(node,0,0,Sw_grad_I);
	}
	VIter_delete(vit);
}

void printSaturationGradient(pMesh theMesh, PointerStruct* pStruct, GeomData* pGCData){
	double Sw_grad_I[3] = {.0,.0,.0};
	int dom_counter = 0;
	int row_I;
	char tag[4]; sprintf(tag,"%d",dom_counter);

	static int steps = 0;
	char file[256]; sprintf(file,"Sw_grad__step-%d.txt",++steps);
	ofstream fid;
	fid.open(file);
	VIter vit = M_vertexIter(theMesh);
	while (pEntity node = VIter_next(vit)){
		pStruct->pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
		pStruct->pPPData->get_Sw_Grad(node,dom_counter,row_I,Sw_grad_I);
		double vol = .0;
		SIter_const dom = pStruct->pSimPar->setDomain_begin();
		for (; dom != pStruct->pSimPar->setDomain_end(); dom++){
			vol += pGCData->getVolume(node,*dom);
		}

		fid << "ID: " << EN_id(node)
							<< "  pointer: " << node
							<< "\tSw = " << pStruct->pPPData->getSaturation(node)
							<< "\tV = " << vol
							<< "\tSw_grad: " << Sw_grad_I[0] << " " << Sw_grad_I[1] << endl;
	}
	VIter_delete(vit);
	fid.close();
}
}
