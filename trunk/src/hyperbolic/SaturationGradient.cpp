/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic.h"

namespace PRS{
	void EBFV1_hyperbolic::calculateSaturationGradient(){
		CPU_Profile::Start();

		int dim = pGCData->getMeshDim();
		resetSaturationGradient();
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calc_Sw_grad_1(dom,dim);
			calc_Sw_grad_2(dom,dim);
		}
		calc_Sw_grad_3(dim);
		calc_Sw_grad_4(dim);

		CPU_Profile::End("SaturationGradient");
	}

	void EBFV1_hyperbolic::calc_Sw_grad_1(int dom, int dim){
		const double* Cij = NULL;
		double Sw_grad_I[3], Sw_grad_J[3], Sw_I, Sw_J, val;
		int i,j,nedges, edge, idx0, idx1, idx0_global, idx1_global, id0, id1, dom_flag;
		nedges = pGCData->getNumEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getCij(dom,edge,Cij);
			pGCData->getEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pGCData->getID(dom,idx0,idx1,id0,id1);
			if (id0 > id1){
				std::swap(idx0_global,idx1_global);
				std::swap(idx0,idx1);
			}
			pPPData->getSaturation(idx0_global,Sw_I);
			pPPData->getSaturation(idx1_global,Sw_J);
			pPPData->get_Sw_Grad(idx0,Sw_grad_I);
			pPPData->get_Sw_Grad(idx1,Sw_grad_J);
			val = 0.5*(Sw_I + Sw_J);
			for (i=0; i<dim; i++){
				Sw_grad_I[i] += val*Cij[i];
				Sw_grad_J[i] += -val*Cij[i];
			}
			pPPData->set_Sw_Grad(idx0,Sw_grad_I);
			pPPData->set_Sw_Grad(idx1,Sw_grad_J);
		}
	}

	void EBFV1_hyperbolic::calc_Sw_grad_2(int dom, int dim){
		const double *Dij = NULL;
		double Sw_grad_I[3], Sw_grad_J[3], Sw_grad_K[3], Sw_I, Sw_J, Sw_K;
		int i, j, edge, idx0, idx1, idx2, idx0_global, idx1_global, idx2_global;

		if (dim==2){
			for (edge=0; edge<pGCData->getNumBDRYEdgesPerDomain(dom); edge++){
				pGCData->getDij(dom,edge,Dij);
				pGCData->getBdryEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
				pPPData->getSaturation(idx0_global,Sw_I);
				pPPData->getSaturation(idx1_global,Sw_J);
				pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);
				for (int i=0; i<dim; i++){
					Sw_grad_I[i] += ((5.*Sw_I + Sw_J)/6.)*Dij[i];
					Sw_grad_J[i] += ((Sw_I + 5.*Sw_J)/6.)*Dij[i];
				}
				pPPData->set_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->set_Sw_Grad(idx1_global,Sw_grad_J);
			}
		}
		else{
			double line1[3] = {6., 1., 1.};
			double line2[3] = {1., 6., 1.};
			double line3[3] = {1., 1., 6.};
			double dot1, dot2, dot3, Sw_vec[3], aux[3];
			for (int face=0; face<pGCData->getNumBdryFacesPerDomain(dom); face++){
				pGCData->getDij(dom,face,Dij);
				pGCData->getBdryFace(dom,face,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
				pPPData->getSaturation(idx0_global,Sw_I);
				pPPData->getSaturation(idx1_global,Sw_J);
				pPPData->getSaturation(idx2_global,Sw_K);
				pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);
				pPPData->get_Sw_Grad(idx2_global,Sw_grad_K);
				dot1 = inner_product(line1,Sw_vec,3);
				dot2 = inner_product(line2,Sw_vec,3);
				dot3 = inner_product(line3,Sw_vec,3);

				for (int i=0; i<3; i++){
					aux[i] = Dij[i]/8.0;
				}

				for (int i=0; i<3; i++){
					Sw_grad_I[i] += aux[i]*dot1;
					Sw_grad_J[i] += aux[i]*dot2;
					Sw_grad_K[i] += aux[i]*dot3;
				}

				pPPData->set_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->set_Sw_Grad(idx1_global,Sw_grad_J);
				pPPData->set_Sw_Grad(idx2_global,Sw_grad_K);
			}
		}
	}

	void EBFV1_hyperbolic::calc_Sw_grad_3(int dim){
		double Sw_grad_tmp[3], Sw_grad[3], volume;
		int i, dom, ndom, nnodes, node, idx;
		ndom = (int)pSimPar->setOfDomains.size();
		// performe sw_grad accumulation for multidomains
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pGCData->getNodeIdx_Global(dom,node,idx);
				pPPData->get_Sw_Grad(idx,Sw_grad);
				pPPData->get_Sw_Grad(node,Sw_grad_tmp);
				for (i=0; i<dim; i++){
					Sw_grad[i] += Sw_grad_tmp[i];
				}
				pPPData->set_Sw_Grad(idx,Sw_grad);
			}
		}
		// weight cumulated sw_grad per node volume
		nnodes = pGCData->getNumNodes();
		for (node=0; node<nnodes; node++){
			pPPData->get_Sw_Grad(node,Sw_grad);
			pGCData->getVolume(node,volume);
			for (i=0; i<dim; i++){
				Sw_grad[i] /= volume;
			}
			pPPData->set_Sw_Grad(node,Sw_grad);
			//cout << "Sw_grad: " << Sw_grad[0] << " " << Sw_grad[1] << " " << Sw_grad[2] << "\n";
		}
		//STOP();
	}

	void EBFV1_hyperbolic::calc_Sw_grad_4(int dim){
		double Sw_grad_I[3], Sw_grad_J[3], Sw_grad_K[3], versor[3], innerp1, innerp2, Dij[3];
		int i,nedges, edge, idx0_global, idx1_global, idx2_global, flag1, flag2, flag3;

		if (dim==2){
			nedges = pGCData->getNumEBE();
			for (edge=0; edge<nedges; edge++){
				pGCData->getVersor_ExternalBdryElement(edge,versor);
				pGCData->getExternalBdryEdges(edge,idx0_global,idx1_global,flag1,flag2);
				pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
				pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);

				// Parallel Component of The Saturation Gradient
				innerp1 = inner_product(Sw_grad_I,versor,dim);
				innerp2 = inner_product(Sw_grad_J,versor,dim);

				for (i=0; i<dim; i++){
					Sw_grad_I[i] = innerp1*versor[i];
					Sw_grad_J[i] = innerp2*versor[i];
				}

				// exclude injection wells
				if ( !pSimPar->isInjectionWell(flag1)  && !pPPData->getProjectedSwgrad(idx0_global)){
					pPPData->set_Sw_Grad(idx0_global,Sw_grad_I);
					pPPData->setProjectedSwgrad(idx0_global,true);
				}
				if ( !pSimPar->isInjectionWell(flag2)  && !pPPData->getProjectedSwgrad(idx1_global) ){
					pPPData->set_Sw_Grad(idx1_global,Sw_grad_J);
					pPPData->setProjectedSwgrad(idx1_global,true);
				}
			}
		}
		else{
			int idx0,idx1,idx2;
			int ndom = (int)pSimPar->setOfDomains.size();
			for (int dom=0; dom<ndom; dom++){
				int nfaces = pGCData->getNumBdryFacesPerDomain(dom);
				double line1[3] = {6., 1., 1.};
				double line2[3] = {1., 6., 1.};
				double line3[3] = {1., 1., 6.};
				double dot1, dot2, dot3, Sw_vec[3], aux[3];
				for (int face=0; face<nfaces; face++){
					pGCData->getDij(dom,face,Dij);
					pGCData->getBdryFace(dom,face,idx0,idx1,idx2,idx0_global,idx1_global,idx2_global);
					pPPData->get_Sw_Grad(idx0_global,Sw_grad_I);
					pPPData->get_Sw_Grad(idx1_global,Sw_grad_J);
					pPPData->get_Sw_Grad(idx2_global,Sw_grad_K);

					// project each nodal gradient on face
					double norma = sqrt( inner_product(Dij,Dij,dim) );
					double n[3] = {Dij[0]/norma, Dij[1]/norma, Dij[2]/norma};

					double scalar[3] = {.0,.0,.0};
					for (i=0; i<3; i++){
						scalar[0] += Sw_grad_I[i]*n[i];
						scalar[1] += Sw_grad_J[i]*n[i];
						scalar[2] += Sw_grad_K[i]*n[i];
					}

					for (i=0; i<3; i++){
						Sw_grad_I[i] = Sw_grad_I[i] - scalar[0]*n[i];
						Sw_grad_J[i] = Sw_grad_J[i] - scalar[1]*n[i];
						Sw_grad_K[i] = Sw_grad_K[i] - scalar[2]*n[i];
					}

					// exclude injection wells
					if ( !pSimPar->isInjectionWell(flag1)  && !pPPData->getProjectedSwgrad(idx0_global)){
						pPPData->set_Sw_Grad(idx0_global,Sw_grad_I);
						pPPData->setProjectedSwgrad(idx0_global,true);
					}
					if ( !pSimPar->isInjectionWell(flag2)  && !pPPData->getProjectedSwgrad(idx1_global) ){
						pPPData->set_Sw_Grad(idx1_global,Sw_grad_J);
						pPPData->setProjectedSwgrad(idx1_global,true);
					}
				}
			}
		}
	}

	void EBFV1_hyperbolic::resetSaturationGradient(){
		double Sw_grad[3]={.0,.0,.0};
		int nnodes = pGCData->getNumNodes();
		for (int node=0; node<nnodes; node++){
			pPPData->set_Sw_Grad(node,Sw_grad);
			pPPData->setProjectedSwgrad(node,false);
		}
	}
}
