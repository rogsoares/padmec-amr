/*
 * gradientSaturation.cpp
 *
 *  Created on: 14/05/2009
 *      Author: rogerio
 */

#include "EBFV1_elliptic.h"

namespace PRS{
	void EBFV1_elliptic::calculatePressureGradient(){
		int dim = pGCData->getMeshDim();
		resetPressureGradient();
		int ndom = (int)pSimPar->setOfDomains.size();
		for (int dom=0; dom<ndom; dom++){
			calc_p_grad_1(dom,dim);
			calc_p_grad_2(dom,dim);
		}
		calc_p_grad_3(dim);
	}

	void EBFV1_elliptic::calc_p_grad_1(int dom, int dim){
		double Cij[3], p_grad_I[3], p_grad_J[3], p_I, p_J, val;
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
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad(dom,idx0,p_grad_I);
			pPPData->get_pw_Grad(dom,idx1,p_grad_J);
			val = 0.5*(p_I + p_J);
			for (i=0; i<dim; i++){
				p_grad_I[i] += val*Cij[i];
				p_grad_J[i] += -val*Cij[i];
			}
			pPPData->set_pw_Grad(dom,idx0,p_grad_I);
			pPPData->set_pw_Grad(dom,idx1,p_grad_J);
		}
	}

	void EBFV1_elliptic::calc_p_grad_2(int dom, int dim){
		double Dij[3], p_grad_I[3], p_grad_J[3], p_I, p_J;
		int i,j,nedges, edge, idx0, idx1,idx0_global, idx1_global;
		nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
		for (edge=0; edge<nedges; edge++){
			pGCData->getDij(dom,edge,Dij);
			pGCData->getBdryEdge(dom,edge,idx0,idx1,idx0_global,idx1_global);
			pPPData->getPressure(idx0_global,p_I);
			pPPData->getPressure(idx1_global,p_J);
			pPPData->get_pw_Grad(dom,idx0,p_grad_I);
			pPPData->get_pw_Grad(dom,idx1,p_grad_J);
			for (int i=0; i<dim; i++){
				p_grad_I[i] += ((5.*p_I + p_J)/6.)*Dij[i];
				p_grad_J[i] += ((p_I + 5.*p_J)/6.)*Dij[i];
			}
			pPPData->set_pw_Grad(dom,idx0,p_grad_I);
			pPPData->set_pw_Grad(dom,idx1,p_grad_J);
		}
	}

	void EBFV1_elliptic::calc_p_grad_3(int dim){
		double p_grad_tmp[3], p_grad[3], vol;
		int i, dom, ndom, nnodes, node, idx;
		ndom = (int)pSimPar->setOfDomains.size();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->get_pw_Grad(dom,node,p_grad);
				pGCData->getVolume(dom,node,vol);
				for (i=0; i<dim; i++){
					p_grad[i] /= vol;
				}
				pPPData->set_pw_Grad(dom,node,p_grad);
			}
		}
	}

	void EBFV1_elliptic::resetPressureGradient(){
		double p_grad[3]={.0,.0,.0};
		int dom, ndom, nnodes, node;
		ndom = (int)pSimPar->setOfDomains.size();
		for (dom=0; dom<ndom; dom++){
			nnodes = pGCData->getNumNodesPerDomain(dom);
			for (node=0; node<nnodes; node++){
				pPPData->set_pw_Grad(dom,node,p_grad);
			}
		}
	}
}
