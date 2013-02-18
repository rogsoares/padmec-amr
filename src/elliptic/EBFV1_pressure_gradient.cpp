#include "EBFV1_elliptic.h"

namespace PRS{

void calculateGradient__bdryFaces(pMesh theMesh, SimulatorParameters *pSimPar, PhysicPropData* pPPData, GeomData *pGCData, int dom, int dom_counter, char* tag);

double EBFV1_elliptic::pressureGradient(){

#ifdef _SEEKFORBUGS_
	bool check1 = false;
	bool check2 = false;
#endif

	int i, row_I, row_J;//, row_K;
	int dom_counter = 0;
	int dim = pGCData->getMeshDim();
	pEntity node, edge;//, face;
	dblarray Cij(dim,.0), Dij(dim,.0);
	double pw_grad_I[3], pw_grad_J[3];//,  pw_grad_K[3];
	double pressure_I, pressure_J;//, pressure_K;

	// loop over domains
	SIter_const iter=pSimPar->setDomain_begin();
	for (;iter!=pSimPar->setDomain_end();iter++){
		// tag to get local node ID per domain
		char tag[4]; sprintf(tag,"%d",dom_counter);

		int dom = *iter;	// domain's flag
		// before any calculation, reset all previous values
		resetPressureGradient(dom_counter,tag);

		// loop over all edges
		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			if (!theMesh->getRefinementDepth(edge) && pGCData->edgeBelongToDomain(edge,dom)){

				// get nodes I and J
				pEntity I = (pVertex)edge->get(0,0);
				pEntity J = (pVertex)edge->get(0,1);

				// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
				if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
					std::swap(I,J);
				}

				// local node IDs numbering
				pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
				pSimPar->getLocalNodeIDNumbering(J,tag,row_J);

				// get Cij vector (normal to control volume surface)
				pGCData->getCij(edge,dom,Cij);

				// get nodal pressure gradient
				pPPData->get_pw_Grad(dom_counter,row_I,pw_grad_I);
				pPPData->get_pw_Grad(dom_counter,row_J,pw_grad_J);

				// get nodal pressure
				pressure_I = pPPData->getPressure(I);
				pressure_J = pPPData->getPressure(J);

#ifdef _SEEKFORBUGS_
				if ( fabs(pressure_I) > 0.0 || fabs(pressure_J) > 0.0 ){
					check1=true;
				}
#endif

				double nrc = (double)pGCData->getNumRC(theMesh,edge) + 1.0;
				double val = 0.5*(pressure_I + pressure_J)/nrc;
				for (i=0; i<dim; i++){
					pw_grad_I[i] +=  val*Cij[i];
					pw_grad_J[i] += -val*Cij[i];
				}

				// update nodal pressure gradient
				pPPData->set_pw_Grad(dom_counter,row_I,pw_grad_I);
				pPPData->set_pw_Grad(dom_counter,row_J,pw_grad_J);
			}
		}
		EIter_delete(eit);


#ifdef _SEEKFORBUGS_
		if (!check1) {
			char msg[256]; sprintf(msg,"Pressure field null for domain %d.\n",dom);
			throw Exception(__LINE__,__FILE__,msg);
		}
#endif

		// edge (boundary) contribution to pressure gradient (2D meshes)
		if (dim==2){
			eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				// get edge flag
				//int flag = EN_getFlag(edge);

				// get only edges on domain's boundaries
				if ( !theMesh->getRefinementDepth(edge) && pGCData->belongsToBoundary(edge) )
					if ( pGCData->edgeBelongToDomain(edge,dom) ){
						Dij[0] = .0; Dij[1] = .0;
						pGCData->getDij(edge,dom,Dij);

						// get nodes I and J
						pEntity I = (pVertex)edge->get(0,0);
						pEntity J = (pVertex)edge->get(0,1);
						// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
						if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
							std::swap(I,J);
						}
						// get nodal pressure
						pressure_I = pPPData->getPressure(I);
						pressure_J = pPPData->getPressure(J);
						// local node IDs numbering
						pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
						pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
						// get nodal pressure gradient
						pPPData->get_pw_Grad(dom_counter,row_I,pw_grad_I);
						pPPData->get_pw_Grad(dom_counter,row_J,pw_grad_J);
						for (i=0; i<dim; i++){
							pw_grad_I[i] += ((5.*pressure_I + pressure_J)/6.0)*Dij[i];
							pw_grad_J[i] += ((pressure_I + 5.*pressure_J)/6.0)*Dij[i];
						}
						// update nodal pressure gradient
						pPPData->set_pw_Grad(dom_counter,row_I,pw_grad_I);
						pPPData->set_pw_Grad(dom_counter,row_J,pw_grad_J);
					}
			}
			EIter_delete(eit);
		}// end of loop over bdry edges
		else{
			// faces contribution to pressure gradient (3D meshes)
			calculateGradient__bdryFaces(theMesh,pSimPar,pPPData,pGCData,dom,dom_counter,tag);

		} // end of loop over bdry triangles

		// weighting by domain's volume
		double vol;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			if ( pGCData->nodeBelongToDomain(node,dom) ){
				//cout << EN_id(node) << ":\t";
				vol = pGCData->getVolume(node,dom);
				//pPPData->get_pw_Grad(node,dom,pw_grad_I);
				pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
				pPPData->get_pw_Grad(dom_counter,row_I,pw_grad_I);
				for (i=0; i<dim; i++) {
					pw_grad_I[i] /= vol;
					//cout << pw_grad_I[i] << "\t";
#ifdef _SEEKFORBUGS_
					if ( fabs(pw_grad_I[i]) > 0.0 ) check2 = true;
#endif
				}
				//cout << endl;
				pPPData->set_pw_Grad(dom_counter,row_I,pw_grad_I);
			}
		}
		VIter_delete(vit);
		//STOP();

		/*
		 *  Calculate pressure gradient on nodes on partition boundaries.
		 *  Only for parallel.
		 */
		//pMData->unifyVectorsOnMeshNodes(pPPData->get_pw_Grad2,pPPData->set_pw_Grad2,pGCData,dom);

#ifdef _SEEKFORBUGS_
		if (!check1) throw Exception(__LINE__,__FILE__,"Pressure field null!\n");
		if (!check2) throw Exception(__LINE__,__FILE__,"Gradient null!\n");
#endif

		dom_counter++;
	} // end of loop over domains
	return 0;
}

int EBFV1_elliptic::resetPressureGradient(int dom, char *tag){
	pEntity node;
	int row_I;
	double pw_grad_I[3] = {.0,.0,.0};
	std::set<int> setNodes;
	pSimPar->getNodesDomain(dom,setNodes);
	std::set<int>::iterator iter = setNodes.begin();
	for(;iter!=setNodes.end();iter++){
		node = theMesh->getVertex(*iter);
		pSimPar->getLocalNodeIDNumbering(node,tag,row_I);
		pPPData->set_pw_Grad(dom,row_I,pw_grad_I);
	}
	return 0;
}

void calculateGradient__bdryFaces(pMesh theMesh, SimulatorParameters *pSimPar, PhysicPropData* pPPData, GeomData *pGCData, int dom, int dom_counter, char* tag){
	int row_I, row_J, row_K;
	double pw_grad_I[3], pw_grad_J[3],  pw_grad_K[3], Dij[3];
	pEntity face;
	FIter fit = M_faceIter(theMesh);
	while ( (face = FIter_next(fit)) ){
		if ( !theMesh->getRefinementDepth(face) && pGCData->getDij(face,dom,Dij) ){

			// get nodes I, J and K
			pEntity I = (pVertex)face->get(0,0);
			pEntity J = (pVertex)face->get(0,1);
			pEntity K = (pVertex)face->get(0,2);

			pSimPar->getLocalNodeIDNumbering(I,tag,row_I);
			pSimPar->getLocalNodeIDNumbering(J,tag,row_J);
			pSimPar->getLocalNodeIDNumbering(K,tag,row_K);

			// get nodal pressure
			double pressure_I = pPPData->getPressure(I);
			double pressure_J = pPPData->getPressure(J);
			double pressure_K = pPPData->getPressure(K);

#ifdef _SEEKFORBUGS_
			//if ( fabs(pressure_I) > 0.0 || fabs(pressure_J) > 0.0 || fabs(pressure_K) > 0.0 ) check1=true;
#endif

			pPPData->get_pw_Grad(dom_counter,row_I,pw_grad_I);
			pPPData->get_pw_Grad(dom_counter,row_J,pw_grad_J);
			pPPData->get_pw_Grad(dom_counter,row_K,pw_grad_K);
			for (int i=0; i<3; i++){
				pw_grad_I[i] += ((6.*pressure_I + pressure_J + pressure_K)/8.0)*Dij[i];
				pw_grad_J[i] += ((pressure_I + 6.*pressure_J + pressure_K)/8.0)*Dij[i];
				pw_grad_K[i] += ((pressure_I + pressure_J + 6.*pressure_K)/8.0)*Dij[i];
			}
			pPPData->set_pw_Grad(dom_counter,row_I,pw_grad_I);
			pPPData->set_pw_Grad(dom_counter,row_J,pw_grad_J);
			pPPData->set_pw_Grad(dom_counter,row_K,pw_grad_K);
		}
	}
	FIter_delete(fit);
}
}
