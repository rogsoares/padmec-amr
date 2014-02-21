/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{
	void GeomData::calculateNumFaces(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		int dim = theMesh->getDim();
		numDomFaces = new int[ndom];
		pEntity face;

		for (int i=0; i<ndom; i++){
			int counter = 0;
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				int flag1 = getFaceFlag(face);
				int flag2 = domainList[i];
				if (flag1==flag2){
					counter++;
				}
			}
			FIter_delete(fit);
			numDomFaces[i] = counter;
		}
	}

	void GeomData::calculateNumFacesTmp(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		int dim = theMesh->getDim();
		numDomFaces_tmp = new int[ndom];
		pEntity face;

		for (int i=0; i<ndom; i++){
			int counter = 0;
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domainList[i]){
					counter++;
				}
			}
			FIter_delete(fit);
			numDomFaces_tmp[i] = counter;
		}
	}

	void GeomData::calculateNumNodes(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		int dim = theMesh->getDim();
		setMeshNodes(M_numVertices(theMesh));
		numNodesPerDomain = new int[ndom];
		pEntity face, tetra;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					if (getFaceFlag(face)==domainList[i]){
						for (int j = 0; j<3; j++){
							nodeSet.insert(face->get(0,j));
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<4; j++){
							nodeSet.insert(tetra->get(0,j));
						}
					}
				}
				RIter_delete(rit);
			}
			numNodesPerDomain[i] = (int)nodeSet.size();
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumBdryNodes(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		int dim = theMesh->getDim();
		numBdryNodesPerDomain = new int[ndom];
		pEntity face, tetra;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					int faceflag = getFaceFlag(face);
					if (faceflag==domainList[i]){
						for (int j = 0; j<3; j++){
							if (getVertexFlag(face->get(0,j))!=faceflag){
								nodeSet.insert(face->get(0,j));
							}
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<4; j++){
							if (getVertexFlag(tetra->get(0,j))!=tetra_flag){
								nodeSet.insert(tetra->get(0,j));
							}
						}
					}
				}
				RIter_delete(rit);
			}
			numBdryNodesPerDomain[i] = (int)nodeSet.size();
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumEdges(pMesh theMesh){
		int dim = theMesh->getDim();
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		numDomEdges = new int[ndom];
		//this->domainList = new int[ndom];
		//		for (int i=0; i<ndom; i++){
			//			this->domainList[i] = domainList[i];
			//		}
		pEntity face, tetra;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			if (dim==2){
				FIter fit = M_faceIter(theMesh);
				while ( (face = FIter_next(fit)) ){
					if (getFaceFlag(face)==domainList[i]){
						for (int j = 0; j<3; j++){
							edgeList.insert(face->get(1,j));
						}
					}
				}
				FIter_delete(fit);
			}
			else{
				RIter rit = M_regionIter(theMesh);
				while ( (tetra = RIter_next(rit)) ){
					int tetra_flag = getTetraFlag(tetra);
					if (tetra_flag==domainList[i]){
						for (int j = 0; j<4; j++){
							edgeList.insert(tetra->get(1,j));
						}
					}
				}
				RIter_delete(rit);
			}
			numDomEdges[i] = (int)edgeList.size();
			edgeList.clear();
		}
	}

	void GeomData::calculateNumBDRYEdges(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		numDomBDRYEdges = new int[ndom];
		pEntity edge, face;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domainList[i]){
					for (int j = 0; j<3; j++){
						edge = (pEdge)face->get(1,j);
						if ( getEdgeFlag(edge)!= domainList[i]){
							edgeList.insert(edge);
						}
					}
				}
			}
			FIter_delete(fit);
			numDomBDRYEdges[i] = (int)edgeList.size();
			edgeList.clear();
		}
		numExtBdryEdges = 0;
		EIter eit = M_edgeIter(theMesh);
		while ( (edge=EIter_next(eit)) ){
			if (E_numFaces(edge)==1){
				numExtBdryEdges++;
			}
		}
		EIter_delete(eit);
	}
	void GeomData::calculateEdgeProperties(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domlist = getDomainList();
		double c1[3], c2[3], vec[3];
		double elength, Lij[3];
		double Icoord[3], Jcoord[3];
		double delta_x = 1e30;
		int i,j;
		int dim = theMesh->getDim();
		pEntity edge;
		for (int dom = 0; dom<ndom; dom++){
			int row = 0;
			int dom_flag = domainList[dom];
			EIter eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				if (edgeBelongToDomain(edge,dom_flag)){
					E_getVerticesCoord(edge,Icoord,Jcoord);
					for (j=0; j<dim; j++){
						Lij[j] = .0;
					}
					makeVector(Jcoord,Icoord,Lij);
					for (i=0; i<dim; i++){
						vec[i] = Lij[i];
					}
					elength = .0;
					for (i=0; i<dim; i++){
						elength += vec[i]*vec[i];
					}
					elength = sqrt(elength);
					if (elength == .0){
						char msg[256]; sprintf(msg,"Edge [%d %d] has null length!",EN_id(edge->get(0,0)),EN_id(edge->get(0,1)));
						throw Exception(__LINE__,__FILE__,msg);
					}
					edge_length[dom].setValue(row,elength);
					for (i=0; i<dim; i++){
						vec[i] /= elength;
					}
					edge_versor[dom].setValue(row,0,vec[0]);
					edge_versor[dom].setValue(row,1,vec[1]);
					edge_versor[dom].setValue(row,2,vec[2]);
					if (delta_x > elength){
						delta_x = elength;
					}
					row++;
				}
			}
			EIter_delete(eit);
		}
		setSmallestEdgeLength(P_getMinDbl(delta_x));
		cout << setprecision(7) << "smalest = " << delta_x << endl;

		int tnedges = 0;
		int row = 0;
		EIter eit = M_edgeIter(theMesh);
		while ( (edge = EIter_next(eit)) ){
			tnedges++;
			if (E_numFaces(edge)==1){
				E_getVerticesCoord(edge,Icoord,Jcoord);
				for (j=0; j<dim; j++){
					Lij[j] = .0;
				}
				makeVector(Jcoord,Icoord,Lij);
				for (i=0; i<dim; i++){
					vec[i] = Lij[i];
				}
				elength = .0;
				for (i=0; i<dim; i++){
					elength += vec[i]*vec[i];
				}
				elength = sqrt(elength);

				for (i=0; i<dim; i++){
					vec[i] /= elength;
				}
				EBE_1[0].setValue(row,0,vec[0]);
				EBE_1[0].setValue(row,1,vec[1]);
				EBE_1[0].setValue(row,2,vec[2]);
				row++;
			}
		}
		EIter_delete(eit);
		this->setTotalNumberOfEdges(tnedges);
	}
}
