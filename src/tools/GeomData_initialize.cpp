/*
 * GeomData_initialize.cpp
 *
 *  Created on: Feb 11, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{

	/*
	 * Based on pre-processor data, matrix structure is created and data transferred
	 */
	void GeomData::initilize(pMesh theMesh, const std::set<int> &setOfDomain){
		// number of domains stay the same no matter how many adaptations happen
		int i = 0;
		int ndom = (int)setOfDomain.size();
		int *domainList = new int[ndom];
		for(set<int>::iterator iter = setOfDomain.begin(); iter!=setOfDomain.end(); iter++){
			domainList[i++] =  *iter;
		}
		setNumDomains(ndom);
		setDomainList(domainList);
		delete[] domainList; domainList = 0;
		//every new mesh adaptation, mesh data structure changes
		initilize(theMesh);
	}

	void GeomData::initilize(pMesh theMesh){
		calculateNumEdges(theMesh);					// calculate number of data to be stored
		calculateNumFaces(theMesh);
		calculateNumBDRYEdges(theMesh);
		calculateNumNodes(theMesh);
		calculateNumBdryNodes(theMesh);
		allocatePointers(M_numVertices(theMesh));	// allocate storage
	}

	void GeomData::dataTransfer(pMesh theMesh){
		calculateEdgeProperties(theMesh);			// fill storage
		transferCijData(theMesh);
		transferDijData(theMesh);
		transferVolData(theMesh);
		mappingNodesIds(theMesh);					// map data to find them quickly
	}

	void GeomData::transferCijData(pMesh theMesh){
		dblarray Cij(3,.0);
		double cij[3], norm;
		int ndom = getNumDomains();
		for (int i=0; i<ndom; i++){
			int dom = domainList[i];
			int row = 0;
			EIter eit = M_edgeIter(theMesh);
			while ( pEdge edge = EIter_next(eit) ){
				if ( edgeBelongToDomain(edge,dom) ){
					getCij(edge,dom,Cij);
					norm = getCij_norm(edge,dom);
					cij[0] = Cij[0];
					cij[1] = Cij[1];
					cij[2] = Cij[2];
					setCij(i,row,cij);
					setCij_norm(i,row,norm);
					row++;
				}
			}
			EIter_delete(eit);
		}
	}

	void GeomData::transferDijData(pMesh theMesh){
		dblarray Dij(3,.0);
		double dij[3];
		int ndom = getNumDomains();
		for (int i=0; i<ndom; i++){
			int dom = domainList[i];
			int row = 0;
			EIter eit = M_edgeIter(theMesh);
			while ( pEdge edge = EIter_next(eit) ){
				if ( edgeBelongToDomain(edge,dom) ){
					if (belongsToBoundary(edge)){
						getDij(edge,dom,Dij);
						dij[0] = Dij[0];
						dij[1] = Dij[1];
						dij[2] = Dij[2];
						setDij(i,row,dij);
						row++;
					}
				}
			}
			EIter_delete(eit);
		}
	}

	void GeomData::transferVolData(pMesh theMesh){
		int idx = 0;
		pEntity node;
		double volume;
		int ndom = getNumDomains();
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int id = EN_id(node);
			volume = .0;
			for (int i=0; i<ndom; i++){
				volume += getVolume(node,domainList[i]);
			}
			setVolume(idx,volume);
			idx++;
		}
		VIter_delete(vit);
	}

	void GeomData::allocatePointers(int nnodes_global){
		int ndom = getNumDomains();
		Cij = new Matrix<double>[ndom];
		Cij_norm = new Matrix<double>[ndom];
		Dij = new Matrix<double>[ndom];
		volume = new Matrix<double>[ndom];
		edges = new Matrix<int>[ndom];
		faces = new Matrix<int>[ndom];
		nodes = new Matrix<int>[ndom];
		ID = new Matrix<int>[ndom];
		volume_bdry = new Matrix<double>[ndom];
		edges_bdry = new Matrix<int>[ndom];
		faces_bdry = new Matrix<int>[ndom];
		ID_bdry = new Matrix<int>[ndom];
		edge_length = new Matrix<double>[ndom];
		edge_versor = new Matrix<double>[ndom];
		volume_global = new Matrix<double>[1];
		volume_global[0].allocateMemory(nnodes_global);
		volume_global[0].initialize(.0);
		EBE_1 = new Matrix<double>[1];
		EBE_2 = new Matrix<int>[1];
		EBE_1[0].allocateMemory(numExtBdryEdges,3);
		EBE_2[0].allocateMemory(numExtBdryEdges,4);
		for (int k=0; k<ndom; k++){
			int nedges = this->numDomEdges[k];
			int nbedges = this->numDomBDRYEdges[k];
			int nnodes = this->numNodesPerDomain[k];
			int nbnodes = this->numBdryNodesPerDomain[k];
			int nfaces = this->numDomFaces[k];

			Cij[k].allocateMemory(nedges,3);
			Cij[k].initialize(.0);
			edges[k].allocateMemory(nedges,6);
			edges[k].initialize(0);
			faces[k].allocateMemory(nfaces,6);
			faces[k].initialize(0);
			nodes[k].allocateMemory(nnodes);
			volume[k].allocateMemory(nnodes);
			volume[k].initialize(.0);
			ID[k].allocateMemory(nnodes);
			ID[k].initialize(0);

			edges_bdry[k].allocateMemory(nbedges,6);
			edges_bdry[k].initialize(0);
			volume_bdry[k].allocateMemory(nbnodes);
			volume_bdry[k].initialize(.0);
			ID_bdry[k].allocateMemory(nbnodes);
			ID_bdry[k].initialize(0);

			Cij_norm[k].allocateMemory(nedges);
			Cij_norm[k].initialize(.0);
			edge_length[k].allocateMemory(nedges);
			edge_length[k].initialize(.0);
			edge_versor[k].allocateMemory(nedges,3);
			edge_versor[k].initialize(.0);

			Dij[k].allocateMemory(nbedges,3);
			Dij[k].initialize(.0);
		}
	}

	void GeomData::deallocatePointers(){
		int ndom = getNumDomains();
		volume_global[0].freeMemory();
		EBE_1[0].freeMemory();
		EBE_2[0].freeMemory();
		delete[] volume_global; volume_global = 0;
		for (int k=0; k<ndom; k++){
			Cij[k].freeMemory();
			Cij_norm[k].freeMemory();
			edge_length[k].freeMemory();
			Dij[k].freeMemory();
			volume[k].freeMemory();
			edges[k].freeMemory();
			faces[k].freeMemory();
			nodes[k].freeMemory();
			ID_bdry[k].freeMemory();
			volume_bdry[k].freeMemory();
			edges_bdry[k].freeMemory();
			ID[k].freeMemory();
			edge_versor[k].freeMemory();
		}
		delete[] Cij; Cij = 0;
		delete[] numDomEdges; numDomEdges = 0;
		delete[] numDomBDRYEdges; numDomBDRYEdges = 0;
		delete[] Cij_norm; Cij_norm = 0;
		delete[] Dij; Dij = 0;
		delete[] volume; volume = 0;
		delete[] edges; edges = 0;
		delete[] faces; faces = 0;
		delete[] nodes; nodes = 0;
		delete[] ID; ID = 0;
		delete[] volume_bdry; volume_bdry = 0;
		delete[] edges_bdry; edges_bdry = 0;
		delete[] ID_bdry; ID_bdry = 0;
		delete[] edge_length; edge_length = 0;
		delete[] edge_versor; edge_versor = 0;
		delete[] EBE_1; EBE_1=0;
		delete[] EBE_2; EBE_2=0;
	}

	void GeomData::mappingNodesIds(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domainList = getDomainList();
		std::map<int,int> mapIDtoIndex, mapBdryIDtoIndex, mapIDtoIndex_global;
		std::set<int> setID, setBdryID;
		int i,j,k,id0,id1,id2;
		pEntity node, edge, face;

		i = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit))){
			mapIDtoIndex_global[EN_id(node)] = i++;
		}
		VIter_delete(vit);

		// STEP1: collect all node IDs from domain dom and node Ids from boudanry domain dom
		for (k=0; k<ndom; k++){
			int dom = domainList[k];
			FIter fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				int faceflag = getFaceFlag(face);
				if (faceflag==domainList[k]){
					// nodes domain
					for(i=0; i<3; i++){
						setID.insert( EN_id(face->get(0,i)) );
					}
					// boundary nodes
					for(i=0; i<3; i++){
						if (getVertexFlag(face->get(0,i))!=faceflag){
							setBdryID.insert(EN_id(face->get(0,i)));
						}
					}
				}
			}
			FIter_delete(fit);

			// get external boundary edges only
			i = 0;
			pEntity v1,v2;
			EIter eit = M_edgeIter(theMesh);
			while ( (edge=EIter_next(eit)) ){
				if (E_numFaces(edge)==1){
					v1 = edge->get(0,0);
					v2 = edge->get(0,1);
					EBE_2[0].setValue(i,0,mapIDtoIndex_global[EN_id(v1)]);
					EBE_2[0].setValue(i,1,mapIDtoIndex_global[EN_id(v2)]);
					EBE_2[0].setValue(i,2,GEN_tag(v1->getClassification()));
					EBE_2[0].setValue(i,3,GEN_tag(v2->getClassification()));
					i++;
				}
			}
			EIter_delete(eit);

			// STEP2: map IDs from step1 like this:
			//		id0 = 0;
			//		id1 = 1;
			//		id2 = 2;
			//		...
			//		idn = n;

			if (setID.size() > this->getNumNodesPerDomain(k)){
				char msg[256]; sprintf(msg,"NUmber of elements collected [%d] is greater than the max [%d]",(int)setID.size(),this->getNumNodesPerDomain(k));
				throw Exception(__LINE__,__FILE__,msg);
			}

			// domain
			i = 0;
			pVertex node;
			for (std::set<int>::iterator iter=setID.begin(); iter!= setID.end(); iter++){
				int id = *iter;
				mapIDtoIndex[id] = i;
				ID[k].setValue(i,id);
				node = theMesh->getVertex(id);
				volume[k].setValue(i,getVolume(node,dom) );
				nodes[k].setValue(i,mapIDtoIndex_global[id]);
				i++;
			}
			setID.clear();
			// boundary
			i = 0;
			for (std::set<int>::iterator iter=setBdryID.begin(); iter!= setBdryID.end(); iter++){
				int id = *iter;
				mapBdryIDtoIndex[id] = i;
				ID_bdry[k].setValue(i,id);
				node = theMesh->getVertex(id);
				volume_bdry[k].setValue(i,getVolume(node,dom) );
				i++;
			}
			setBdryID.clear();

			// STEP3: use mapped ID to initialize edge struct, where:
			//		edge_id0 = index_mapped
			//		edge_id1 = index_mapped
			i = 0;
			eit = M_edgeIter(theMesh);
			while ( (edge=EIter_next(eit)) ){
				if ( edgeBelongToDomain(edge,domainList[k]) ){
					v1 = edge->get(0,0);
					v2 = edge->get(0,1);
					id0 = EN_id( v1 );
					id1 = EN_id( v2 );
					edges[k].setValue(i,0,mapIDtoIndex[id0]); // edge stores the index and not node ID
					edges[k].setValue(i,1,mapIDtoIndex[id1]); // edge stores the index and not node ID
					edges[k].setValue(i,2,mapIDtoIndex_global[id0]); // edge stores the index and not node ID
					edges[k].setValue(i,3,mapIDtoIndex_global[id1]); // edge stores the index and not node ID
					edges[k].setValue(i,4,GEN_tag(v1->getClassification()));
					edges[k].setValue(i,5,GEN_tag(v2->getClassification()));
					i++;
				}
			}
			i = 0;
			eit = M_edgeIter(theMesh);
			while ( (edge=EIter_next(eit)) ){
				if ( edgeBelongToDomain(edge,domainList[k]) ){
					if (belongsToBoundary(edge)){
						id0 = EN_id( edge->get(0,0) );
						id1 = EN_id( edge->get(0,1) );
						edges_bdry[k].setValue(i,0,mapBdryIDtoIndex[id0]); // edge stores the index and not node ID
						edges_bdry[k].setValue(i,1,mapBdryIDtoIndex[id1]); // edge stores the index and not node ID
						edges_bdry[k].setValue(i,2,mapIDtoIndex[id0]); // edge stores the index and not node ID
						edges_bdry[k].setValue(i,3,mapIDtoIndex[id1]); // edge stores the index and not node ID
						edges_bdry[k].setValue(i,4,mapIDtoIndex_global[id0]); // edge stores the index and not node ID
						edges_bdry[k].setValue(i,5,mapIDtoIndex_global[id1]); // edge stores the index and not node ID
						i++;
					}
				}
			}

			i = 0;
			fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				int faceflag = getFaceFlag(face);
				if (faceflag==domainList[k]){
					for(j=0;j<3;j++){
						id0 = EN_id( face->get(0,j) );
						faces[k].setValue(i,j,mapIDtoIndex[id0]); // edge stores the index and not node ID
						faces[k].setValue(i,j+3,mapIDtoIndex_global[id0]); // edge stores the index and not node ID
					}
					i++;
				}
			}
			FIter_delete(fit);

			mapBdryIDtoIndex.clear();
			mapIDtoIndex.clear();
		}
		mapIDtoIndex_global.clear();
	}

	void GeomData::mappingNodesIds_Tmp(pMesh theMesh){
		int ndom = getNumDomains();
		const int* domainList = getDomainList();
		std::map<int,int> mapIDtoIndex, mapBdryIDtoIndex, mapIDtoIndex_global;
		std::set<int> setID, setBdryID;
		int i,j,k,id0,id1,id2;
		pEntity node, face;

		// global mapping
		i = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit))){
			mapIDtoIndex_global[EN_id(node)] = i++;
		}
		VIter_delete(vit);

		// calculate number of faces for new mesh (adapted mesh)
		numDomFaces_tmp = new int[ndom];
		for (k=0; k<ndom; k++){
			numDomFaces_tmp[k] = 0;
			FIter fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				if ( getFaceFlag(face)==domainList[k] ){
					numDomFaces_tmp[k]++;
				}
			}
			FIter_delete(fit);
		}

		// allocate faces_tmp
		faces_tmp = new Matrix<int>[ndom];
		for (k=0; k<ndom; k++){
			// number of faces per domain
			faces_tmp[k].allocateMemory(numDomFaces_tmp[k]);
			faces_tmp[k].initialize(0);

			// preparing for local mapping (by domain)
			int dom = domainList[k];
			FIter fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				int faceflag = getFaceFlag(face);
				if (faceflag==domainList[k]){
					// nodes domain
					for(i=0; i<3; i++){
						setID.insert( EN_id(face->get(0,i)) );
					}
				}
			}
			FIter_delete(fit);

			// creates mapping
			i = 0;
			pVertex node;
			for (std::set<int>::iterator iter=setID.begin(); iter!= setID.end(); iter++){
				int id = *iter;
				mapIDtoIndex[id] = i;
				i++;
			}
			setID.clear();

			// use mapping
			i = 0;
			fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				int faceflag = getFaceFlag(face);
				if (faceflag==domainList[k]){
					for(j=0;j<3;j++){
						id0 = EN_id( face->get(0,j) );
						faces_tmp[k].setValue(i,j,mapIDtoIndex[id0]); // edge stores the index and not node ID
						faces_tmp[k].setValue(i,j+3,mapIDtoIndex_global[id0]); // edge stores the index and not node ID
					}
					i++;
				}
			}
			FIter_delete(fit);
			mapIDtoIndex.clear();
		}
		mapIDtoIndex_global.clear();
	}
}
