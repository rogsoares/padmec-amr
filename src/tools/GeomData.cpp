#include "GeomData.h"

namespace PRS{
	GeomData::GeomData(){
		setMeshDataId("GC_Data");

		/*
		 * Set reservoir height as unitary for 3-D domains
		 */
		reservoirHeight = 1.0;
	}

	GeomData::~GeomData(){
	}

	void GeomData::getCij(pEntity edge, const int &dom, double *Cij){
		DataArray cij(3);
		getCij(edge,dom,cij);
		Cij[0] = cij[0];
		Cij[1] = cij[1];
		Cij[2] = cij[2];
	}

	void GeomData::getCij(pEntity edge, const int &dom, DataArray &Cij){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		MapDataArray::const_iterator MDA_citer = pCoeffnt->Cij.find(dom);
		if ( MDA_citer != pCoeffnt->Cij.end() )
			Cij = MDA_citer->second;
		else
			std::fill(Cij.begin(),Cij.end(),0.0);
	}

	//	Dij can belong to one or two domains. If between two domains it must point to
	//	outside of the domain 'dom'. From pre-processor file, each Dij is followed by two
	//	integers flags corresponding to domains associated to it. Dij ALWAYS points to outside
	//	of the domain corresponding to the FIRST flag. Every time user calls getDij(),
	//	dom argument will be compared to FIRST flag to give Dij the right orientation.
	bool GeomData::getDij(pEntity edge, const int &dom, DataArray &Dij){
		// says if Dij belong to dom
		bool Dij_exist_on_domain = false;
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		//	printf("doms: %d %d  flag: %d\n",pCoeffnt->dom1,pCoeffnt->dom2,GEN_tag( edge->getClassification() ));
		if (dom == pCoeffnt->dom1 || dom == pCoeffnt->dom2){
			// compare first flag to dom
			double signal = (dom == pCoeffnt->dom1)?-1.0:1.0;
			// get Dij
			Dij = pCoeffnt->Dij;
			// give Dij right orientation
			for (int i=0; i<(int)Dij.size(); i++) Dij[i] = -signal*Dij[i];
			Dij_exist_on_domain = true;
		}
		pCoeffnt = 0;
		return Dij_exist_on_domain;
	}

	bool GeomData::getDij(pEntity face, int dom, double *Dij){
		int dom1=0, dom2=0;
		EN_getDataInt(face,MD_lookupMeshDataId("dom1"),&dom1);
		EN_getDataInt(face,MD_lookupMeshDataId("dom2"),&dom2);
		bool Dij_exist_on_domain = false;

		if (dom == dom1 && dom == dom2)
			return false;

		//printf("dom1: %d   dom2: %d  dom: %d\n",dom1,dom2,dom);
		if (dom == dom1 || dom == dom2){
			double signal = (dom == dom1)?-1.0:1.0;						// check Dij orientation
			EN_getDataDbl(face,MD_lookupMeshDataId("Dij_x"),&Dij[0]);	// Dij x coordenate
			EN_getDataDbl(face,MD_lookupMeshDataId("Dij_y"),&Dij[1]);	// Dij y coordenate
			EN_getDataDbl(face,MD_lookupMeshDataId("Dij_z"),&Dij[2]);	// Dij z coordenate
			for (int i=0; i<3; i++) Dij[i] = -signal*Dij[i];			// give Dij right orientation
			Dij_exist_on_domain = true;
		}
		return Dij_exist_on_domain;
	}

	void GeomData::setDij(pEntity face, int dom1, int dom2, double *Dij){
		EN_attachDataInt(face,MD_lookupMeshDataId("dom1"),dom1);
		EN_attachDataInt(face,MD_lookupMeshDataId("dom2"),dom2);
		EN_attachDataDbl(face,MD_lookupMeshDataId("Dij_x"),Dij[0]);	// Dij x coordenate
		EN_attachDataDbl(face,MD_lookupMeshDataId("Dij_y"),Dij[1]);	// Dij y coordenate
		EN_attachDataDbl(face,MD_lookupMeshDataId("Dij_z"),Dij[2]);	// Dij z coordenate
	}

	void GeomData::setDij(pEntity edge, const int &dom1, const int &dom2, const DataArray &Dij){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->Dij = Dij;
		pCoeffnt->dom1 = dom1;
		pCoeffnt->dom2 = dom2;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	double GeomData::getVolume(pEntity node, const int &dom){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		MData_Iter iter = pCoeffnt->volume.find(dom);
		return ( iter != pCoeffnt->volume.end() )?pCoeffnt->volume[dom]:0.0;
	}

	double GeomData::getWeightedVolume(pEntity node){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		return pCoeffnt->weightedvolume;
	}

	void GeomData::setCij(pEntity edge, const int &dom, DataArray Cij){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->Cij[dom] = Cij;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	void GeomData::setCij_norm(pEntity edge, const int &dom, double norm){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->Cij_norm = norm;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	double GeomData::getCij_norm(pEntity edge, const int &dom){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		return pCoeffnt->Cij_norm;
	}

	void GeomData::setVolume(pEntity node, const int &dom, const  double &v){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		pCoeffnt->volume[dom] = v;
		setAttachedData_pointer(node,pCoeffnt);
	}

	int GeomData::getDomainFlag(pEntity node){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		std::map<int,double>:: iterator iter = pCoeffnt->volume.begin();
		return iter->first;
	}

	void GeomData::setWeightedVolume(pEntity node, const  double &v){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		pCoeffnt->weightedvolume = v;
		setAttachedData_pointer(node,pCoeffnt);
	}

	void GeomData::getVolume(pEntity node, int dom, double &vol){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		MData_Iter iter = pCoeffnt->volume.find(dom);
		if ( iter != pCoeffnt->volume.end() ) vol = iter->second;
	}

	bool GeomData::edgeBelongToDomain(pEntity edge, const int &dom){
		//		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		//		MapDataArray::const_iterator MDA_citer = pCoeffnt->Cij.find(dom);
		//		return ( MDA_citer != pCoeffnt->Cij.end() );
		int belong = 0;
		char dom_string[256]; sprintf(dom_string,"dom_str_%d",dom);
		EN_getDataInt(edge,MD_lookupMeshDataId( dom_string ),&belong);
		return belong;
	}

	bool GeomData::nodeBelongToDomain(pEntity node, const int &dom){
		//		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		//		MapData::const_iterator MDA_citer = pCoeffnt->volume.find(dom);
		//		return ( MDA_citer != pCoeffnt->volume.end() );
		int belong = 0;
		char dom_string[256]; sprintf(dom_string,"dom_str_%d",dom);
		EN_getDataInt(node,MD_lookupMeshDataId( dom_string ),&belong);
		return belong;
	}

	void GeomData::setEdgeLength(pEntity edge, double length){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->length = length;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	double GeomData::getEdgeLength(pEntity edge){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		return pCoeffnt->length;
	}

	void GeomData::setEdgeVector(pEntity edge, dblarray vec){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->edgeVector = vec;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	void GeomData::getEdgeVector(pEntity edge, dblarray& vec){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		vec = pCoeffnt->edgeVector;
	}

	void GeomData::setEdgeVec_Unitary(pEntity edge, dblarray vec){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->edgeVectorUnitary = vec;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	void GeomData::getEdgeVec_Unitary(pEntity edge, dblarray& vec){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		vec = pCoeffnt->edgeVectorUnitary;
	}

	void GeomData::setTotalReservoirVolume(double V_local){
		reservoirVolume = P_getSumDbl(V_local);
	}

	double GeomData::getReservoirVolume() const{
		return reservoirVolume;
	}

	void GeomData::setNumRemoteCopies(pEntity node, int nrc){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		pCoeffnt->numRemoteCopies = nrc;
		setAttachedData_pointer(node,pCoeffnt);
	}

	int GeomData::getFlag(pEntity e){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(e);
		return pCoeffnt->flag;
	}

	void GeomData::setFlag(pEntity e, int flag){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(e);
		pCoeffnt->flag = flag;
		setAttachedData_pointer(e,pCoeffnt);
	}

	void GeomData::set_belongsToBoundary(pEdge edge, bool k){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		pCoeffnt->isBoundary = k;
		setAttachedData_pointer(edge,pCoeffnt);
	}

	bool GeomData::belongsToBoundary(pEdge edge){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		return pCoeffnt->isBoundary;
	}

	void GeomData::setNumRC(pEntity ent, int dom, int nrc){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(ent);
		pCoeffnt->numRC[dom] = nrc;
		setAttachedData_pointer(ent,pCoeffnt);
	}

	void GeomData::setSumIJ(pEntity node, DataArray sumIJ){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		pCoeffnt->sumIJ = sumIJ;
		setAttachedData_pointer(node,pCoeffnt);
	}

	void GeomData::getSumIJ(pEntity node, DataArray& sumIJ){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		sumIJ = pCoeffnt->sumIJ;
	}

	void GeomData::calculateNumNodes(pMesh theMesh, int ndom, int* domList){
		numNodesPerDomain = new int[ndom];
		pEntity face;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domList[i]){
					for (int j = 0; j<3; j++){
						nodeSet.insert(face->get(0,j));
					}
				}
			}
			FIter_delete(fit);
			numNodesPerDomain[i] = (int)nodeSet.size();
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumBdryNodes(pMesh theMesh, int ndom, int* domList){
		numBdryNodesPerDomain = new int[ndom];
		pEntity face;
		std::set<pEntity> nodeSet;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				int faceflag = getFaceFlag(face);
				if (faceflag==domList[i]){
					for (int j = 0; j<3; j++){
						if (getVertexFlag(face->get(0,j))!=faceflag){
							nodeSet.insert(face->get(0,j));
						}
					}
				}
			}
			FIter_delete(fit);
			numBdryNodesPerDomain[i] = (int)nodeSet.size();
			nodeSet.clear();
		}
	}

	void GeomData::calculateNumEdges(pMesh theMesh, int ndom, int* domList){
		numDomEdges = new int[ndom];
		this->domList = new int[ndom];
		for (int i=0; i<ndom; i++){
			this->domList[i] = domList[i];
		}
		pEntity face;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domList[i]){
					for (int j = 0; j<3; j++){
						edgeList.insert(face->get(1,j));
					}
				}
			}
			FIter_delete(fit);
			numDomEdges[i] = (int)edgeList.size();
			edgeList.clear();
		}
	}

	void GeomData::calculateNumBDRYEdges(pMesh theMesh, int ndom, int* domList){
		numDomBDRYEdges = new int[ndom];
		pEntity edge, face;
		std::set<pEntity> edgeList;
		for (int i=0; i<ndom; i++){
			FIter fit = M_faceIter(theMesh);
			while ( (face = FIter_next(fit)) ){
				if (getFaceFlag(face)==domList[i]){
					for (int j = 0; j<3; j++){
						edge = (pEdge)face->get(1,j);
						if ( getEdgeFlag(edge)!= domList[i]){
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

	void GeomData::allocatePointers(pMesh theMesh, int ndom){
		Cij = new Matrix<double>[ndom];
		Cij_norm = new Matrix<double>[ndom];
		Dij = new Matrix<double>[ndom];
		volume = new Matrix<double>[ndom];
		edges = new Matrix<int>[ndom];
		nodes = new Matrix<int>[ndom];
		ID = new Matrix<int>[ndom];
		volume_bdry = new Matrix<double>[ndom];
		edges_bdry = new Matrix<int>[ndom];
		ID_bdry = new Matrix<int>[ndom];
		edge_length = new Matrix<double>[ndom];
		edge_versor = new Matrix<double>[ndom];
		int nnodes_global = M_numVertices(theMesh);
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

			Cij[k].allocateMemory(nedges,3);
			Cij[k].initialize(.0);
			edges[k].allocateMemory(nedges,6);
			edges[k].initialize(0);
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

	void GeomData::deallocatePointers(int ndom){
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

	void GeomData::getCij(int dom, int row, double* cij){
		cij[0] = Cij[dom].getValue(row,0);
		cij[1] = Cij[dom].getValue(row,1);
		cij[2] = Cij[dom].getValue(row,2);
	}
	void GeomData::setCij(int dom, int row, double* cij){
		Cij[dom].setValue(row,0,cij[0]);
		Cij[dom].setValue(row,1,cij[1]);
		Cij[dom].setValue(row,2,cij[2]);
	}

	void GeomData::getDij(int dom, int row, double* dij){
		dij[0] = Dij[dom].getValue(row,0);
		dij[1] = Dij[dom].getValue(row,1);
		dij[2] = Dij[dom].getValue(row,2);
	}
	void GeomData::setDij(int dom, int row, double* dij){
		Dij[dom].setValue(row,0,dij[0]);
		Dij[dom].setValue(row,1,dij[1]);
		Dij[dom].setValue(row,2,dij[2]);
	}

	void GeomData::cleanData(pMesh theMesh){
		EIter eit = M_edgeIter(theMesh);
		while ( pEdge edge = EIter_next(eit) ){

		}
		EIter_delete(eit);
	}

	void GeomData::mappingNodesIds(pMesh theMesh, int ndom, int* domlist){
		std::map<int,int> mapIDtoIndex, mapBdryIDtoIndex, mapIDtoIndex_global;
		std::set<int> setID, setBdryID;
		int i,k,id0,id1;
		pEntity node, edge, face;

		i = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit))){
			mapIDtoIndex_global[EN_id(node)] = i++;
		}
		VIter_delete(vit);


		// STEP1: collect all node IDs from domain dom and node Ids from boudanry domain dom
		for (k=0; k<ndom; k++){
			int dom = domlist[k];
			FIter fit = M_faceIter(theMesh);
			while ( (face=FIter_next(fit) )){
				int faceflag = getFaceFlag(face);
				if (faceflag==domlist[k]){
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
				if ( edgeBelongToDomain(edge,domlist[k]) ){
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
				if ( edgeBelongToDomain(edge,domlist[k]) ){
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
			mapBdryIDtoIndex.clear();
			mapIDtoIndex.clear();
		}
		mapIDtoIndex_global.clear();
	}

	void GeomData::calculateEdgeProperties(pMesh theMesh,int ndom, int* domList){
		double c1[3], c2[3], vec[3];
		double elength, Lij[3];
		double Icoord[3], Jcoord[3];
		double delta_x = 1e30;
		int i,j;
		int dim = theMesh->getDim();
		pEntity edge;
		for (int dom = 0; dom<ndom; dom++){
			int row = 0;
			int dom_flag = domList[dom];
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
		cout << "smalest = " << delta_x << endl;

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
