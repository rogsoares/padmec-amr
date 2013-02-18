#include "GeomData.h"
//
//Matrix<double> Cij_matrix;
//Matrix<double> Cij_norm_matrix;
//Matrix<double> Dij_matrix;
//Matrix<double> volume_matrix;


namespace PRS
{

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
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(edge);
		MapDataArray::const_iterator MDA_citer = pCoeffnt->Cij.find(dom);
		return ( MDA_citer != pCoeffnt->Cij.end() );
	}

	bool GeomData::nodeBelongToDomain(pEntity node, const int &dom){
		Coefficients* pCoeffnt = getAttachedData_pointer<Coefficients>(node);
		MapData::const_iterator MDA_citer = pCoeffnt->volume.find(dom);
		return ( MDA_citer != pCoeffnt->volume.end() );
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
}

