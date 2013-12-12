#include "PhysicPropData.h"

Matrix<double> *pGrad_matrix;
Matrix<double> SwGrad_matrix;
Matrix<double> nonvisc_matrix;

namespace PRS{

	PhysicPropData::PhysicPropData(){
		setMeshDataId("PP_Data");

		steady_state = false;

		pGetGradArray = new GetPFuncGrad[2];
		pGetGradArray[0] = get_pw_Grad;
		pGetGradArray[1] = get_Sw_Grad;
	}

	PhysicPropData::~PhysicPropData(){
	}

	void PhysicPropData::initialize(MeshData *pMData, SimulatorParameters *pSimPar, pMesh theMesh, bool isUpdate){
		Swr = pSimPar->Swr();				// Irreducible water saturation
		Sor = pSimPar->Sor();				// Residual oil saturation
		mi_w = pSimPar->waterViscosity();	// water viscosity
		mi_o = pSimPar->oilViscosity();		// oil viscosity

		// do not set initial saturation field for update (mesh adaptation)
		if (!isUpdate){
			setInitialSaturation(theMesh,pSimPar);
		}

		// TODO: bacalho. saturacoes muito pequenas <E-200
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			double Sw = this->getSaturation(node);
			if ( Sw < 1e-8 ){
				Sw = .0;
			}
			setSaturation(node,Sw);
		}
		VIter_delete(vit);

		ksModel = pSimPar->ksModel();
		//int dim = theMesh->getDim();
		//int numEdges = M_numEdges(theMesh);
		int numNodes = M_numVertices(theMesh);
		int ndom = pSimPar->getNumDomains();
		const int* pNumNodesDom = pSimPar->getNumNodesDomain();		// let me know how many nodes belong to ech domain

#ifdef __SEEKFORBUGS__
		if (!numNodes || !numEdges || !ndom)
			throw Exception(__LINE__,__FILE__,"Null value! Exiting....\n");
#endif

		pGrad_matrix = new Matrix<double>[ndom];
		for (int k=0; k<ndom; k++){
			int nrows = pNumNodesDom[k];
			pGrad_matrix[k].allocateMemory(nrows,3);
			pGrad_matrix[k].initialize(.0);
			//pGrad_matrix[k].printNumRowsCols();
		}
		SwGrad_matrix.allocateMemory(numNodes,3);
		SwGrad_matrix.initialize(.0);
	}

	void PhysicPropData::deallocateData(SimulatorParameters *pSimPar){
		int ndom = pSimPar->getNumDomains();
		for (int k=0; k<ndom; k++){
			pGrad_matrix[k].freeMemory();
		}
		delete[] pGrad_matrix; pGrad_matrix = 0;
		SwGrad_matrix.freeMemory();
	}

	void PhysicPropData::setInitialVelocity(pMesh theMesh, SimulatorParameters *pSimPar){
		pEdge edge;
		for (SIter sit=pSimPar->setDomain_begin(); sit!=pSimPar->setDomain_end(); sit++)	{
			const int dom = *sit;
			EIter eit = M_edgeIter(theMesh);
			while ( (edge = EIter_next(eit)) ){
				dblarray v(3,.0);
				setVelocity_new(edge,dom,v);
				setVelocity_old(edge,dom,v);
			}
			EIter_delete(eit);
		}
	}

	void PhysicPropData::setVelocity_new(pEdge edge, const int &dom, dblarray v){
		setVelocity(edge,dom,true,v);
	}

	void PhysicPropData::setVelocity_old(pEdge edge, const int &dom, dblarray v){
		setVelocity(edge,dom,false,v);
	}

	void PhysicPropData::setVelocity(pEdge edge, const int &dom, bool time, dblarray v){
		EdgePhysicalProperties* pEdgeCoeff = getAttachedData_pointer<EdgePhysicalProperties>(edge);
		if (pEdgeCoeff){
			if ( time )
				pEdgeCoeff->v_new[dom] = v;
			else
				pEdgeCoeff->v_old[dom] = v;
			setAttachedData_pointer(edge,pEdgeCoeff);
		}
		else
			throw Exception(__LINE__,__FILE__,"Null pointer.\n");
	}

	void PhysicPropData::getVelocity_new(pEdge edge, const int &dom, dblarray &v){
		getVelocity(edge,dom,true,v);
	}

	void PhysicPropData::getVelocity_old(pEdge edge, const int &dom, dblarray &v){
		getVelocity(edge,dom,false,v);
	}

	void PhysicPropData::getVelocity(pEdge edge, const int &dom, bool time, dblarray &v){
		EdgePhysicalProperties *pEdgeCoeff = getAttachedData_pointer<EdgePhysicalProperties>(edge);
		if ( time )
			v = pEdgeCoeff->v_new[dom];
		else
			v = pEdgeCoeff->v_old[dom];
	}

	void PhysicPropData::setInitialSaturation(pMesh theMesh, SimulatorParameters *simPar){

		/*
		 * IF MESH IS 2-D
		 * */
		pVertex node;
		bool hasInjectionWell = false;
		double Sw;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			Sw = simPar->getInitialSaturation(node);
			if ( Sw > .0 ){
				printf("Injection well located in node %d Sw = %f flag: %d\n",EN_id(node),Sw,flag);
				hasInjectionWell = true;
			}
			setSaturation(node,Sw);
		}
		VIter_delete(vit);

		vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			Sw = simPar->getInitialSaturation(node);
		}
		VIter_delete(vit);

		// rank with injection well must say to all other ranks that it's OK, a injection well
		// was being informed!
		int well = (hasInjectionWell)?1:0;
		well = P_getMaxInt(well);

		if (!well){
			//throw Exception(__LINE__,__FILE__,"Injection wells are missing!");
		}

		/*
		 * IF MESH IS 3-D
		 * */
		// search for flagged edges
		// =========================================================================
//		pEntity edge;
//		EIter eit = M_edgeIter( theMesh );
//		while ( (edge = EIter_next(eit)) ){
//			int flag = GEN_tag(edge->getClassification());
//			if (simPar->isInjectionWell(flag)){
//				setSaturation(edge->get(0,0),1.0);
//				setSaturation(edge->get(0,1),1.0);
//				edge->get(0,0)->classify( edge->getClassification() );
//				edge->get(0,1)->classify( edge->getClassification() );
//			}
//		}
//		EIter_delete(eit);
//
//		// search for flagged faces (on boundary only)
//		// =========================================================================
//		pEntity face;
//		FIter fit = M_faceIter( theMesh );
//		while ( (face = FIter_next(fit)) ){
//			int flag = GEN_tag(face->getClassification());
//			if (simPar->isInjectionWell(flag)){
//				for (int i=0; i<3; i++){
//					setSaturation(face->get(0,i),1.0);
//					face->get(0,i)->classify( face->getClassification() );
//				}
//			}
//		}
//		FIter_delete(fit);
		//throw 1;
	}


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * In modified IMPES implementation, it's desired to keep saturation field stored
	 * at the begin of every new implicit time-step. When breakthrough is reached it
	 * will be necessary to step-back to this point and advance in old-fashion way
	 * until breakthrough be reached again.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/
	void PhysicPropData::storeSwField(pMesh theMesh){
		pVertex node;
		double Sw;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit)) ){
			Sw = getSaturation(node);
			setSaturation_Old(node,Sw);
		}
		VIter_delete(vit);
	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * If breakthrough has been reached, step back to the beginning of implicit time-
	 * step and advance in time in old-fashion way. Retrieve saturation field.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/
	void PhysicPropData::retrieveSwField(pMesh theMesh){
		pVertex node;
		double Sw;
		VIter vit = M_vertexIter(theMesh);
		while ( (node=VIter_next(vit)) ){
			Sw = getSaturation_Old(node);
			setSaturation(node,Sw);
		}
		VIter_delete(vit);
	}

	void PhysicPropData::setSw_min(pEntity node, double Sw){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		pNode->Sw_min = Sw;
		setAttachedData_pointer(node,pNode);
	}

	void PhysicPropData::setSw_max(pEntity node, double Sw){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		pNode->Sw_max = Sw;
		setAttachedData_pointer(node,pNode);
	}

	void PhysicPropData::setS_Limit(pEntity node, double S_Limit){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		pNode->S_Limit = S_Limit;
		setAttachedData_pointer(node,pNode);
	}

	double PhysicPropData::getSw_min(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		return pNode->Sw_min;
	}

	double PhysicPropData::getSw_max(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		return pNode->Sw_max;
	}

	double PhysicPropData::getS_Limit(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		return pNode->S_Limit;
	}

/**
 * check if node gradient was projected (boundary nodes only!)
 */
	bool PhysicPropData::isProjected(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		return pNode->projected;
	}

	void PhysicPropData::setAsProjected(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		pNode->projected = true;
	}

	void PhysicPropData::setAsNOTProjected(pEntity node){
		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
		pNode->projected = false;
	}

	bool PhysicPropData::isVelProjected(pEntity edge){
		EdgePhysicalProperties* pEdgeCoeff = getAttachedData_pointer<EdgePhysicalProperties>(edge);
		return pEdgeCoeff->projected;
	}

	void PhysicPropData::setVelAsProjected(pEntity edge){
		EdgePhysicalProperties* pEdgeCoeff = getAttachedData_pointer<EdgePhysicalProperties>(edge);
		pEdgeCoeff->projected = true;
	}

	void PhysicPropData::setVelAsNOTProjected(pEntity edge){
		EdgePhysicalProperties* pEdgeCoeff = getAttachedData_pointer<EdgePhysicalProperties>(edge);
		pEdgeCoeff->projected = false;
	}

	double PhysicPropData::getTotalMobility(pEntity vertex){
		if (steady_state){
			return 1.0;
		}
		double Sw = getSaturation(vertex);
		double krw = get_ksw(Sw);
		double kro = get_kso(Sw);
		return krw/mi_w + kro/mi_o;
	}

	double PhysicPropData::getFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);//get_ksw(Sw);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);//get_kso(Sw);
		return krw/( krw + (mi_w/mi_o)*kro );
	}

	double PhysicPropData::getOilFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);//get_ksw(Sw);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);//get_kso(Sw);
		return (kro/mi_o)/( krw/mi_w + kro/mi_o );
	}

	double PhysicPropData::get_ksw(const double &Sw){
		switch(ksModel){
		case 1:
			return Sw;
		case 2:
			return pow(Sw,2);
		case 3:
			return pow((Sw - Swr)/(1. - Swr - Sor),2);
		default:
			throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
		}
	}

	double PhysicPropData::get_kso(const double &Sw){
		switch(ksModel){
		case 1:
			return 1. - Sw;
		case 2:
			return pow(1. - Sw,2);
		case 3:
			return pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
		default:
			throw Exception(__LINE__,__FILE__,"Unknown model for relative permeability.\n");
		}
	}
}
