#include "PhysicPropData.h"

Matrix<double> *pGrad_matrix;	// each pointer is a matrix nx3
Matrix<double> SwGrad;			// a unique matrix nnodesx3 for whole mesh
Matrix<double> *SwGrad_dom;		// for domain
Matrix<double> Sw;		// for domain

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

	void PhysicPropData::initialize(MeshData *pMData, SimulatorParameters *pSimPar, pMesh theMesh, bool isUpdate, GeomData* pGCData){
		int numNodes = M_numVertices(theMesh);
		int ndom = pSimPar->getNumDomains();
		const int* pNumNodesDom = pSimPar->getNumNodesDomain();		// let me know how many nodes belong to ech domain

		pGrad_matrix = new Matrix<double>[ndom];
		velocity = new Matrix<double>[ndom];
		SwGrad_dom = new Matrix<double>[ndom];
		for (int k=0; k<ndom; k++){
			int nrows = pNumNodesDom[k];
			int nedges = pGCData->getNumEdgesPerDomain(k);
			pGrad_matrix[k].allocateMemory(nrows,3);
			pGrad_matrix[k].initialize(.0);
			velocity[k].allocateMemory(nedges,6);
			velocity[k].initialize(.0);
			SwGrad_dom[k].allocateMemory(nrows,3);
			SwGrad_dom[k].initialize(.0);
		}
		int nnodes = M_numVertices(theMesh);
		SwGrad.allocateMemory(nnodes,3);
		SwGrad.initialize(.0);
		Sw.allocateMemory(nnodes);
		Sw.initialize(0);
		injectionWell.allocateMemory(nnodes);
		projectedSw_grad.allocateMemory(nnodes);
		nonvisc.allocateMemory(nnodes);

		Swr = pSimPar->Swr();				// Irreducible water saturation
		Sor = pSimPar->Sor();				// Residual oil saturation
		mi_w = pSimPar->waterViscosity();	// water viscosity
		mi_o = pSimPar->oilViscosity();		// oil viscosity

		// do not set initial saturation field for update (mesh adaptation)
		if (!isUpdate){
			setInitialSaturation(theMesh,pSimPar);
		}

		// TODO: bacalho. saturacoes muito pequenas <E-200
		int idx = 0;
		int k = 0;
		nfree = 0, nneumann = 0;
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			double Sw = this->getSaturation(idx);
			if ( Sw < 1e-8 ){
				Sw = .0;
			}
			setSaturation(node,Sw);
			this->setSaturation(idx,Sw);
			// todo: TIRAR ESSE BACALHO DAQUI!!!
			if ( pSimPar->isProductionWell(flag) ){
				nneumann++;
				injectionWell.setValue(k++,false);
			}
			if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
				nfree++;
				injectionWell.setValue(k++,false);
			}
			idx++;
		}
		VIter_delete(vit);


		pWellsFree_index = new int[nfree];
		pWellsNeumann_index = new int[nneumann];
		idx = 0;
		int free_idx = 0;
		int neumann_idx = 0;
		vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			if ( pSimPar->isProductionWell(flag) ){
				this->pWellsNeumann_index[neumann_idx] = idx;
				neumann_idx++;
			}
			if ( !pSimPar->isInjectionWell(flag) && !pSimPar->isProductionWell(flag) ){
				this->pWellsFree_index[free_idx] = idx;
				free_idx++;
			}
			idx++;
		}
		VIter_delete(vit);

		ksModel = pSimPar->ksModel();

#ifdef __SEEKFORBUGS__
		if (!numNodes || !numEdges || !ndom)
			throw Exception(__LINE__,__FILE__,"Null value! Exiting....\n");
#endif


	}

	void PhysicPropData::deallocateData(SimulatorParameters *pSimPar){
		int ndom = pSimPar->getNumDomains();
		for (int k=0; k<ndom; k++){
			pGrad_matrix[k].freeMemory();
			velocity[k].freeMemory();
			//SwGrad_matrix[k].freeMemory();
		}
//		SwGrad_matrix[ndom].freeMemory();
//		delete[] SwGrad_matrix; SwGrad_matrix = 0;
		delete[] pGrad_matrix; pGrad_matrix = 0;
		delete[] velocity; velocity = 0;
	}

	void PhysicPropData::setInitialSaturation(pMesh theMesh, SimulatorParameters *simPar){
		pVertex node;
		bool hasInjectionWell = false;
		double Sw;
		int idx = 0;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int flag = GEN_tag(node->getClassification());
			Sw = simPar->getInitialSaturation(node);
			if ( Sw > .0 ){
				printf("Injection well located in node %d Sw = %f flag: %d\n",EN_id(node),Sw,flag);
				hasInjectionWell = true;
			}
			setSaturation(node,Sw);
			this->setSaturation(idx,Sw);
			idx++;
		}
		VIter_delete(vit);

		// rank with injection well must say to all other ranks that it's OK, a injection well was being informed!
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
////				edge->get(0,0)->classify( edge->getClassification() );
////				edge->get(0,1)->classify( edge->getClassification() );
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
//	void PhysicPropData::storeSwField(pMesh theMesh){
//		pVertex node;
//		double Sw;
//		VIter vit = M_vertexIter(theMesh);
//		while ( (node=VIter_next(vit)) ){
//			Sw = getSaturation(node);
//			setSaturation_Old(node,Sw);
//		}
//		VIter_delete(vit);
//	}

	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * If breakthrough has been reached, step back to the beginning of implicit time-
	 * step and advance in time in old-fashion way. Retrieve saturation field.
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **/
//	void PhysicPropData::retrieveSwField(pMesh theMesh){
//		pVertex node;
//		double Sw;
//		VIter vit = M_vertexIter(theMesh);
//		while ( (node=VIter_next(vit)) ){
//			Sw = getSaturation_Old(node);
//			setSaturation(node,Sw);
//		}
//		VIter_delete(vit);
//	}

//	void PhysicPropData::setSw_min(pEntity node, double Sw){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		pNode->Sw_min = Sw;
//		setAttachedData_pointer(node,pNode);
//	}
//
//	void PhysicPropData::setSw_max(pEntity node, double Sw){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		pNode->Sw_max = Sw;
//		setAttachedData_pointer(node,pNode);
//	}
//
//	void PhysicPropData::setS_Limit(pEntity node, double S_Limit){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		pNode->S_Limit = S_Limit;
//		setAttachedData_pointer(node,pNode);
//	}
//
//	double PhysicPropData::getSw_min(pEntity node){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		return pNode->Sw_min;
//	}
//
//	double PhysicPropData::getSw_max(pEntity node){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		return pNode->Sw_max;
//	}

//	double PhysicPropData::getS_Limit(pEntity node){
//		NodePhysicalProperties* pNode = getAttachedData_pointer<NodePhysicalProperties>(node);
//		return pNode->S_Limit;
//	}

	double PhysicPropData::getTotalMobility(pEntity vertex){
		if (steady_state){
			return 1.0;
		}
		double Sw = getSaturation(vertex);
		double krw = get_ksw(Sw);
		double kro = get_kso(Sw);
		return krw/mi_w + kro/mi_o;
	}

	double PhysicPropData::getTotalMobility(double Sw){
		double krw = get_ksw(Sw);
		double kro = get_kso(Sw);
		return krw/mi_w + kro/mi_o;
	}

	double PhysicPropData::getFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
		return krw/( krw + (mi_w/mi_o)*kro );
	}

	double PhysicPropData::getOilFractionalFlux(const double &Sw){
		double krw = pow((Sw - Swr)/(1. - Swr - Sor),2);
		double kro = pow((1. - Sw - Swr)/(1. - Swr - Sor),2);
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
