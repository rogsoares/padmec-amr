#include "SimulatorParameters.h"
#include <sstream>

namespace PRS{
	
	SimulatorParameters::SimulatorParameters(){
	}
	
	SimulatorParameters::SimulatorParameters(pMesh mesh):theMesh(mesh){
		TimeStepCounter = 0;
		exportIter = 0;
		accSimTime = .0;
		stop_simulation = true;
		useHOApp = false;
		restart = false;
		vtk_step = 0;
		
		PVI_increment = 0.01;	// every 5% of total simulation time a new VTK file will be printed
		PVI_accumulated = .0;	// summation of all PVI_increments
		
		allowPrintingVTK = false;
		pctype = PCNONE;
		firstVTKupdate = true;
		EBFV1_pressureSolver_scheme = false;
		_doAdaptation = false;
	}
	
	SimulatorParameters::~SimulatorParameters(){
		MIter_RockProperties mit = mapRockProp.begin();
		for (; mit != mapRockProp.end(); mit++){
			double *K = mit->second->K;
			delete[] K; K = 0;
		}
		mapRockProp.clear();
		mapBC.clear();
	}
	
	void SimulatorParameters::checkPermeabilityTensor(){
		const double* K = getPermeability(3300);//BACALHO
		// Kxx*Kyy >= Kxy*Kyx
		if (K[0]*K[3] < K[1]*K[2]) throw Exception(__LINE__,__FILE__,"Permeability tensor must obey the following relation: K_xx*K_yy >= K_xy*K_yx\n");
		K_Isotropic = (K[1]*K[2] == .0)?true:false;
	}
	
	const double* SimulatorParameters::getPermeability(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()) return mit->second->K;
		cout << "Warning: no permeability tensor associated to domain " << dom << ".\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}
	
	void SimulatorParameters::getDomains(){
		for (MIter_RockProperties mit=mapRockProp.begin(); mit!=mapRockProp.end(); mit++)
			setOfDomains.insert(mit->first);
	}
	
	double SimulatorParameters::getBC_Value(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter != mapBC.end()) return mIter->second->val;
		cout << "Warning: getBC_Value() return null value\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		return 0;
	}
	
	bool SimulatorParameters::hasNodeWell(const int &flag){
		// if node has a prescribed saturation it's supposed that exist a well
		if ( mapSaturation.size() ){
			MapFlagIter MIter = mapSaturation.find(flag);
			return (MIter != mapSaturation.end());
		}
		else{
			isNodeFree(flag);
			return well;
		}
	}
	
	bool SimulatorParameters::isNodeFree(const int &flag){
		MapFlagIter mIter = mapBC.find(flag);
		if (mIter == mapBC.end()) return true;
		else{
			BdryConditionType *bct = mIter->second;
			return ( !bct->type.compare("dirichlet") )?false:true;
		}
	}
	
	bool SimulatorParameters::finishSimulation(){
		std:: cout << setprecision(2) << fixed;
		if (!stop_simulation || (getAccumulatedSimulationTime() >= getSimTime()) ){
			if (!P_pid()){
				std::cout << "#################################\n";
				std::cout << "Simulation " << (double)100.0*getAccumulatedSimulationTime()/getSimTime() << "% concluded.\n";
				std::cout << "End of Simulation\n";
				std::cout << "#################################\n";
			}
			return 1;
		}
		else{
			if (!P_pid()) std::cout << "Simulation " << (double)100.0*getAccumulatedSimulationTime()/getSimTime() << "% concluded.\n";
			return 0;
		}
	}
	
	// physical parameters
	double SimulatorParameters::getPorosity(const int &dom){
		MIter_RockProperties mit = mapRockProp.find(dom);
		if (mit != mapRockProp.end()) return mit->second->porosity;
		
		cout << "Warning: porosity() return null value.";
		cout << "Domain "<< dom << " not found.\n";
		cout << __FILE__ << "\t at line " << __LINE__ << endl;
		
		return 0;
	}
	
	// For a specific well, the flow rate for each node will be weighted by nodal volume.
	// Qt = (V1/Vt)Qt + (V2/Vt)Qt + ... + (Vn/Vt)Qt
	// where:
	// 		Qt 				- total flow rate (read from file)
	// 		Q1 = (V1/Vt)Qt 	- flow rate in node 1
	// 		Q2 = (V2/Vt)Qt 	- flow rate in node 2
	// 		Qn = (Vn/Vt)Qt 	- flow rate in node n
	void SimulatorParameters::weightWellFlowRateByVolume(pMesh theMesh, GeomData *pGCData){
		double Vi, Vt;
		pVertex node;
		map<int,set<int> >::iterator miter;	// for each well flag
		SIter_const siter;			// get all flagged node IDs for that well
		
		#ifdef _SEEKFORBUGS_
		if (!M_numVertices(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of vertices NULL");
		}
		if (!mapNodesOnWell.size()){
			throw Exception(__LINE__,__FILE__,"mapNodesOnWell size NULL");
		}
		#endif 
		
		// FOR EACH PRODUCTION WELL
		for (miter=mapNodesOnWell.begin(); miter!=mapNodesOnWell.end(); miter++){
			int flag = miter->first;
			// FOR EACH NODE ON PRODUCTION WELL
			for (siter = miter->second.begin(); siter!=miter->second.end();siter++){
				int id = *siter;
				Vt = .0;
				// to which domain node belongs get its volume
				for (SIter_const sit=setOfDomains.begin(); sit!=setOfDomains.end(); sit++)	{
					int dom = *sit;	// domains flag
					node = (mEntity*)theMesh->getVertex( id );
					
					#ifdef _SEEKFORBUGS_
					if (!node){
						char msg[256]; sprintf(msg,"NULL node with ID %d",id);
						throw Exception(__LINE__,__FILE__,msg);
					}
					#endif
					pGCData->getVolume(node,dom,Vi);
					Vt += Vi;
				}

				well_info = MWells[flag];
				well_info.wellVolume = Vt;
				MWells[flag] = well_info;
				cout << "flag = " << flag << endl;
				cout << "Vt   = " << Vt << endl;
				cout << "ID   = " << id << endl;
			}
		}
		//STOP();
		if (!MWells.size()){
			throw Exception(__LINE__,__FILE__,"MWells size = 0\n");
		}
	}
	
	// get node ids associated to wells
	void SimulatorParameters::getWells(pMesh theMesh, int dim){
		// if user has provided some well, count how many nodes are flagged for each one
		//cout << __LINE__ << "\t" << __FILE__ << endl;
		if ( MWells.size() ){
			pEntity node, edge;
			
			// associate a set container to each well flag
			// EX.: two injection and one production wells provided
			// flag		-	nodes IDS associated to that flag
			//  10		-	12,34,76,285,99,43,19,...
			//  11		-	1,3114,56,2285,399,643,189,...
			//  51		-	133,36,7,2855,929,463,149,...
			// initialize map container with a set container for each flag
			//cout << __LINE__ << "\t" << __FILE__ << endl;
			for (MWIter miter = MWells.begin();miter!=MWells.end();miter++){
				setNodesOnWell setNOW;
				mapNodesOnWell[miter->first] = setNOW;
			}
			//cout << __LINE__ << "\t" << __FILE__ << endl;
			int flag;
			//cout << __LINE__ << "\t" << __FILE__ << endl;
			// loop over nodes searching for wells
			VIter vit = M_vertexIter( theMesh );
			while ( (node = VIter_next(vit)) ){
				// get node's flag
				flag = (!node->getClassification())?0:GEN_tag( node->getClassification() );
				
				// check for INJECTION/PRODUCTION wells
				if ( this->isInjectionWell(flag) || this->isProductionWell(flag) ){
					MNOWIter = mapNodesOnWell.find(flag);
					// insert node ID for flag
					if ( MNOWIter!=mapNodesOnWell.end() )
						MNOWIter->second.insert( EN_id(node) );
				}
			}
			VIter_delete(vit);
			//cout << __LINE__ << "\t" << __FILE__ << endl;
			// loop edges searching for nodes flagged as wells
// 			EIter eit = M_edgeIter( theMesh );
// 			while ( (edge = EIter_next(eit)) ){
// 				if (!theMesh->getRefinementDepth(edge)){				// get only elements without children
// 					if (!edge->getClassification()){
// 						pEdge root = edge->root();
// 						flag = GEN_tag( root->getClassification());
// 					}
// 					else{
// 						flag = GEN_tag( edge->getClassification());
// 					}
// 					
// 					// check for INJECTION/PRODUCTION wells
// 					if ( ( flag>=10 && flag<=90) || (flag==2001) || (flag==2002)){
// 						MNOWIter = mapNodesOnWell.find(flag);
// 						if ( MNOWIter!=mapNodesOnWell.end() ){
// 							// insert node ID for flag
// 							MNOWIter->second.insert( EN_id(edge->get(0,0)) );
// 							MNOWIter->second.insert( EN_id(edge->get(0,1)) );
// 						}
// 					}
// 				}
// 			}
// 			EIter_delete(eit);
			//cout << __LINE__ << "\t" << __FILE__ << endl;
		}
		else{
			throw Exception(__LINE__,__FILE__,"MWells size = 0");
		}
	}
	
	bool SimulatorParameters::isInjectionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isInjection:false;
	}
	
	double SimulatorParameters::getInitialSaturation(pEntity node){
		int flag = GEN_tag(node->getClassification());
		MapFlagIter MIter = mapSaturation.find(flag);
		return (MIter == mapSaturation.end())?Sw_initial():(1.0 - Sor());
	}
	
	bool SimulatorParameters::isProductionWell(int flag) const{
		MWCIter mit = MWells.find( flag );
		return (mit!=MWells.end())?mit->second.isProduction:false;
	}
	
	double SimulatorParameters::getFlowrateValue(int flag) const{
		MWCIter mit = MWells.find( flag );
		return mit->second.flowRate;
	}
	
	double SimulatorParameters::getTotalInjectionFlowRate() const{
		double Q = 0.0;
		for (MWCIter mwiter=MWells.begin(); mwiter!=MWells.end(); mwiter++){
			//	printf("well flag %d. Flow rate: %f.  Volume: %f\n",mwiter->first,mwiter->second.flowRate,mwiter->second.wellVolume);
			if ( isInjectionWell( mwiter->first ) ) Q += mwiter->second.flowRate;
		}
		return Q;
	}
	
	double SimulatorParameters::getWellVolume(int well_flag) const{
		MWCIter mit = MWells.find( well_flag );
		return mit->second.wellVolume;
	}
	
	void SimulatorParameters::checkIfRankHasProductionWell(){
		VIter vit = M_vertexIter(theMesh);
		while (pVertex node = VIter_next(vit)){
			if ( isProductionWell(node) ){
				_rankHasProductionWell = true;
				break;
			}
		}
		VIter_delete(vit);
	}
	
	/*
	 * double GeomData::getReservoirVolume() returns the total reservoir volume.
	 * To compute the initial oil volume it's necessary to take account rock po-
	 * rosity that can vary from domain to domain. So, setInitialOilVolume() function
	 * make a loop over all control volumes from all partitions (if there are two
	 * or more)
	 */
	void SimulatorParameters::setInitialOilVolume(pMesh theMesh, GeomData *pGCData){
		pVertex node;
		double TPV = .0;	// Total Porous Volume
		for (SIter_const dom = setDomain_begin(); dom!=setDomain_end(); dom++){
			double TPV_tmp = .0;	// Total Porous Volume
			VIter vit = M_vertexIter(theMesh);
			while ( (node = VIter_next(vit)) ){
				if ( pGCData->nodeBelongToDomain(node,*dom) )
					TPV_tmp += pGCData->getVolume(node,*dom);///(pGCData->getNumRemoteCopies(node) + 1.0);
			}
			VIter_delete(vit);
			TPV_tmp *= getPorosity(*dom);
			TPV += TPV_tmp;
			//printf("TPV_tmp: %f TPV: %f\n",TPV_tmp,TPV);;
		}
		
		double Sw = 1.0 - Sw_initial();
		//printf("TPV: %f   TPV*Sw: %f\n",TPV,TPV*Sw); throw 1;
		_IOV = P_getSumDbl(TPV)*Sw;
	}
	
	
	/*
	 * This function is called every time a new time-step is calculated.
	 * We want to know how many steps there are within a fraction of PVI or how
	 * many step there are within every time allowPrintingVTK is set true.
	 * TSCountingList stores this values.
	 */
	void SimulatorParameters::correctTimeStep(double &timeStep){
		static int TScounter = 0; TScounter++;
		
		if (firstVTKupdate){
			updatePrintOutVTKFrequency();  
		}
		double timeFrequency = getPrintOutVTKFrequency();
		double accST = timeStep + getAccumulatedSimulationTime();
		if ( accST > timeFrequency ){
			timeStep = timeStep - (accST - timeFrequency);
			accST = timeFrequency;
			allowPrintingVTK = true;
			
			TSCountingList.push_back(TScounter); TScounter = 0;
			
			#ifdef __PRINTDEBUGGING__
			if (!P_pid()) std::cout << "VTK file output simulation time: " << accST << endl;
			#endif
		}
	}
	
	void SimulatorParameters::printOutVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, pFunc_PrintVTK printVTK){
		
		// 	#ifdef _SEEKFORBUGS_
		// 		allowPrintingVTK = true;
		// 	#endif
		// 
		// 	#ifdef CRUMPTON_EXAMPLE
		// 		allowPrintingVTK = true;
		// 	#endif
		
		// 	#ifdef TESTINGADAPTATION
		// 		allowPrintingVTK = true;
		// 	#endif
		
//		allowPrintingVTK=true;
//		allowPrintingVTK = true;
		if (allowPrintingVTK){
			static int theStep = getStepOutputFile();
			char fname[256];
			++theStep;
			sprintf(fname,"%s__%d-of-%d__step-%d.vtk",expofName.c_str(),P_pid(),P_size(),theStep);
			PetscPrintf(PETSC_COMM_WORLD,"VTK Output: step file #%d\n",theStep);
			printVTK(theMesh,pData1,pData2,pData3,fname);
			updatePrintOutVTKFrequency();
			allowPrintingVTK = false;
		}
	}
	
	void SimulatorParameters::updatePrintOutVTKFrequency(){
		firstVTKupdate = false;
		
		// a new VTK file will be printed at each 0.05 PVI (5% of total simulation time)
		PVI_accumulated += getPVIincrement();
		
		// when next vtk file must be print out
		vtk_time_frequency = PVI_accumulated*getSimTime();
	}
}
