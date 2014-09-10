/*
 * LoadSimulationParameters.cpp
 *
 *  Created on: 19/01/2011
 *      Author: rogerio
 */

#include "SimulatorParameters.h"
#include "EBFV1__pre-processors.h"
#include <sstream>

namespace PRS{

	void SimulatorParameters::inputDataFromFiles(GeomData *pGCData){
		getParametersDirectoryPath();
		loadNumericParameters(pGCData);
		load_preprocessorData(pGCData);
		loadPhysicalParameters(pGCData->getMeshDim());
		loadSolverParameters();
		loadHighOrderParameter();
		loadMeshAdaptationParameters();
	}

	void SimulatorParameters::getParametersDirectoryPath(){
		if (!P_pid()) std::cout << "getting parameters directory path... ";

		// read numeric parameters
		ifstream fid;
		fid.open("simulation-parameters/filenameParameters.dat");
		if ( !fid.is_open() ){
			throw Exception(__LINE__,__FILE__,"filenameParameters.dat could not be opened or it does not exist\n");
		}

		setPositionToRead(fid,"path:");
		fid >> parametersPath;

		// check if numeric.dat physical.dat slope_limiters.dat and solve.dat files are presented
		string checkForNumeric(parametersPath + "/numeric.dat");
		string checkForPhysical(parametersPath + "physical.dat");
		string checkForSlope_limiters(parametersPath + "/slope_limiters.dat");
		string checkForSolve(parametersPath + "/solve.dat");

		ifstream fid1, fid2, fid3, fid4;
		fid1.open(checkForNumeric.c_str());
		fid2.open(checkForPhysical.c_str());
		fid3.open(checkForSlope_limiters.c_str());
		fid4.open(checkForSolve.c_str());

		string msg1("could not be opened or it does not exist\n");
		char msg[256];

		if ( !fid1.is_open() ) sprintf(msg,"numeric.dat %s",msg1.c_str());
		if ( !fid2.is_open() ) sprintf(msg,"physical.dat %s",msg1.c_str());
		if ( !fid3.is_open() ) sprintf(msg,"slope-limiters.dat %s",msg1.c_str());
		if ( !fid4.is_open() ) sprintf(msg,"solve.dat %s",msg1.c_str());

		string failMsg(msg);

		if (!failMsg.length()){
			throw Exception(__LINE__,__FILE__,msg);
		}
		else{
			fid1.close();
			fid2.close();
			fid3.close();
			fid4.close();
		}
		if (!P_pid()) std::cout << "OK.\n";
	}

	void SimulatorParameters::initialize(GeomData *pGCData, pMesh theMesh){
#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: SimulatorParameters::initialize\tIN\n";
#endif
		if (!M_numFaces(theMesh)){
			throw Exception(__LINE__,__FILE__,"Number of mesh elements null!");
		}

		if (!this->useRestart()){
			setStepOutputFile(0);
		}

		getWells(theMesh,pGCData->getMeshDim());
		weightWellFlowRateByVolume(theMesh,pGCData);
		checkIfRankHasProductionWell();
		setInitialOilVolume(theMesh,pGCData);
		setLocalNodeIDNumbering(theMesh);
#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: SimulatorParameters::initialize\tOUT\n";
#endif
	}

	void SimulatorParameters::deallocateData(){
		mapNodesOnWell.clear();
		_mapNodesDomain.clear();
		delete[] numNodesDom; numNodesDom = 0;
	}

	void SimulatorParameters::load_preprocessorData(void *pData){
		switch ( getEllipticSolver() ){
		case 1:{
			int ndom;
			cout << prepFilename().c_str() << endl;

			ifstream fid;
			fid.open(prepFilename().c_str());
			if ( !fid.is_open() ){
				char msg[256]; sprintf(msg,"File: %s \ncould not be opened or it doesn't exist. Verify if path is correct.\n",prepFilename().c_str());
				throw Exception(__LINE__,__FILE__,msg);
			}
			fid.close();

			// load mesh using FMDB
			M_load(theMesh,prepFilename().c_str());
			
			if (theMesh->getDim()==2){
				EBFV1_preprocessor_2D(theMesh,pData,ndom);
			}
			else{
				EBFV1_preprocessor_3D(theMesh,pData,ndom);
			}

			if (ndom != getNumDomains()){
				char msg[256]; sprintf(msg,"Number of domains do not match."
						" Domains: %d (physical.dat) %d (pre-processor)\n",getNumDomains(),ndom);
			}
			break;
		}
		default:
			throw Exception(__LINE__,__FILE__,"Only case 1 for elliptic solver is implemented.\n");
		}
	}

	void SimulatorParameters::loadNumericParameters(GeomData *pGCData){
		if (!P_pid()) std::cout << "loading numeric parameters... ";
		string str;
		ifstream fid;

		// read numeric parameters
		string numeric(parametersPath + "/numeric.dat");
		fid.open(numeric.c_str());
		if ( !fid.is_open() ){
			char msg[256]; sprintf(msg,"File: %s \ncould not be opened or it doesn't exist\n",numeric.c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}


		// start reading file
		// -------------------------------------------------------------------------
		setPositionToRead(fid,"type method ID right above:");
		fid >> ellipticSolver;

		setPositionToRead(fid,"type method ID right above:");
		fid >> hyperbolicSolver;

		setPositionToRead(fid,"Use high order approximation:");
		fid >> str;
		if (!str.compare("yes")) setUseHOApproximation();

		setPositionToRead(fid,"C.F.L (courant-friedilich-lendroff number):");
		fid >> _CFL;

		setPositionToRead(fid,"Discretization time on production well:");
		fid >> _timeStepWell;

		if (!str.compare("yes") && _CFL > .5)
			throw Exception(__LINE__,__FILE__,"High order approximation cannot use CFL condition number greater than 0.5!\n");

		setPositionToRead(fid,"mesh filename:");
		fid >> prepfName;

		setPositionToRead(fid,"export filename (to VTK format without extension .vtk):");
		fid >> expofName;

		string strYes;
		setPositionToRead(fid,"Start simulation from <restart> file?");
		getline(fid,strYes);

		// RESTART: check if simulation starts from the beginning or from where it stopped.
		string::size_type pos = strYes.find("yes",0);
		if (pos != string::npos){
			restart = true;

			ifstream rfid;
			// seek for the last VTK file generated
			for (int i=1; i<=100; i++){
				char strtmp[512]; sprintf(strtmp,"%s__0-of-1__step-%d.vtk",expofName.c_str(),i);
				rfid.open(strtmp);
				if (rfid.is_open()){
					cout << strtmp << endl;
					rfid.close();
				}
				else{
					if (i==1){
						throw Exception(__LINE__,__FILE__,"Any VTK file has been generated using the provided case file. Please, check carefully for mistypes\n");
					}
					lastpvi = i-2;
					PVI_cumulative = double(lastpvi/100.0);
					char strtmp2[512]; sprintf(strtmp2,"%s__0-of-1__step-%d.vtk",expofName.c_str(),lastpvi-1);
					restartFilename = strtmp2;
					rfid.close();

					// set cumulative simulation time variable: cumTS = summation of all time steps (from the beginning to PVI i-1
					string strline;
					string data[8];
					char strtmp4[512]; sprintf(strtmp4,"%s_simulation-monitor-%d.csv",expofName.c_str(),P_size());
					rfid.open(strtmp4);
					if (!rfid.is_open()){
						throw Exception(__LINE__,__FILE__,"Simulation monitor file could not be opened.\n");
					}
					int pvi = 0;
					getline(rfid,strline);
					while(pvi!=lastpvi){
						rfid >> stepCounter >> data[0] >> pvi  >> cumTS >> data[2] >> data[3] >> data[4] >> data[5] >> data[6] >> data[7];
						//cout << "PVI = " << pvi << " \tstep#: " << stepCounter << "\tcumulative ST: " << cumTS << endl;
					}
					//exit(1);
					break;
				}
			}
			this->setStepOutputFile(lastpvi+1);
			cout << "Step counter: "  << stepCounter << setprecision(8) << "\tCumulate simulation time: " << (double)cumTS << endl;
			cout << "\n\nRestart file name: " << restartFilename << "\n\n";
			//exit(1);
		}

		setPositionToRead(fid,"Dirichlet (Specify: flags and respective pressure values):");
		mapFlag(fid,"dirichlet");

		setPositionToRead(fid,"Injection Wells: 10 - 50");
		getWells(fid);

		setPositionToRead(fid,"Production Wells: 51 - 90");
		getWells(fid);

		setPositionToRead(fid,"Specify water saturation in injection wells:");
		mapFlag(fid,"saturation");

		setPositionToRead(fid,"water saturation at t=0:");
		fid >> _Sw_initial;

		setPositionToRead(fid,"simulation-time:");
		fid >> _ST;

		setPositionToRead(fid,"dvtol:");
		fid >> _dvtol;

		setPositionToRead(fid,"Reservoir Dimensions: L - H:");
		fid >> reservoir_Length >> reservoir_Height;

		pGCData->setReservoirHeight(reservoir_Height);

		fid.close();
		if (!P_pid()) std::cout << "done.\n";
	}

	void SimulatorParameters::loadPhysicalParameters(int dim){
		if (!P_pid()) std::cout << "loading physical parameters... ";
		ifstream fid;
		string str;
		// read physical parameters
		string physical(parametersPath + "/physical.dat");
		fid.open(physical.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"physical.dat could not be opened or it doesn't exist.\n");

		// start reading file
		// ---------------------------------------------------------------------
		setPositionToRead(fid,"Oil density");
		fid >> str; _oil_density = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Water density:");
		fid >> str; _water_density = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Oil viscosity:");
		fid >> str; _oil_viscosity = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Water viscosity:");
		fid >> str; _water_viscosity = strtod(str.c_str(), 0);

		setPositionToRead(fid,"# dom, K and phi:");
		readRockProperties(fid,dim);

		setPositionToRead(fid,"Residual Oil Satutarion");
		fid >> str; _Sor = strtod(str.c_str(), 0);

		setPositionToRead(fid,"Irreducible Water Satutarion");
		fid >> str; _Swr = strtod(str.c_str(), 0);

		setPositionToRead(fid,"type model ID right above:");
		fid >> _ksModel;

		//		printf("_ksModel = %d\n",_ksModel); throw 1;

		fid.close();
		getDomains();
		if (!P_pid()) std::cout << "done.\n";
	}

	void SimulatorParameters::loadSolverParameters(){
		if (!P_pid()) std::cout << "Loading solver parameters... ";
		ifstream fid;
		string str;
		// read solver parameters
		string solver(parametersPath + "/solver.dat");
		fid.open(solver.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"solver.dat could not be opened or it doesn't exist.\n");

		// start reading file
		// ---------------------------------------------------------------------
		setPositionToRead(fid,"absolute convergence tolerance");
		fid >> _abstol;

		setPositionToRead(fid,"relative convergence tolerance for KSP solver");
		fid >> _rtol;

		setPositionToRead(fid,"relative convergence tolerance for external iteration");
		fid >> _rtol2;

		setPositionToRead(fid,"convergence tolerance in terms of the norm of the change in the solution between steps");
		fid >> _dtol;

		setPositionToRead(fid,"maximum number of iterations");
		fid >> _maxit;

		setPositionToRead(fid,"Sets number of iterations at which GMRES, FGMRES and LGMRES restarts:");
		fid >> _Krylov_restart;

		if (!P_pid()){
			std::cout << "\n abstol : " << _abstol << endl;
			std::cout << "   rtol : " << _rtol << endl;
			std::cout << "  rtol2 : " << _rtol2 << endl;
			std::cout << "   dtol : " << _dtol << endl;
			std::cout << "restart : " << _Krylov_restart << endl;
		}

		/*
		 * Load a 'database' to set a preconditioner chosen by user. <private>
		 */
		setPreconditionerDataBase();

		/*
		 * Read solver.dat and check which preconditioner user chose.
		 */
		string::size_type pos;
		int numpc = (int)mapPC.size(); // number of preconditioners presented on data base.
		setPositionToRead(fid,"Define a preconditioner (Default: none)");
		for (int i=0; i<numpc; i++) {
			getline(fid,str,'\n');

			// search for chosen option
			pos = str.find("(x)",0);

			// if marked with 'x', search in 'str' for pre-conditioner available
			if (pos != string::npos){
				std::map<std::string,PCType>::iterator miter = mapPC.begin();
				while (miter != mapPC.end()){
					pos = str.find(miter->first,0);
					if (pos != string::npos){				// miter->first: string name
						pctype = miter->second;				// miter->second preconditioner
						if (!P_pid()) std::cout << "     PC : " << miter->first << std::endl;
					}
					miter++;
				}
			}
		}
		if (!P_pid()){
			std::cout << "done.\n";
			if (getPCType() == PCNONE) std::cout << "\n**** WARNING *** No preconditioner chosen.\n\n" ;
		}

		// read pressure solver scheme for EBFV1 formulation
		str.clear();
		setPositionToRead(fid,"Use <defect-correction> scheme for EBFV1 pressure solver?: (default: matrix-free scheme)");
		getline(fid,str);
		pos = str.find("(x)",0);
		if (pos != string::npos) {
			setUseDefectCorrection();
			if (!P_pid()) printf("Using defect-correction scheme for EBFV1 pressure solver.\n");
		}

		fid.close();
	}

	void SimulatorParameters::loadHighOrderParameter(){
		if (!P_pid()) std::cout << "Loading high order parameters... ";
		ifstream fid;
		string str;
		// read solver parameters
		string slimiters(parametersPath + "/slope_limiters.dat");
		fid.open(slimiters.c_str());
		if ( !fid.is_open() ) throw Exception(__LINE__,__FILE__,"solve-limiters.dat could not be opened or it doesn't exist.\n");

		// start reading file
		// ---------------------------------------------------------------------
		const int HOM = 5;
		std::string theString[HOM];

		/*
		 * Read from file all slope limiter function options
		 */
		setPositionToRead(fid,"Nodal sloper limiter functions (default: MUSCL):");
		for (int i=0; i<HOM; i++) getline(fid,theString[i],'\n');

		/*
		 * Define a "data bank" for node slope limiters function
		 */
		mapSL_node["(x) - Superbee"] = node_Superbee;
		mapSL_node["(x) - Minmod"] = node_Minmod;
		mapSL_node["(x) - Van_Albada"] = node_Van_Albada;
		mapSL_node["(x) - Osher"] = node_Osher;
		mapSL_node["(x) - WoodField"] = node_WoodField;

		std::map<string,NSLF>::iterator mit;

		// define node SL default as MUSCL
		slf_method_Node = node_MUSCL;

		/*
		 * Find which option user has chosen
		 */
		for (int i=0; i<HOM; i++){
			mit = mapSL_node.find(theString[i]);
			if ( mit != mapSL_node.end() ) slf_method_Node = mit->second;
			theString[i].clear();
		}

		/*
		 * Read from file all slope limiter function options
		 */
		setPositionToRead(fid,"Edge sloper limiter functions (default: MUSCL):");
		for (int i=0; i<HOM; i++) getline(fid,theString[i],'\n');

		/*
		 * Define a "data bank" for edge slope limiters function
		 */
		mapSL_edge["(x) - Superbee"] = SUPERBEE;
		mapSL_edge["(x) - Minmod"] = MINMOD;
		mapSL_edge["(x) - Van Albada"] = VAN_ALBADA;
		mapSL_edge["(x) - Van Albada"] = OSHER;
		mapSL_edge["(x) - WoodField"] = WOODFIELD;

		// define edge SL default as MUSCL
		slf_method_Edge = MUSCL;
		/*
		 * Find which option user has chosen
		 */
		std::map<string,ESLF>::iterator mite;
		for (int i=0; i<HOM; i++){
			mite = mapSL_edge.find(theString[i]);
			if ( mite != mapSL_edge.end() ) slf_method_Edge = mite->second;
		}

		/*
		 * Modified Taylor expansion series extrapolating nodal saturation value on volume
		 * control interfaces. Mark an 'x' (without quotes) inside the parentheses to define
		 * the coefficient 'k' for the saturation extrapolation.
		 */

		// define koef as 0 <default>
		koef = .0;
		map_koef["(x) -  -1,   método de ponderação a montante de segunda ordem"] = -1.;
		map_koef["(x) - 1/3, método de ponderação a montante de terceira ordem"] = 1./3.;
		map_koef["(x) -   1,   método de diferenças centradas de três pontos"] = 1.;

		/*
		 * Read from file the options
		 */
		setPositionToRead(fid,"Taylor expansion coefficient - (default: koef = 0, método de Fromm):");
		for (int i=0; i<4; i++) getline(fid,theString[i],'\n');

		/*
		 * Find which option the use has chosen
		 */
		std::map<string,double>::iterator mit_koef;
		for (int i=0; i<HOM; i++){
			mit_koef = map_koef.find(theString[i]);
			if ( mit_koef != map_koef.end() ){
				koef = mit_koef->second;
				if (!P_pid()) cout << "koef = " << koef << endl;
			}
		}


		// Delta value for woodfield
		_WFdelta = 0.2;
		setPositionToRead(fid,"Delta value for Woodfield:");
		fid >> _WFdelta;

		fid.close();
		if (!P_pid()) std::cout << "done.\n";
	}

	void SimulatorParameters::setPositionToRead(ifstream &fid, string str2){
		int count = 0;
		string str1;
		do{
			getline(fid,str1);
			if (++count > 100){
				if (!P_pid()){
					cout << "\n\nWARNING: Could not find <" << str2
							<<">.\nProbably you are using an out of date input file."
							" Default value will be used instead.\n\n";
				}break;
			}
		}while( str1.compare(str2) );
	}

	void SimulatorParameters::readRockProperties(ifstream &fid, int dim){
		RockProperties *rockprop;
		double num;
		string str;
		istringstream stream1;
		std::list<double> K_tmp;

		/*
		 * This is a way to use the simulator to evaluate elliptic equation without screw-up
		 * the input data procedure.
		 */
	#ifdef CRUMPTON_EXAMPLE
		RockProperties *rockprop1 = new RockProperties;
		RockProperties *rockprop2 = new RockProperties;

		rockprop1->K = new double[4];	// get permeability tensor
		rockprop1->K[0] = 1; rockprop1->K[1] = 0;
		rockprop1->K[2] = 0; rockprop1->K[3] = 1;
		mapRockProp[3300] = rockprop1;

		rockprop2->K = new double[4];	// get permeability tensor
		rockprop2->K[0] = (double)2*ALPHA; rockprop2->K[1] = (double)1*ALPHA;
		rockprop2->K[2] = (double)1*ALPHA; rockprop2->K[3] = (double)2*ALPHA;
		mapRockProp[3301] = rockprop2;

	#else
		while( true ){
			getline(fid,str);
			if (!str.compare("end")) break;
			rockprop = new RockProperties;
			const int dom = atoi(str.c_str());	// get domain flag
			rockprop->K = new double[dim*dim];	// get permeability tensor
			// read permeability tensor
			getline(fid,str);
			stream1.str(str);
			int count =0;
			// check if file data size match to mesh dimensions
			// 2-D -> K[2][2] or 3-D -> K[3][3]
			if (!P_pid()) cout << "\nK = ";
			while( stream1 >> num ){
				K_tmp.push_back(num);
				count++;
				if (!P_pid()) cout << num << " ";
			}
			if ( count != dim*dim )
				throw Exception(__LINE__,__FILE__,"Mesh dimension and permeability tensor size do not match!\n");

			std::copy(K_tmp.begin(),K_tmp.end(),rockprop->K);

			// get rotation tensor
			//			getline(fid,str);
			//			double gamma = strtod(str.c_str(), 0)*pi/180.0;
			//			cout << "rotation tensor = "  << gamma << endl; STOP();

			// get porosity
			getline(fid,str);
			rockprop->porosity = strtod(str.c_str(), 0);
			mapRockProp[dom] = rockprop;
			K_tmp.clear();
			stream1.clear();
		}
		if (!P_pid()) cout << endl;
		if (mapRockProp.size()==0)
			throw Exception(__LINE__,__FILE__,"You must provide rock properties\n");
	#endif
	}

	void SimulatorParameters::mapFlag(ifstream &fid, string whatmap){
		string str;
		BdryConditionType *bct;

		while( 1 ){
			fid >> str;
			if (!str.compare("end")) break;
			cout << str << endl;
			int flag = atoi(str.c_str());
			fid >> str;
			double val = strtod(str.c_str(), 0);
			cout << str << endl;

			// for every flag number create a new BdryConditionType pointer and add inform type of condition and its value.
			bct = new BdryConditionType;
			bct->type = whatmap;
			bct->val = val;

			// two map container are used: one for saturation and other to dirichlet/neumann boundary conditions
			if (!whatmap.compare("saturation")){
				printf("mapSaturation[flag] = %d\n",flag);
				mapSaturation[flag] = bct;
				printf("val = %f\n",val);
			}
			else{
				mapBC[flag] = bct;
				printf("flag: %d has boundary condition: %f\n",flag,bct->val);
			}
		}
	}

	void SimulatorParameters::getWells(ifstream &fid){
		string str;
		while( true ){
			fid >> str;
			if (!str.compare("end")) break;
			int flag = atoi(str.c_str());

			fid >> str;
			double val = strtod(str.c_str(), 0);
			// associate a vector of size 2 to each well flag
			// data[0] = Q, total flow rate
			// data[1] = V, total well volume (to be calculated later)
			//			dblarray data(2);
			//			data[0] = val;
			//			MWells[flag] = data;

			well_info.isInjection = false;
			well_info.isProduction = true;
			well_info.flowRate = val;
			if (val>.0){
				well_info.isInjection = true;
				well_info.isProduction = false;
			}
			MWells[flag] = well_info;
			//cout << "Flag: " << flag << "\tFlow: " << val << endl;
		}
	}

	string SimulatorParameters::getFilename(string meshName){
		char s1[256], s4[256];
		generateFilename(meshName,s1);
		sprintf(s4,"%s_CFL%.2f-%d-of-%d",s1,CFL(),P_pid(),P_size());
		string s5(s4);
		return s5;
	}

	void SimulatorParameters::generateFilename(const string &meshFilename, char *filename){
		string str(meshFilename);
		string::size_type start = str.find_last_of("/");
		start = (start == string::npos)?0:start;
		string::size_type end = str.find_last_of(".");
		char buffer[256];
		memset( buffer, '\0', 256 );
		str.copy(buffer,end-start,start);
		sprintf(filename,"%s",buffer);
	}

	void SimulatorParameters::setPreconditionerDataBase(){
		mapPC["(x) JACOBI"]    = PCJACOBI;
		mapPC["(x) SOR"]       = PCSOR;
		mapPC["(x) LU"]        = PCLU;
		mapPC["(x) SHELL"]     = PCSHELL;
		mapPC["(x) BJACOBI"]   = PCBJACOBI;
		mapPC["(x) MG"]        = PCMG;
		mapPC["(x) EISENSTAT"] = PCEISENSTAT;
		mapPC["(x) ILU"]       = PCILU;
		mapPC["(x) ICC"]       = PCICC;
		mapPC["(x) ASM"]       = PCASM;
	}

	void SimulatorParameters::setLocalNodeIDNumbering(pMesh theMesh){
//		int k = 0;
//		pEntity elem;
//		std::set<pEntity> nodesDomain;
//		std::set<int>::iterator iter = setOfDomains.begin();
//		for(;iter!=setOfDomains.end();iter++){
//			if (theMesh->getDim()==2){
//			FIter fit = M_faceIter( theMesh );
//			while ( (elem = FIter_next(fit)) ) {
//				if ( !theMesh->getRefinementDepth(elem) ){
//					int flag = EN_getFlag(elem);
//					if ( flag==*iter ){
//						for (int i=0;i<3;i++){
//							nodesDomain.insert(elem->get(0,i));
//						}
//					}
//				}
//			}
//			FIter_delete(fit);}
//			else{
//				RIter rit = M_regionIter( theMesh );
//				while ( (elem = RIter_next(rit)) ) {
//					if ( !theMesh->getRefinementDepth(elem) ){
//						int flag = EN_getFlag(elem);
//						if ( flag==*iter ){
//							for (int i=0;i<4;i++){
//								nodesDomain.insert(elem->get(0,i));
//							}
//						}
//					}
//				}
//				RIter_delete(rit);
//			}
//
//			int i = 0;
//			char tag[4]; sprintf(tag,"%d",k++);	// flag node with as k-th domain
//			std::set<pEntity>::iterator iter_global = nodesDomain.begin();
//			for(;iter_global!=nodesDomain.end();iter_global++){
//				EN_attachDataInt(*iter_global,MD_lookupMeshDataId(tag),i++);
//			}
//			nodesDomain.clear();
//		}
	}

}