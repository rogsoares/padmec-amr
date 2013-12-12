#include "SIMULATION_core.h"

namespace PRS           // PRS: Petroleum Reservoir Simulator
{
SIMULATION_core::SIMULATION_core(){
	pElliptic_eq = 0;
	pHyperbolic_eq = 0;
	theMesh = MS_newMesh(0);
	//pMeshAdapt = NULL;
}

SIMULATION_core::~SIMULATION_core(){
}

int SIMULATION_core::initialize(int argc, char **argv){
	/*
	 * Initialize parallel libraries
	 */
	ParUtil::Instance()->init(argc,argv);
	PetscInitialize(&argc,&argv,(char *)0,(char *)0);

	if ( argc==1 ){
		throw Exception(__LINE__,__FILE__, Exception::INIT_ERROR );
	}
	simFlag = atoi(argv[1]);

	printSimulationHeader();

	/*
	 * Initialize simulation pointers
	 */
	pPPData = new PhysicPropData;
	pGCData = new GeomData;
	pSimPar = new SimulatorParameters(theMesh);
	pMData = new MeshData(pSimPar,theMesh);

	/*
	 * Load data from files
	 */
	pSimPar->inputDataFromFiles(pGCData);

	/*
	 * If adaptation required, initialize:
	 * 		1- Adaptation Pointer (h-Refinement or adaptative remeshing)
	 * 		2- Error Analysis Pointer
	 * 		3- Interpolation Function
	 */
	//if (pSimPar->userRequiresAdaptation()){
		// 1- Adaptation Pointer (h-Refinement or adaptative remeshing)
		// ----------------------------------------------------------------------------------
		int dim = pGCData->getMeshDim();
		if (dim<2 || dim>3){
		  throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
		}
		if (dim==3){
			throw Exception(__LINE__,__FILE__,"Only 2-D adaptation allowed for while.");
		}

		// 2- Error Analysis Pointer
		// ----------------------------------------------------------------------------------
		pErrorAnalysis = new ErrorAnalysis_2D;
		pIData = new InterpolationDataStruct;

		pIData->getLevelOfRefinement = ErrorAnalysis::getLevelOfRefinement;
		pIData->numFields = 2;
		pIData->pGetDblFunctions = new GetDblFunction[2];
		pIData->pSetDblFunctions = new SetDblFunction[2];
		pIData->pGetDblFunctions[0] = pPPData->getPressure;
		pIData->pGetDblFunctions[1] = pPPData->getSaturation;
		pIData->pSetDblFunctions[0] = pPPData->setPressure;
		pIData->pSetDblFunctions[1] = pPPData->setSaturation;
		pIData->m2 = MS_newMesh(0);

		//switch( pSimPar->getRefStrategy() ){
		switch( pSimPar->getRefStrategy() ){
		case H_REFINEMENT:
			pMeshAdapt = new H_Refinement_2D;
			pIData->isElementSpecial = H_Refinement_2D::isElementSpecial;
			break;
		case ADAPTIVE_REMESHING:
			pMeshAdapt = new AdaptiveRemeshing(argc, argv);
			break;
 		case RH_REFINEMENT:
 			//pMeshAdapt = new RH_Refinement(argc, argv);
 			break;
		default:
			throw Exception(__LINE__,__FILE__,"Unknown adaptation strategy.");
		}

//		// 3- Interpolation Function
//		// ----------------------------------------------------------------------------------
//				switch ( pSimPar->getInterpolationMethod() ){
//				case h_REFINEMENT:
//					pInterpolateData = hRefinement;
//					break;
//				case LINEAR:
//					pInterpolateData = Linear;
//					break;
//				case QUADRATIC:
//					pInterpolateData = Quadratic;
//					break;
//		//		case ADAPTATIVE:
//		//			pInterpolateData = Adaptative;
//		//			break;
//		//		case CONSERVATIVE:
//		//			pInterpolateData = Conservative;
//		//			break;
//		//		case PURE_INJECTION:
//		//			pInterpolateData = PureInjection;
//		//			break;
//		//		case HALF_WEIGHTING:
//		//			pInterpolateData = HalfWeighting;
//		//			break;
//		//		case FULL_WEIGHTING:
//		//			pInterpolateData = FullWighting;
//		//			break;
//				default:
//					throw Exception(__LINE__,__FILE__,"Interpolation method unknown. Exiting....");
//				}
	//}


	/*
	 *  Initialization procedure based on previous loaded data from file:
	 *  	- well flow
	 *  	- physical properties
	 *  	- boundary and initial conditions, mapping
	 */
	pSimPar->initialize(pGCData,theMesh);
	pPPData->initialize(pMData,pSimPar,theMesh,false);
	pMData->initialize(theMesh,pGCData);

	if (!P_pid()) printf("Number of processes required: %d\n",P_size());

	// If restart required, physical properties come from a specified vtk file provided by user.
// 	if (pSimPar->useRestart()){
// 		restartSimulation(pSimPar,pPPData,theMesh);
// 	} 

	/*
	 * Oil production output
	 */
	if ( pSimPar->rankHasProductionWell() ) {
		string path = pSimPar->getOutputPathName();
		char tmp[256]; sprintf(tmp,"%s_oil-production-%d.xls",path.c_str(),P_size());
		string FileName(tmp);
		pOilProduction = new OilProductionManagement(FileName,pSimPar->getInitialOilVolume(),pSimPar->getTotalInjectionFlowRate());
	}

	/*
	 *  Initialize elliptic and hyperbolic solver pointers
	 */
	pElliptic_eq = init_EllipticSolverPointer( pSimPar->getEllipticSolver() );
	pHyperbolic_eq = init_HyperbolicSolverPointer( pSimPar->getHyperbolicSolver() );
	return 0;
}


Elliptic_equation* SIMULATION_core::init_EllipticSolverPointer(int elliptic_method){
	switch ( elliptic_method ){
	case 1:
		return new EBFV1_elliptic(theMesh,pPPData,pSimPar,pGCData,pMData);
	default:
		throw Exception(__LINE__,__FILE__,"Could not initialize a pointer to pElliptic_eq. Unknown method.\n");
	}
}

Hyperbolic_equation* SIMULATION_core::init_HyperbolicSolverPointer(int hyperbolic_method){
	switch ( hyperbolic_method ){
	case 1:
		return new EBFV1_hyperbolic(theMesh,pPPData,pSimPar,pGCData,pMData,pOilProduction,pErrorAnalysis);
	case 2:
		return  new EBFV1_hyperbolic_adaptative(theMesh,pPPData,pSimPar,pGCData,pMData,pOilProduction,pErrorAnalysis);
	default:
		throw Exception(__LINE__,__FILE__,"Could not initialize a pointer to pHiperbolic_eq. Unknown method.\n");
	}
}
///#define TRACKING_PROGRAM_STEPS
void SIMULATION_core::updatePointersData(pMesh theMesh){
	cout<< "UPDATEPOINTERS"<<endl;
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: updating Pointers\tIN\n";
#endif
	// starting deallocating data related to simulation pointers
	cout << __LINE__ << endl;
	
	pSimPar->deallocateData(); cout << __LINE__ << endl;
	pPPData->deallocateData(pSimPar); cout << __LINE__ << endl;
	FIter fit = M_faceIter( theMesh );
	FIter_delete(fit);
	cout << __LINE__ << endl;
	pMData->deallocateData(); cout << __LINE__ << endl;
	
	
	
	// initialize them once more
	pSimPar->initialize(pGCData,theMesh); cout << __LINE__ << endl;
	pPPData->initialize(pMData,pSimPar,theMesh,true); cout << __LINE__ << endl;
	pMData->initialize(theMesh,pGCData); cout << __LINE__ << endl;
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: updating Pointers\tOUT\n";
#endif
	
}

int SIMULATION_core::finalize(){
	/*
	 * Write to file oil production output. Only rank 0 is in charge of it.
	 */
	string path = pSimPar->getOutputPathName();

	char tmp[256]; sprintf(tmp,"%s_PETSc_summary_nproc%d.log",path.c_str(),P_size());
	PetscErrorCode ierr = PetscLogPrintSummary(MPI_COMM_WORLD,tmp); CHKERRQ(ierr);

	// free memory
	delete pElliptic_eq;
	delete pHyperbolic_eq;
	delete pPPData;
	delete pSimPar;
	delete pGCData;
	delete pMData;

	// finalize MPI
	ParUtil::Instance()->Finalize();
	return 0;
}
}
