#ifndef SIMULATORPARAMETERS_H_
#define SIMULATORPARAMETERS_H_

#include "GeomData.h"
#include "pre-processorList.h"
#include "SimParStrucs.h"


typedef std::set<int> setNodes;
enum RefinementStrategies{H_REFINEMENT, ADAPTATIVE_REMESHING};
enum INTERPOLATION_OPTIONS {h_REFINEMENT,LINEAR, QUADRATIC, ADAPTATIVE, CONSERVATIVE, PURE_INJECTION, HALF_WEIGHTING, FULL_WEIGHTING};

namespace PRS{

/**
 * SimulatorParameters class provides all IO data stream during simulation.
 * Numeric, physical, pre-processor data input are made through this class.
 * When a specific class object is not applied to do so within SimulatorPara-
 * meters class, it performs the data reading itself.
 * VTK printing files are also made using SimulatorParameters.
 */
class SimulatorParameters{
public:

	SimulatorParameters();
	SimulatorParameters(pMesh theMesh);
	~SimulatorParameters();

	// read data from files.
	void inputDataFromFiles(GeomData *);
	void initialize(GeomData *, pMesh);
	void deallocateData();
	void loadNumericParameters(GeomData*);
	void loadPhysicalParameters(int);
	void loadSolverParameters();
	void loadHighOrderParameter();
	void loadMeshAdaptationParameters();

	/** load_preprocessor data is based on user decision from numeric.dat.
	 * This function uses a void pointer object. The main reason is to use a
	 * pointer to GeomData which has been designed to set/get
	 * geometric data associated to mesh entities. But it not a rule and any
	 * kind of object may be passed through load_preprocessorData avoiding
	 * additional code implementation. Preprocessor file name is defined on
	 * numeric.dat file.
	 */
	void load_preprocessorData(void*);

	// returns true if fimulation time was reached
	bool finishSimulation();

	// finishes simulation for some specific reason
	void stopSimulation() { stop_simulation = false; }

	// called inside loadParameters() to read rock properties on the same file (weird!)
	void readRockProperties(ifstream &fid, int);

	int getEllipticSolver() const { return ellipticSolver; }
	int getHyperbolicSolver() const { return hyperbolicSolver; }
	string exportFilename() const { return expofName; }
	string getOutputPathName() const { return expofName; }
	string prepFilename() const { return prepfName; }

	/*
	 * numeic parameters
	 * --------------------------------------------------------------------
	 */
	double CFL() const { return _CFL; }
	double getBC_Value(const int &flag);	// return value associated to flag
	bool isNodeFree(const int &flag);		// returns false if dirichlet
	double getWellTimeDiscretion() const { return _timeStepWell; } // Returns a number that is used to divide
	// the current time step (fractional time-step)


	/**
	 *  solver settings parameters
	 *  -------------------------------------------------------------------------
	 */
	double abstol() const { return _abstol; }
	double rtol() const { return _rtol; }
	double rtol2() const { return _rtol2; }
	double dtol() const { return _dtol; }
	int maxit() const { return _maxit;	}
	bool useDefectCorrection() const { return EBFV1_pressureSolver_scheme; }
	void setUseDefectCorrection() { EBFV1_pressureSolver_scheme = true; }
	PCType getPCType() const { return pctype; }
	double getPVIincrement() const { return PVI_increment; }
	int getKrylov_restart() const { return _Krylov_restart; }


	// physical parameters
	// -------------------------------------------------------------------------
	double oilDensity() const { return _oil_density; }
	double waterDensity() const { return _water_density; }
	double oilViscosity() const { return _oil_viscosity; }
	double waterViscosity() const { return _water_viscosity; }
	double getPorosity(const int &dom);
	const double* getPermeability(const int &dom);
	bool is_K_Isotropic() const { return K_Isotropic; }
	int ksModel() const { return _ksModel; }	// return relative permeability model

	// initial conditions
	// -------------------------------------------------------------------------
	double Sw_initial() const { return _Sw_initial; }
	double Sor() const { return _Sor; }
	double Swr() const { return _Swr; }
	double getInitialOilVolume() const { return _IOV; }


	/**
	 * Well management
	 * --------------------------------------------------------------------
	 */
	MapWells MWells;
	// for each well flag will be several nodes flagged
	// as edge loop is performed a node can be counted twice. So a set
	// container is used to store the ID from node
	typedef set<int> setNodesOnWell;
	map<int,setNodesOnWell> mapNodesOnWell;
	map<int,setNodesOnWell>::iterator MNOWIter;
	bool isInjectionWell(pEntity node){ return isInjectionWell(GEN_tag(node->getClassification())); }
	bool isProductionWell(pEntity node){ return isProductionWell(GEN_tag(node->getClassification())); }
	bool isInjectionWell(int flag) const;
	bool isProductionWell(int flag) const;
	double getTotalInjectionFlowRate() const;
	bool hasNodeWell(const int &flag);	// returns true if node has a flag for well
	double getWellVolume(int well_flag) const;	// returns the sum of nodal volumes of the nodes well
	void getWells(pMesh theMesh, int);
	bool rankHasProductionWell() const { return _rankHasProductionWell; }
	void checkIfRankHasProductionWell();
	void weightWellFlowRateByVolume(GeomData *);
	double getFlowrateValue(int flag) const;
	double getInitialSaturation(pEntity);

	// reservoir geometric dimension used to make some physical properties dimensionless
	//void getReservoirGeometricDimensions(double &L, double &H) const;
	void getReservoirGeometricDimensions(double &L, double &H) const {
		H = reservoir_Height; L = reservoir_Length;
	}
	double dimensionlessFactorForDeltaTImplicit(double dom);
	double getDVTOL() const { return _dvtol; }
	void setInitialOilVolume(GeomData*);	// this work for multi-domains

	set<int> setOfDomains;
	int getNumDomains() const { return setOfDomains.size(); }

	/*
	 * Counts how many nodes belong to each domain.
	 */
	void setNumElementDomain(pMesh theMesh);

	const int *getNumNodesDomain()const{
		return numNodesDom;
	}
	void getNumEdgesDomain(std::set<int> &setNumEdgesDom){
		int numDomains = setOfDomains.size();
		for(int i=0;i<numDomains;i++){
			setNumEdgesDom.insert(numEdgesDom[i]);
		}
	}

	SIter_const setDomain_begin() const { return setOfDomains.begin(); }
	SIter_const setDomain_end() const { return setOfDomains.end(); }

	/*
	 * High order managemet
	 * --------------------------------------------------------------------
	 */
	bool useHOApproximation() const { return useHOApp; }
	NSLF getNodeSlopeLimitFunc() const { return slf_method_Node; }
	ESLF getEdgeSlopeLimitFunc() const { return slf_method_Edge; }
	double get_koef() const{ return koef; }
	double get_WoodfieldDelta() const { return _WFdelta; }

	// skip comments or any other unnecessary string until a specific point be reached
	void setPositionToRead(ifstream &fid, string str);

	/*
	 * Restart stuff
	 * --------------------------------------------------------------------
	 */
	bool useRestart() const { return restart; }
	string getRestartFilename() const { return restartFilename; }
	double getAccumulatedSimulationTime() const { return accSimTime; }
	void setAccumulatedSimulationTime(double time_increment) { accSimTime += time_increment; }
	int getStepOutputFile() const { return vtk_step; }
	void setStepOutputFile(int s) { vtk_step = s; }
	void setCPU_time(double cput) { cpu_time=cput; }
	double getCPU_time() const { return cpu_time; }
	void getSimParFiles();


	/*
	 * Simulation time monitoring
	 * --------------------------------------------------------------------
	 */
	double getPVIaccumulated() const { return PVI_accumulated; }
	void setPVIaccumulated(double pviacc) { PVI_accumulated = pviacc; }
	void setTStepNumber(int tsn) { tsnumber=tsn; }
	int getTStepNumber() const { return tsnumber; }
	double getSimTime() const { return _ST; }
	double getPVI() const { return _PVI; }


	/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	 * If it's desired to print out a VTK file every fraction of PVI or every
	 * n days, timeStep must be corrected to reproduce the exact time evolution.
	 * For example, if every 0.1PVI a new VTK file must be generated it has:
	 *
	 * sum(dt) = n*PVI
	 * sum(dt) = dt_1 + dt_2 + ... + dt_n
	 *
	 * if ( sum(dt)>n*PVI ) dt_n = dt_n - (ST - n*PVI)
	 *
	 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
	void correctTimeStep(double&);
	void printOutVTK(pMesh theMesh, void *pData1, void *pData2, void *pData3, pFunc_PrintVTK printVTK);
	double getPrintOutVTKFrequency() const { return vtk_time_frequency; }
	void setPrintOutVTKFrequency(double vtf) { vtk_time_frequency = vtf; }
	void updatePrintOutVTKFrequency();
	string getFilename(string f);
	void generateFilename(const string &meshFilename, char *filename);


	/*
	 * Adaptation parameters
	 */
	void setAdaptation() { _doAdaptation = true; }
	bool userRequiresAdaptation() const { return _doAdaptation; }

	void setTol1(double tol) { tol1 = tol; }
	double getToleranceForAllElements() const {  return tol1; }

	void setTol2(double tol) { tol2 = tol; }
	double getToleranceForAllElements_excludingSingularities() const {  return tol2; }

	void setMax_2D(int max) { numMax_2D = max; }
	// avoid excessive element subdivisions due to singularities
	int getMax2D() const {  return numMax_2D; }

	void setMax_3D(int max) { numMax_3D = max; }
	// avoid excessive element subdivisions due to singularities
	int getMax3D() const {  return numMax_3D; }

	void setNumSubdivision_perStep(int max) { numMaxSubdivisions = max; }

	// avoid excessive element subdivisions to get a less refined mesh
	int getNumSubdivision_perStep() const {  return numMaxSubdivisions; }

	// return a set of nodes ID from domain dom
	void getNodesDomain(int dom, setNodes& setnodes){
		setnodes = _mapNodesDomain[dom];
	}

	void setLocalNodeIDNumbering(pMesh);

	// return a set of nodes ID from domain dom
	void getLocalNodeIDNumbering(pEntity node, char* tag, int &local_ID){
		EN_getDataInt(node,MD_lookupMeshDataId(tag),&local_ID);
	}

	int getLocalNodeIDNumbering(pEntity node, char* tag){
		int id;
		EN_getDataInt(node,MD_lookupMeshDataId(tag),&id);
		return id;
	}

	RefinementStrategies getRefStrategy() const{
		return refstrategy;
	}

	void setInterpolationMethod(INTERPOLATION_OPTIONS im){
		_intpmethod = im;
	}

	INTERPOLATION_OPTIONS getInterpolationMethod() const{
		return _intpmethod;
	}

private:

	RefinementStrategies refstrategy;

	// stores how many time steps are performed within every fraction of PVI
	std::list<int> TSCountingList;

	/*
	 * The following three variables below store user's decision to use, from
	 * a set of options, a specific implementation. The chosen option will be
	 * executed through a switch/case coding.
	 */
	NSLF slf_method_Node;		// related to high order method
	ESLF slf_method_Edge;		// related to high order method
	FRACTIONALFLUX fw_method; // related to fractional flux

	std::map<string,NSLF> mapSL_node;
	std::map<string,ESLF> mapSL_edge;
	std::map<string,double> map_koef;

	/*
	 *  koef: coefficient for Sw extrapolation used in high order implementation
	 *  Darlan's thesis (Carvalho, 2005), Eqs. (4.78) and (4.79), p92.
	 */
	double koef;

	// delta for woodfield high order method
	double _WFdelta;

	/*
	 * Get from file number of physical domains (heterogeneous media).
	 */
	void getDomains();

	// store flags and their values into a specific map container
	void mapFlag(ifstream &fid, string whatmap);
	//void checkIfRankHasProductionWell();
	MapRockProperties mapRockProp;
	MapFlag mapBC;
	MapFlag mapSaturation;
	int numNodesOnInjectWell;

	// numeric parameters
	// -------------------------------------------------------------------------
	double _CFL;
	unsigned int _timeStepWell;
	unsigned int TimeStepCounter;

	PCType pctype;
	double _abstol;	// the absolute convergence tolerance
	double _rtol;	// the relative convergence tolerance
	double _rtol2;	// the relative convergence tolerance for external iteration
	double _dtol;	// the divergence tolerance
	int _maxit;		// maximum number of iterations
	bool EBFV1_pressureSolver_scheme;

	// physical parameters
	// -------------------------------------------------------------------------
	double _phi;
	double _oil_density;
	double _water_density;
	double _oil_viscosity;
	double _water_viscosity;
	double _Sw_initial;	// water saturation at t=0
	double _Sor;		// residual oil saturation
	double _Swr;		// irreducible water saturation
	int _ksModel;		// relative permeability model
	double _IOV;			// Initial Oil Volume

	bool _rankHasProductionWell;


	double _PVI;	// pore volume injected
	double _ST;		// simulation time
	bool well;		// inside isNodeFree() well is set to inform if node has a well or not
	int exportIter;

	void printParameters();

	/*
	 * Get from file and mesh wells (injection and production)
	 */
	void getWells(ifstream &fid);

	bool stop_simulation;
	int _Krylov_restart;
	double _dvtol;
	double reservoir_Height;
	double reservoir_Length;
	double initialOilVolume;
	int ellipticSolver;
	int hyperbolicSolver;

	// pointer for mesh data structure
	pMesh theMesh;

	// strings for input/output files
	string expofName;
	string prepfName;

	/*
	 * If high order approximation is required by user, set useHOApp=true
	 */
	void setUseHOApproximation() { useHOApp = true; }
	bool useHOApp;

	/*
	 * Restart stuff
	 */
	bool restart;
	string restartFilename;
	double accSimTime;
	int vtk_step;
	double cpu_time;
	int tsnumber;

	/*
	 * string variable for file name
	 */
	string numeric_Filename;
	string physical_Filename;
	string slopeLimiter_Filename;
	string solver_Filename;

	/// Says when the next vtk file must be printed out
	double vtk_time_frequency;
	double PVI_increment;
	double PVI_accumulated;
	bool allowPrintingVTK;
	bool firstVTKupdate;

	// map preconditioners
	void setPreconditionerDataBase();
	std::map< std::string, PCType> mapPC;

	// says where numeric.dat physical.dat slope_limiters.dat solve.dat are located
	void getParametersDirectoryPath();
	string parametersPath;

	void checkPermeabilityTensor();
	bool K_Isotropic;

	WellInfo well_info;

	/*
	 * Adaptation parameters
	 */
	bool _doAdaptation;		// Use mesh adaptation:
	double tol1; 			// Error tolerance for all mesh elements:
	double tol2; 			// Error tolerance for all mesh elements excluding those on singularities regions:
	int numMax_2D;			// Maximum number of 2D element subdivisions (recommended: 4)
	int numMax_3D;			// Maximum number of 3D element subdivisions (recommended: 3)
	int numMaxSubdivisions; // Maximum number of refinement steps:


	// domains
	int *numNodesDom;
	int *numEdgesDom;


	std::map<int,setNodes> _mapNodesDomain;	      	// for each domain -> there is a set o of global node IDs
	//std::map<int,setNodes> _mapLocalID_perDomain;	// for each domain -> there is a set o of local node IDs

	INTERPOLATION_OPTIONS _intpmethod;	// interpolation method for adaptative mesh refinement
};
}
#endif /*SimulatorParameters_H_*/
