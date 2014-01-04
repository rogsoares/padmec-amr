#ifndef PHYSICALPROPERTIESDATA_H_
#define PHYSICALPROPERTIESDATA_H_

#include "MeshData.h"
#include "Matrix.h"
#include "GeomData.h"

// NASTY THING
// I need a variable to store pressure gradients to be used into a static member function
// If this variable is a class member, it will not compile!
// Global variable <argh!>. I swear I would not do that!
extern Matrix<double> *pGrad_matrix;
extern Matrix<double> *SwGrad_matrix;
extern Matrix<double> nonvisc_matrix;


namespace PRS		// PRS: Petroleum Reservoir Simulator
{

// properties associated to edge: EX.: velocity
// =========================================================================
struct EdgePhysicalProperties{
	map<int,dblarray> v_new;
	map<int,dblarray> v_old;
	bool projected;
};

// properties associated to edge: EX.: pressure, ssaturation, grad(u), ...
// =========================================================================
struct NodePhysicalProperties{
	map<int,dblarray> pw_Grad;
	map<int,dblarray> Sw_Grad;
	double pw;
	double Sw;
	double Sw_old;
	double Sw_min;
	double Sw_max;
	double S_Limit;
	double nonvisc;
	double volume;
	bool projected;
};


/**
 * PhysicPropData class provides methods to set/get physical proper-
 * ties such as pressure, saturation and their respective gradients as well
 * as rock-fluid properties (fractional flow and total mobilities). Fluid ve-
 * locity is also managed by this class. All properties are associated to a
 * mesh entity (node or edge) through the class AttachData that is in charge
 * of store all data (physical, numeric, geometric,...) to mesh entities.
 */
class PhysicPropData : public AttachData{
public:

	PhysicPropData();
	~PhysicPropData();

	void initialize(MeshData *, SimulatorParameters *, pMesh, bool, GeomData* pGCData);
	void deallocateData(SimulatorParameters *pSimPar);

	/*
	 *  set/get gradient pressure in new format. As pGrad_matrix is declared outside of class definition, it will not give errors in a static
	 *  member function. If pGrad_matrix was declared as a member of class, a class object will be required do be declared inside function, which
	 *  will call the constructor everytime set/get gradient was called.
	 *
	 *  dom - domain flag. As computation is performed by domain, it follows crescent order (0,1,2,3,...). Then, domain's flag is not necessary
	 *  to know, only the order (if it is the first, second, ....)
	 *  row - local node numbering. Do not use global ID or it will not work!!!!
	 *  grad - pointer to an array of three positions. Ex.: double p[3];
	 */
	static void set_pw_Grad(pVertex node, int dom, int row, const double* grad){
		pGrad_matrix[dom].setValue(row,0,grad[0]);
		pGrad_matrix[dom].setValue(row,1,grad[1]);
		pGrad_matrix[dom].setValue(row,2,grad[2]);
	}

	static void get_pw_Grad(pVertex node, int dom, int row, double* grad){
		grad[0] = pGrad_matrix[dom].getValue(row,0);
		grad[1] = pGrad_matrix[dom].getValue(row,1);
		grad[2] = pGrad_matrix[dom].getValue(row,2);
	}

	/*
	 * The same explanation for the functions below. Those are for saturation.
	 */
	static void set_Sw_Grad(pVertex node, int dom, int row, const double* grad){
//		row =  1;
//		static double time = .0;
//		double t1 = MPI_Wtime();
//		SwGrad_matrix[0].setValue(row,0,grad[0]);
//		SwGrad_matrix[0].setValue(row,1,grad[1]);
//		SwGrad_matrix[0].setValue(row,2,grad[2]);
//		double t2 = MPI_Wtime();
//		time += t2-t1;
//		cout << "time: " << time << endl;

		for(int i=0;i<3;i++){
			char string[256]; sprintf(string,"Swgrad_%d_%d",dom,i);
			EN_attachDataDbl(node,MD_lookupMeshDataId(string),grad[i]);
		}
	}

	static void get_Sw_Grad(pVertex node, int dom, int row, double* grad){
//		static double time = .0;
//		double t1 = MPI_Wtime();
//		grad[0] = SwGrad_matrix[0].getValue(row,0);
//		grad[1] = SwGrad_matrix[0].getValue(row,1);
//		grad[2] = SwGrad_matrix[0].getValue(row,2);
//		double t2 = MPI_Wtime();
//		time += t2-t1;
//		cout << "time: " << time << endl;
		for(int i=0;i<3;i++){
			char string[256]; sprintf(string,"Swgrad_%d_%d",dom,i);
			EN_getDataDbl(node,MD_lookupMeshDataId(string),&grad[i]);
		}
	}

	static void setNonViscTerm(int row, double nonvisc){
		nonvisc_matrix.setValue(row,nonvisc);
	}

	static double getNonViscTerm(int row){
		return nonvisc_matrix.getValue(row);
	}


	static void setNonViscTerm(pEntity node, double nonvisc){
		EN_attachDataDbl(node,MD_lookupMeshDataId( "NonViscTerm_id" ),nonvisc);
	}

	static double getNonViscTerm(pEntity node){
		double nonvisc;
		EN_getDataDbl(node,MD_lookupMeshDataId( "NonViscTerm_id" ),&nonvisc);
		return nonvisc;
	}

	static void setPressure(pEntity node, double p){
		EN_attachDataDbl(node,MD_lookupMeshDataId("p_id"),p);
	}

	static double getPressure(pEntity node){
		double p;
		EN_getDataDbl(node,MD_lookupMeshDataId("p_id"),&p);
		return p;
	}
	static void setSaturation(pEntity node, double sat){
		EN_attachDataDbl(node,MD_lookupMeshDataId("sat_id"),sat);
	}

	static void setSaturation_Old(pEntity node, double sat){
		EN_attachDataDbl(node,MD_lookupMeshDataId("satold_id"),sat);
	}

	static double getSaturation(pEntity node){
		double sat;
		EN_getDataDbl(node,MD_lookupMeshDataId("sat_id"),&sat);
		return sat;
	}

	static double getSaturation_Old(pEntity node){
		double sat;
		EN_getDataDbl(node,MD_lookupMeshDataId("satold_id"),&sat);
		return sat;
	}

	void setInitialSaturation(pMesh, SimulatorParameters*);
	void setSw_max(pEntity,double);
	void setSw_min(pEntity,double);
	void setS_Limit(pEntity,double);

	// GET VELOCITIES
//	void setInitialVelocity(pMesh, SimulatorParameters*);
//	void setVelocity_new(pEdge, const int&, std::vector<double>);
//	void setVelocity_old(pEdge, const int&, std::vector<double>);

	void setVelocity_new(int dom, int row, double* vel){
		velocity[dom].setValue(row,0,vel[0]);
		velocity[dom].setValue(row,1,vel[1]);
		velocity[dom].setValue(row,2,vel[2]);
	}
	void setVelocity_old(int dom, int row, double* vel){
		velocity[dom].setValue(row,3,vel[0]);
		velocity[dom].setValue(row,4,vel[1]);
		velocity[dom].setValue(row,5,vel[2]);
	}

	void getVelocity_new(int dom, int row, double* vel){
		vel[0] = velocity[dom].getValue(row,0);
		vel[1] = velocity[dom].getValue(row,1);
		vel[2] = velocity[dom].getValue(row,2);
	}

	void getVelocity_old(int dom, int row, double* vel){
		vel[0] = velocity[dom].getValue(row,3);
		vel[1] = velocity[dom].getValue(row,4);
		vel[2] = velocity[dom].getValue(row,5);
	}

//	void setVelocity(pEdge, const int&, bool, std::vector<double>);
//	void getVelocity_new(pEdge, const int &, std::vector<double>&);
//	void getVelocity_old(pEdge, const int &, std::vector<double>&);
//	void getVelocity(pEdge, const int &, bool, std::vector<double>&);
//	void getVelocity(GeomData*, SimulatorParameters*, pEdge, const int&, dblarray&, double&, double&);

	double getSw_max(pEntity);
	double getSw_min(pEntity);
	double getS_Limit(pEntity);

	// GET Volume/mobility,fractionalflux
	double getVolume(pEntity, const int&);
	double getWeightedVolume(pEntity);
	double getTotalMobility(pEntity);
	double getFractionalFlux(const double&);
	double getOilFractionalFlux(const double &);
	double get_ksw(const double&);
	double get_kso(const double&);

	void retrieveSwField(pMesh);
	void storeSwField(pMesh);

	// that's for Sw gradient
	bool isProjected(pEntity);
	void setAsProjected(pEntity);
	void setAsNOTProjected(pEntity);

	// that's for velocity
	bool isVelProjected(pEntity);
	void setVelAsProjected(pEntity);
	void setVelAsNOTProjected(pEntity);

	// get set/get pointers to arrays of pointer functions
	GetPFuncGrad* get_getPFuncArray() { return pGetGradArray; }
	//	GetPFuncScalar* get_getPFuncScalarArray() const { return pGetScalarArray; }
	//	SetPFuncScalar* get_setPFuncScalarArray() const { return pSetScalarArray; }

	/*! brief For steady-state simulations, total mobility must be equal 1, otherwise it must be calculated.
	 * \param state if true, lambda_total = 1.0;
	 */
	void setSimulationState(bool state){
		steady_state = state;
	}

private:
	bool steady_state;	// check if steady-state or transient simulations
	double Swr;			// Irreducible water saturation
	double Sor;			// Residual oil saturation
	double mi_w;		// water viscosity
	double mi_o;		// oil viscosity
	int ksModel;		// rel. permeability model flag

	// pointers to arrays of pointer functions
	GetPFuncGrad* pGetGradArray;

	Matrix<double>* velocity;
};
}
#endif /*PHYSICALPROPERTIESDATA_H_*/
