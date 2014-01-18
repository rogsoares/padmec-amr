#ifndef PHYSICALPROPERTIESDATA_H_
#define PHYSICALPROPERTIESDATA_H_

#include "MeshData.h"
#include "Matrix.h"
#include "GeomData.h"

extern Matrix<double> *pGrad_dom;
extern Matrix<double> pressure;
extern Matrix<double> SwGrad;
extern Matrix<double> Sw;


namespace PRS{

/**
 * PhysicPropData class provides methods to set/get physical properties such as pressure, saturation and their respective gradients
 * as well as rock-fluid properties (fractional flow and total mobilities). Fluid velocity is also managed by this class. All proper-
 * ties are associated to a mesh entity (node or edge) through the class AttachData that is in charge of store all data (physical,
 * numeric, geometric,...) to mesh entities.
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
	static void set_pw_Grad(int dom, int row, const double* grad){
		pGrad_dom[dom].setValue(row,0,grad[0]);
		pGrad_dom[dom].setValue(row,1,grad[1]);
		pGrad_dom[dom].setValue(row,2,grad[2]);
	}

	static void get_pw_Grad(int dom, int row, double* grad){
		grad[0] = pGrad_dom[dom].getValue(row,0);
		grad[1] = pGrad_dom[dom].getValue(row,1);
		grad[2] = pGrad_dom[dom].getValue(row,2);
	}

	static void setPressure(int idx, double p){
		pressure.setValue(idx,p);
	}

	static double getPressure(int idx, double& p){
		p = pressure.getValue(idx);
		return p;
	}

	static void set_Sw_Grad(pVertex node, int dom, int row, const double* grad){
		for(int i=0;i<3;i++){
			char string[256]; sprintf(string,"Swgrad_%d_%d",dom,i);
			EN_attachDataDbl(node,MD_lookupMeshDataId(string),grad[i]);
		}
	}

	static void get_Sw_Grad(pVertex node, int dom, int row, double* grad){
		for(int i=0;i<3;i++){
			char string[256]; sprintf(string,"Swgrad_%d_%d",dom,i);
			EN_getDataDbl(node,MD_lookupMeshDataId(string),&grad[i]);
		}
	}

	static void set_Sw_Grad(int row, const double* grad){
		SwGrad.setValue(row,0,grad[0]);
		SwGrad.setValue(row,1,grad[1]);
		SwGrad.setValue(row,2,grad[2]);
	}

	static void get_Sw_Grad(int row, double* grad){
		grad[0] = SwGrad.getValue(row,0);
		grad[1] = SwGrad.getValue(row,1);
		grad[2] = SwGrad.getValue(row,2);
	}

	void set_Sw_Grad(int dom, int row, const double* grad){
		SwGrad_dom[dom].setValue(row,0,grad[0]);
		SwGrad_dom[dom].setValue(row,1,grad[1]);
		SwGrad_dom[dom].setValue(row,2,grad[2]);
	}

	void get_Sw_Grad(int dom, int row, double* grad){
		grad[0] = SwGrad_dom[dom].getValue(row,0);
		grad[1] = SwGrad_dom[dom].getValue(row,1);
		grad[2] = SwGrad_dom[dom].getValue(row,2);
	}

	static void setSaturation(int idx, double v){
		Sw.setValue(idx,v);
	}

	static double getSaturation(int idx){
		return Sw.getValue(idx);
	}

	static void setNonViscTerm(pEntity node, double nonvisc){
		EN_attachDataDbl(node,MD_lookupMeshDataId( "NonViscTerm_id" ),nonvisc);
	}

	static double getNonViscTerm(pEntity node){
		double nonvisc;
		EN_getDataDbl(node,MD_lookupMeshDataId( "NonViscTerm_id" ),&nonvisc);
		return nonvisc;
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

//	double getSw_max(pEntity);
//	double getSw_min(pEntity);
//	double getS_Limit(pEntity);

	// GET Volume/mobility,fractionalflux
	double getVolume(pEntity, const int&);
	double getWeightedVolume(pEntity);
	double getTotalMobility(pEntity);
	double getTotalMobility(double Sw);
	double getFractionalFlux(const double&);
	double getOilFractionalFlux(const double &);
	double get_ksw(const double&);
	double get_kso(const double&);

	void retrieveSwField(pMesh);
	void storeSwField(pMesh);


	// get set/get pointers to arrays of pointer functions
	GetPFuncGrad* get_getPFuncArray() { return pGetGradArray; }

	/*! brief For steady-state simulations, total mobility must be equal 1, otherwise it must be calculated.
	 * \param state if true, lambda_total = 1.0;
	 */
	void setSimulationState(bool state){
		steady_state = state;
	}

	int getNumNodeFreeWells() const{
		return nfree;
	}

	int getNumNodesWells() const{
		return nneumann;
	}

	void getNeumannIndex(int i, int &idx) const{
		idx = pWellsNeumann_index[i];
	}

	void getFreeIndex(int i, int &idx) const{
		idx = pWellsFree_index[i];
	}

	void getNonvisc(int idx, double &val){
		val = nonvisc.getValue(idx);
	}

	void setNonvisc(int idx, double val){
		nonvisc.setValue(idx,val);
	}

	void resetNonvisc(double&);

	void initializeNonvisc(){
		nonvisc.initialize(.0);
	}

	bool isInjectionWell(int idx){
		return injectionWell.getValue(idx);
	}

	void setProjectedSwgrad(int k, bool val){
		projectedSw_grad.setValue(k,val);
	}

	bool getProjectedSwgrad(int k){
		return projectedSw_grad.getValue(k);
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
	Matrix<double> nonvisc;
	Matrix<double> *SwGrad_dom;
	Matrix<bool> projectedSw_grad;
	Matrix<bool> injectionWell; // idx = 0 (not injectio/ could be production), idx = 1 free well
	int nfree;
	int nneumann;
	int* pWellsFree_index;
	int* pWellsNeumann_index;

	int nnodes;
};
}
#endif /*PHYSICALPROPERTIESDATA_H_*/
