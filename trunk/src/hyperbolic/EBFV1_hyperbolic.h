/*
 * EBFV1-hyperbolic.h
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#ifndef EBFV1HYPERBOLI_H_
#define EBFV1HYPERBOLI_H_

#include "Hyperbolic_equation.h"
#include "HighOrderApproximation.h"

namespace PRS
{

	/**
	 * EBFV1_hyperbolic class holds all steps needed to solve the saturation e-
	 * quation including the velocity equation. The main functions implemented
	 * are:
	 * 		calculateVelocityField (fluid velocity associated to edges)
	 * 		calculateIntegralAdvectiveTerm (flux through a control volume)
	 * 		calculateExplicitAdvanceInTime (performe: S(n+1) = S(n) + advanceTime
	 */
	class EBFV1_hyperbolic : public Hyperbolic_equation{
	public:

		EBFV1_hyperbolic();
		EBFV1_hyperbolic(pMesh, PhysicPropData*, SimulatorParameters*,GeomData*, MeshData*, OilProductionManagement*, ErrorAnalysis*);
		~EBFV1_hyperbolic();
		virtual double solver(pMesh, double&);

	protected:

		/**
		 * Main functions related to the advective equation presented in order
		 * that they must be called.
		 */
		double calculateVelocityField(pMesh, int,int);
		double calculateIntegralAdvectiveTerm(pMesh, int, int, double&);
		double calculateExplicitAdvanceInTime(pMesh, double);

		/*
		 * For new time step, erase previous nonvisc data
		 */
		void resetNodalNonviscTerms(pMesh);

		/*
		 * For parallel simulation only. Node on partition boundary have the same value
		 */
		int updateNonViscTerm(pMesh);

		double getTimeStep();

		/*
		 * Evaluate Sw_new = Sw_old + DT*nonvisc on node without and with production wells respectively
		 */
		int nodeWithOut_Wells(pMesh, double);
		void nodeWith_Wells(pMesh, double);
		
		/// For mesh adaptation, saturation field interpolation between old and new mesh over Sw_t, and not Sw_t+1
		void saveSwField(pMesh);

		void setRecoveredOil(double val){
			oilRecovered = val;
		}

		double getCumulativeOil() const{
			return _cumulativeOil;
		}

		void setCumulativeOil(double co){
			_cumulativeOil = co;
		}

		double getRecoveredOil() const{
			return oilRecovered;
		}

		double oilRecovered;

		/**
		 * Numeric parameter to compute nonvisc term
		 */
		double alpha_max;

		/**
		 * For multi domains problems, one time step is compute by domain and
		 * the minimum is used to advance saturation
		 */
		std::set<double> timeStepByDomain;

		OilProductionManagement *pOPManager;
		MeshData *pMData;
		GeomData *pGCData;

		/**
		 * Define a pointer to handle high order approximations for saturation
		 * field.
		 */
		HighOrderApproximation *pHOApproximation;

		/**
		 * Saturation gradient is calculated for all domains at once.
		 * Nodes on boundary domains contains one gradient vector for each
		 * domain.
		 */
		double calculateSaturationGradient(pMesh);

//		int *rowToImport;
//		Mat joinNodes;
//		Mat updateValues;
//		map<int,double> mapPB_nodes;	// map partition boundary nodes


	private:
		/// every new time-step, nodal gradient must be set to zero and start a new calculation
		void resetSaturationGradient(pMesh);
		
		/// loop over all edges (omega domain
		void calc_Sw_grad_1(pMesh, int, int, int);
		
		/// loop over all boundary edges (2-D) or all external faces (3-D)
		void calc_Sw_grad_2(pMesh, int, int);
		
		/// Averaging by Total Volume in 3-D (area in 2-D problems)
		void calc_Sw_grad_3(pMesh, int, int);
		
		void calc_Sw_grad_31(pMesh);

		/// Imposition of Homogeneus Neumman Boundary Conditions
		void calc_Sw_grad_4(pMesh, int);
		
		double _cumulativeOil;

		PetscErrorCode ierr;
	};
}


#endif /* EBFV1HYPERBOLI_H_ */
