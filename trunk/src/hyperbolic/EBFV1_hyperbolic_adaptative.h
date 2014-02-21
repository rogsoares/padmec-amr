/*
 * EBFV1_hyperbolic_adaptative.h
 *
 *  Created on: 05/05/2009
 *      Author: rogerio
 */

#ifndef EBFV1_HYPERBOLIC_ADAPTATIVE_H_
#define EBFV1_HYPERBOLIC_ADAPTATIVE_H_

#include "EBFV1_hyperbolic.h"

namespace PRS{
	/**
	 * class: EBFV1_hyperbolic_adaptative
	 *
	 * The goal of this class is to implement an adaptative saturation time advance
	 * based on the velocity norm variation. It's supposed that the velocity field varies
	 * much slower than the saturations field. It suggests that the velocity variable
	 * need not be calculated every time step. Thus, a significantly reduce in CPU ti-
	 * me consumption may be reached without compromise the results accuracy.
	 * It's defined two time steps counters: one related to the velocities field
	 * (vel_ts) and another one related to the saturations field (sat_ts). It's assumed
	 * that vel_ts >= sat_ts and sat_ts is restricted by the CFL condition. The variable
	 * vel_ts is always calculated first than sat_ts. It's assumed an extrapolation in the
	 * velocities filed variation behavior in such way results will not be prejudiced.
	 *
	 * Here, functions defined in EBFV1_hyperbolic class are reused.
	 */
	class EBFV1_hyperbolic_adaptative : public EBFV1_hyperbolic {
	public:

		EBFV1_hyperbolic_adaptative();
		EBFV1_hyperbolic_adaptative(pMesh, PhysicPropData *, SimulatorParameters *,GeomData *, MeshData *, OilProductionManagement *, ErrorAnalysis*);
		~EBFV1_hyperbolic_adaptative();
		double solver(pMesh, double&);

	private:

		int sumNum_p_DT;						// summation of number of pressure timesteps
		int sumNum_Sw_DT;						//  summation of number of Sw timesteps
		double p_timestepOld;
		double cumulative_p_timestep;

		int vel_ts_counter;						/// counts number of implicit time-steps
		int sat_ts_counter;						/// counts number of explicit time-steps
		double DVTOL;							/// velocity variation norm tolerance
		double vel_factor;						/// turn velocity vector dimensionless
		double time_factor;						/// turn sat_ts dimensionless
		double sat_ts;							/// explicit time advance
		double vel_ts;							/// implicit time advance
		double vel_ts_old;						/// last implicit time advance
		bool initial_vel_ts_old;				/// set vel_ts_old as sat_ts at the beginning of simulation

		bool   allowAdaptativeTimeAdvance;		/// controls saturations field time advance
		bool   allowVelocityCalculation;				/// make velocities filed be calculated
		/// just once per explicit (sat_ts) time-step

		double cumulativeExplicitTS;			/// cumulativeExplicitTS = sum(sat_ts)
		/// if (cumulativeExplicitTS == vel_ts)
		/// calculate new pressure field
		//double cumulativeSimulationTime;

		// implicitTS returns a dimensionsless delta T
		double dv_norm;

		/**
		 * Print on file adaptative data:
		 * 	- time (PVI)
		 *  - explicit_counter/implicit_counter ratio
		 *  - velocity norm
		 */
		ofstream fid;
		bool print_ATE_Data;
		void MIMPES_output(double Sw_timestep_sum, double p_timestep, int numCFL_Steps, int sumNum_p_DT, int sumNum_Sw_DT, double DV_norm);


		/**
		 * Compute a new implicit time-step (vel_ts) if necessary. This happens
		 * either at the beginning of simulation or when sum of sat_ts be greater
		 * than DT.
		 */
		void calculateImplicitTS(double &p_timestep, double timeStep, double &DV_norm, bool &go_ImplicitTS, bool &go_MIMPES);

		/**
		 * Compute a new implicit time-step as function of the last implicit time-step,
		 * velocity norm and DVTOL. DVTOL comes from the input data file 'numeric.dat'.
		 */
		void calculateVelocityVariationNorm(double &DV_norm, int dim);

		/**
		 * reset logical variables which control programming flux.
		 */
		void resetVariables();

		/**
		 * time_factor and vel_factor are used to turn implicit time-step and velocity vector
		 * dimensionless variables
		 */
		void createDimensionLessFactors();

		/**
		 * Verify if the explicit time advance (summation of sat_ts) do not exceeded the
		 * limit imposed by the implicit time advance (vel_ts)
		 */
		void correct_p_TS(double &p_timestep, bool &go_MIMPES);
		void correct_Sw_TS(double &timeStep, bool &go_MIMPES);

		void setCumulative_p_TS(double ts){
			cumulative_p_timestep = ts;
		}
		double getCumulative_p_TS() const{
			return cumulative_p_timestep;
		}
	};
}
#endif /* EBFV1_HYPERBOLIC_ADAPTATIVE_H_ */
