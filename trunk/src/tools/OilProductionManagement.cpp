/*
 * OilProductionManagement.cpp
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#include "OilProductionManagement.h"

namespace PRS{

	OilProductionManagement::OilProductionManagement(){}

	OilProductionManagement::OilProductionManagement(string fname, double iov,double TotalInjectionFlowRate){
		IOV = iov;
		TIFR = TotalInjectionFlowRate;
		PetscPrintf(PETSC_COMM_WORLD,"IOV = %f   TIFR = %f\n",IOV,TIFR);
		fid.open(fname.c_str(), ios::ate);
		if (!fid.is_open()) {
			char msg[256]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n",fname.c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}
		fid << "PVI time-step cumulative-time Recovered-Oil Cumulative-Oil Num.Time-steps" << std::endl;
	}

	OilProductionManagement::~OilProductionManagement(){
		fid << "CLOSE" << std::endl;
		fid.close();
	}

	void OilProductionManagement::printOilProduction(double timeStep,
			                                         double cml_time,
			                                         double total_SimTime,
			                                         double rec_oil,
			                                         double cml_oil,
			                                         int timestep_counter){
		fid << fixed << setprecision(2) << (double)(cml_time/total_SimTime) << " " << setprecision(7)
			<< timeStep << " " << cml_time << " " << rec_oil << " " <<  cml_oil/IOV << " " << timestep_counter << endl;
	}
}
