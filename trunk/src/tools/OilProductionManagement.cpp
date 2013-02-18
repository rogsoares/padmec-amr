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
		ttpOilProduction = false;
		printStepSum = getPrintStep();
		output_frequency = getPrintStep();
		PVI_accumulated = .0;

		PetscPrintf(PETSC_COMM_WORLD,"IOV = %f   TIFR = %f\n",IOV,TIFR);

		/*
		 * Open a file and write on it oil recovered and accumulated.
		 * If restart is required, new data will be wrote from the end of file.
		 */
		fid.open(fname.c_str(), ios::ate);

		if (!fid.is_open()) {
			char msg[256]; sprintf(msg,"File '%s' could not be opened or it does not exist.\n",fname.c_str());
			throw Exception(__LINE__,__FILE__,msg);
		}
		fid << "PVI Oil-Recovered Oil-Accumulated" << std::endl;
		fid.precision(10);
	}

	OilProductionManagement::~OilProductionManagement(){
		fid << "CLOSE" << std::endl;
		fid.close();
	}

	void OilProductionManagement::printOilProduction(double timeStep, double accumlatedTime,
			double totalSimulationTime,	double recOilValue){
		double oilRec = recOilValue;
		static double oilAcc = .0;
		oilAcc += (double)oilRec*timeStep/IOV;

		if (oilRec<1.0 && output_frequency == 0.1){
			printf("WARNING: oilRec = %f and oilAcc = %f at ts = %f\n",oilRec,oilAcc,accumlatedTime);
		}

		if (accumlatedTime >= output_frequency){
			/*
			 * Print data on file using commas in place of dots instead.
			 */
			char cString[256];
			sprintf(cString,"%.10f %.10f %.10f",(double)accumlatedTime/totalSimulationTime
			                               ,(double)oilRec/TIFR
			                               ,oilAcc);
			string theString(cString);
			replaceAllOccurencesOnString(theString,1,".",",");
			fid << theString << endl;

			// a new VTK file will be printed at each 0.01 PVI (1% of total simulation time)
			PVI_accumulated += getPrintStep();

			// when next vtk file must be print out
			output_frequency = PVI_accumulated*totalSimulationTime;
		}
	}
}
