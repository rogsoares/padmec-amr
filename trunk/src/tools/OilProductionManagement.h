/*
 * OilProductionManagement.h
 *
 *  Created on: 14/01/2009
 *      Author: rogerio
 */

#ifndef OILPRODUCTIONMANAGEMENT_H_
#define OILPRODUCTIONMANAGEMENT_H_

#include "auxiliar.h"

namespace PRS{

	struct OPData{
		double oilRec;
		double oilAcc;
		double time_step;
	};

	typedef std::list<OPData> OPList;
	typedef std::list<OPData>::iterator LIter;

	class OilProductionManagement{
	public:

		OilProductionManagement();
		OilProductionManagement(string, double, double);
		~OilProductionManagement();

		void update_OilProductionData(double,double,double);

		LIter OP_history_begin();
		LIter OP_history_end();

		void createOilProductionOutput(string);

		// print a new oil production value frequency controlled by getPrintStep
		void printOilProduction(double timeStep,
                double cml_time,
                double total_SimTime,
                double rec_oil,
                double cml_oil,
                int timestep_counter);

		// print frequency control
		double getPrintStep() const { return (double)1./200.; }

	private:

		/*Stores oil recovery for each time step*/
		OPList OP_history;

		double output_frequency;
		double PVI_accumulated;

		double printStepSum;
		bool ttpOilProduction;	// time to print Oil production
		ofstream fid;			// output stream for oil production
		double IOV;				// Initial Oil Volume
		double TIFR;			// Total Injection Flow Rate
	};
}

#endif /* OILPRODUCTIONMANAGEMENT_H_ */
