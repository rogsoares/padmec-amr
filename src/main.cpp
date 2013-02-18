#include "SIMULATION_core.h"

int main(int argc, char **argv){
	/*@
	 * SIMULATION_core managers all simulation steps.
	 * Simulation is divided into three steps:
	 *   	initialization    : pointer creation, data reading, variables defining
	 *    	time-step solving : call equation functions
	 *    	stop simulation   : write output data to files, free memory
	 *
	 * Note.: During the execution of any of these three steps an error signal
	 * 		  (Exception) can be thrown. If it happens, the exception is caught
	 *        and simulation terminates after showing the cause and where it
	 *        occurred.
	 @*/
	PRS::SIMULATION_core sim;
	try{
		sim.initialize(argc,argv);
		sim.solver();
	}
	catch (Exception excp) {
		excp.showExceptionMessage();
	}

	/*
	 * If an exception is thrown, finalize all necessary things before terminate running.
	 */
	sim.finalize();
	return 0;
}
