/*
 * Restart.h
 *
 *  Created on: 22/07/2009
 *      Author: rogerio
 */

#ifndef RESTART_H_
#define RESTART_H_

#include "PhysicPropData.h"

namespace PRS{

/*
 * RestartSimulation allows to continuing simulation from the last vtk file ge-
 * nerated. No matter how many time steps have been performed since that point,
 * the simulation will continue exactly from the last vtk file, which saturation
 * field will be used as initial condition. Remember that each vtk file corres-
 * ponding to 5% of all simulation time.
 */
	void restartSimulation(SimulatorParameters *, PhysicPropData *, pMesh);

	}
#endif /* RESTART_H_ */

/**
 * Restart to do:
 *
 * 1 - vtk files are continued for serial
 * 2 - vtk files are continued for parallel
 * 3 - oil production are continued for serial
 * 4 - oil productionare are continued for parallel
 * 5 - s_AdaptTimeData-DVTOL are continued for serial
 * 6 - s_AdaptTimeData-DVTOL are continued for parallel
 */
