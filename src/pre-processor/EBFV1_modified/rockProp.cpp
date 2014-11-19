/*
 * rockProperties.cpp
 *
 *  Created on: Oct 26, 2014
 *      Author: rogerio
 */

#include "rockProp.h"

void getRockPropertyFuncPointer(RPFP* pFunc){

#ifdef HOMOGENEOUSPERMEABILITY
	pFunc = homogeneousPermeability;
#endif

}

void homogeneousPermeability(pEntity face, int dim, double* K){
	K[0] = 1.0; K[1] = 0.0;
	K[2] = 0.0; K[3] = 1.0;
}

