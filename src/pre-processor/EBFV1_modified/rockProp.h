/*
 * rockProp.h
 *
 *  Created on: Oct 31, 2014
 *      Author: rogerio
 */

#ifndef ROCKPROP_H_
#define ROCKPROP_H_

#include "auxiliar.h"
//This file is intent to define functions which represent rock properties such as:
//
//		1) absolute permeability tensor
//		2) porosity

// RPFP: Rock property function pointer
typedef void (*RockPropFuncPointer)(pEntity,int,double*);

void getRockPropertyFuncPointer(RockPropFuncPointer*);


/*
 * USERS MUST DEFINE THEIR FUNCTION BELOW
 */

void homogeneousPermeability(pEntity,int,double*);

#endif /* ROCKPROP_H_ */
