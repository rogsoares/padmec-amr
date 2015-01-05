/*
 * utilities.h
 *
 *  Created on: Jan 4, 2015
 *      Author: rogerio
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "include.h"

double R_Volume(const double* p1, const double* p2, const double* p3, const double* p4);
void computeCrossProduct(const double* a, const double* b, double* n);
double computeDotProduct(const double* a, const double* b);


#endif /* UTILITIES_H_ */
