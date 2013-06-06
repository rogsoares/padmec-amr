/*
 * IAR_main.cpp
 *
 *  Created on: 14/02/2013
 *      Author: rogsoares
 *
 *      IAR - Interpolation Adaptative Remeshing
 *      Here are called all interpolation function for remeshing procedure
 */


#include "interpolation.h"

void Linear(InterpolationDataStruct* pIntpData){
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tIN\n";
#endif

	int dim = pIntpData->m1->getDim();
	/*
	 * Octree structure is created for searching the element on base mesh which contains data fields to be interpolated to the adaptaded mesh.
	 *
	 * pIntpData->m2:	Base mesh
	 * pIntpData->m1:	To receive interpolated data from pIntpData->m2
	 */
	pIntpData->theOctree = OctreeCreate2<iterall>(pIntpData->m2->beginall(dim),pIntpData->m2->endall(dim),dim);


	calculate_GeometricCoefficients(pIntpData,dim);
	calculate_LinearInterpolation(pIntpData,dim);

#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tOUT\n";
#endif
}

void Quadratic(InterpolationDataStruct* pIntpData){
//	int dim = pIntpData->m1->getDim();
//	Linear(pIntpData);
//	calculate_Gradients(pIntpData);
//	calculate_DerivativesError(pIntpData);
//	calculate_QuadraticInterpolation(pIntpData);
}

void Adaptative(InterpolationDataStruct* pIntpData){
	throw Exception(__LINE__,__FILE__,"Under construction!");
}

void Conservative(InterpolationDataStruct* pIntpData){
	throw Exception(__LINE__,__FILE__,"Under construction!");
}

void PureInjection(InterpolationDataStruct* pIntpData){
	throw Exception(__LINE__,__FILE__,"Under construction!");
}

void HalfWeighting(InterpolationDataStruct* pIntpData){
	throw Exception(__LINE__,__FILE__,"Under construction!");
}

void FullWighting(InterpolationDataStruct* pIntpData){
	throw Exception(__LINE__,__FILE__,"Under construction!");
}

