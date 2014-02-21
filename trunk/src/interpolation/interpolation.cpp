/*
 * interpolation.cpp
 *
 *  Created on: Dec 6, 2013
 *      Author: rogerio
 */

#include "interpolation.h"

void interpolation(InterpolationDataStruct* pIData, INTERPOLATION_OPTIONS opt){
	initialize(pIData);
	int dim = pIData->m1->getDim();
	switch ( opt ){
	case h_REFINEMENT:
		throw Exception(__LINE__,__FILE__,"Underconstruction!");
		break;
	case LINEAR:
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tIN\n";
		#endif
		calculate_GeometricCoefficients(pIData,dim);
		calculate_LinearInterpolation(pIData,dim);
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tOUT\n";
		#endif
		break;
	case QUADRATIC:
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
		#endif
		calculate_GeometricCoefficients(pIData,dim);
		calculate_LinearInterpolation(pIData,dim);
		for (int field=0; field<pIData->numFields; field++){
			calculate_Gradients(pIData,field);
			calculate_DerivativesError(pIData);
			calculate_QuadraticInterpolation(pIData,field);
		}
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
		#endif
		break;
	default:
		throw Exception(__LINE__,__FILE__,"Interpolation method unknown. Exiting....");
	}
	finalize(pIData);
}


void initialize(InterpolationDataStruct* pIData){
	int dim = pIData->m1->getDim();
	pIData->theOctree = OctreeCreate2<iterall>(pIData->m2->beginall(dim),pIData->m2->endall(dim),dim);
	int nrows1 = M_numVertices(pIData->m1) + 1;	// New mesh (to)
	int nrows2 = M_numVertices(pIData->m2) + 1; // Old mesh (from)
	pIData->pNodeValue.allocateMemory(nrows2);
	pIData->pGrad.allocateMemory(nrows2,5);
	pIData->pGeomCoeff.allocateMemory(nrows1,dim+1);
	pIData->pInterpolatedVals.allocateMemory(nrows1);
}

void finalize(InterpolationDataStruct* pIData){
	pIData->pNodeValue.freeMemory();
	pIData->pInterpolatedVals.freeMemory();
	pIData->pGeomCoeff.freeMemory();
	pIData->pGrad.freeMemory();
	pIData->theOctree = 0;
}
