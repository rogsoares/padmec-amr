/*
 * CalculateGlobalError.cpp
 *
 *  Created on: 28/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

void ErrorAnalysis::calculate_GlobalError(pMesh theMesh, GetPFuncGrad pGetGradient){
	double error_sum = this->calculate_ErrorSum(theMesh);
	double smooth_gradNorm = this->getSmoothedGradNorm();
	setGlobalError(sqrt(fabs(error_sum/smooth_gradNorm)));
}

/*
 * Define an average error. It's an distributed error over the mesh.
 * 	//eq. 4.24 e 3.2.1 (p‡g 32) do algoritmo - tese de Filipe
 */
void ErrorAnalysis::calculate_AvgError(pMesh theMesh, double tol, bool excludingSingularities){
	double SGN = (excludingSingularities)?this->getSmoothedGradNorm_singular():this->getSmoothedGradNorm();
	//this->countElements(theMesh,excludingSingularities);
	double numElements = (double)(excludingSingularities)?this->getNumElements_excludingSingularities():this->getNumElements();
	double avgError = tol*( sqrt(SGN/numElements) );
	if (excludingSingularities){
		setAverageError_excludingSingularities(avgError);
	}
	else{
		setAverageError(avgError);
	}
}
