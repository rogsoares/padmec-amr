/*
 * EBFV1_preprocessor.cpp
 *
 *  Created on: 17/08/2012
 *      Author: rogsoares
 */

#include "EBFV1__pre-processors.h"

int EBFV1_preprocessor(pMesh theMesh, void *pData){
	int ndom;
	if (theMesh->getDim()==2){
		EBFV1_preprocessor_2D(theMesh,pData,ndom);
	}
	else{
		EBFV1_preprocessor_3D(theMesh,pData,ndom);
	}
	return 0;
}
