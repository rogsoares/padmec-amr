/*
 * EBFV1_modified.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */

#include "EBFV1_modified.h"

void EBFV1_modified_preprocessor_2D(pMesh theMesh, GeomData* pGCData){

	// loop over triangles elements. For each one calculate:
	// 		Geometric coeficients: Cij, Dij, and nodal volumes;
	// 		matrices E, F and G;

	pEntity face;
	int numElements, dim;
	double K[4];
	double Dij[6], Cij[6], vol;
	double E_ij[8], E_jk[8], E_ik[8];	// E: 2x4
	double F_ij[8], F_jk[8], F_ik[8];	// F: 4x2
	double G_ij[4], G_jk[4], G_ik[4];	// G: 2x2
	double A_ij[6], A_jk[6], A_ik[6];	// A: 2x3
	RockPropFuncPointer* pFUNC_RockProp;

	dim = theMesh->getDim();
	numElements = M_numFaces(theMesh);
	pGCData->initializeElementMatrix(numElements);
	calculateEdgeLength(theMesh,pGCData);
	getRockPropertyFuncPointer(pFUNC_RockProp);					// based on user inputs, define a function pointer to what rock property

	int k = 0;
	FIter fit = M_faceIter(theMesh);
	while ( (face = FIter_next(fit)) ){
		pRockPropFunc(face,dim,K)
		calculateGeomCoefficients(face,Cij,Dij,vol);			// calculates Cij and Dij for all three element edges and volume as well.
		calculateMatrix_E(face,K,pGCData,Cij,E_ij,E_jk,E_ik);		// calculates edge matrix for all three edges
		calculateMatrix_F(face,Cij,Dij,vol,F_ij,F_jk,F_ik);		// calculates edge matrix for all three edges
		calculateMatrix_G(face,K,pGCData,Cij,G_ij,G_jk,G_ik);		// calculates edge matrix for all three edges
		assemblyMatrix_A(E_ij,E_jk,E_ik,F_ij,F_jk,F_ik,G_ij,G_jk,G_ik,A_ij,A_jk,A_ik);
		pGCData->setElementMatrices(k,A_ij,A_jk,A_ik);
	}
	FIter_delete(fit);
}
