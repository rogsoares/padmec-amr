/*
 * EBFV1_modified.cpp
 *
 *  Created on: Oct 17, 2014
 *      Author: rogerio
 */

#include "EBFV1_modified.h"

void EBFV1_modified_preprocessor_2D(pMesh theMesh, GeomData* pGCData){

	// loop over triangles elements. For each one calculate:
	// Geometric coeficients: Cij, Dij, and nodal volumes;
	// matrices E, F and G;
	// for all three edges and the assembly the local edge matrix

	pEntity face;
	double Dij[6], Cij[6], vol;
	double E_ij[8], E_ik[8], E_jk[8];
	double F_ij[8], F_ik[8], F_jk[8];
	double G_ij[4], G_ik[4], G_jk[4];
	double A_ij[6], A_jk[6], A_ik[6];

	pGCData->initializeElementMatrix(M_numFaces(theMesh));

	int k = 0;
	FIter fit = M_faceIter(theMesh);
	while ( (face = FIter_next(fit)) ){
		calculateGeomCoefficients(face,Cij,Dij,vol);
		calculateMatrix_E(face,pGCData,Cij,E_ij,E_jk,E_ik);
		calculateMatrix_F(face,Cij,Dij,vol,F_ij,F_jk,F_ik);
		calculateMatrix_G(face,pGCData,Cij,G_ij,G_jk,G_ik);
		assemblyMatrix_A(E_ij,E_jk,E_ik,F_ij,F_jk,F_ik,G_ij,G_jk,G_ik,A_ij,A_jk,A_ik);
		pGCData->setElementMatrices(k,A_ij,A_jk,A_ik);
		pGCData->setIJK_IDs(k,EN_id(face->get(0,0)),EN_id(face->get(0,1)),EN_id(face->get(0,2)));
	}
	FIter_delete(fit);
}

void getPermeabilityTensor(pEntity elem, double* K){

}
