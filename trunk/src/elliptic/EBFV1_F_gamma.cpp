#include "EBFV1_elliptic.h"

namespace PRS{

	int EBFV1_elliptic::gradient_F_bdry(Mat F, int dom){
		double Dij[3], Fij_column1[4], Fij_column2[4], volumeI, volumeJ, volumeK;
		int ndom, nedges, nfaces, id0, id1, id2, idx_0, idx_1, i, j, pos1, pos2, pos3;

		if (pGCData->getMeshDim()==2){
			nedges = pGCData->getNumBDRYEdgesPerDomain(dom);
			for (j = 0; j<nedges; j++){
				pGCData->getBdryEdge(dom,j,idx_0,idx_1);
				pGCData->getBdryID(dom,idx_0,idx_1,id0,id1);
				pGCData->getBdryVolume(dom,idx_0,idx_1,volumeI,volumeJ);
				pGCData->getDij(dom,j,Dij);
				Fij_column1 [0] = 5.*Dij[0]/(6.*volumeI);
				Fij_column1 [1] = 5.*Dij[1]/(6.*volumeI);
				Fij_column1 [2] = Dij[0]/(6.*volumeJ);
				Fij_column1 [3] = Dij[1]/(6.*volumeJ);
				Fij_column2 [0] = Dij[0]/(6.*volumeI);
				Fij_column2 [1] = Dij[1]/(6.*volumeI);
				Fij_column2 [2] = 5.*Dij[0]/(6.*volumeJ);
				Fij_column2 [3] = 5.*Dij[1]/(6.*volumeJ);
				id0 = pMData->get_AppToPETSc_Ordering(id0);
				id1 = pMData->get_AppToPETSc_Ordering(id1);
				pos1 = 2*(id0-1);
				pos2 = 2*(id1-1);
				int idxm[4] = {pos1,pos1+1,pos2,pos2+1};
				int idxn[2] = {id0-1,id1-1};
				ierr = MatSetValues(F,4,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
				ierr = MatSetValues(F,4,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
			}
		}
		else{
			double aux_vec[3];
			//nfaces = pGCData->getNumBDRYFacesPerDomain(dom);
			for (j = 0; j<nfaces; j++){
				//pGCData->getBdryVolume(dom,idx_0,idx_1,idx2,aux_vec,3);
				for (i=0; i<3; i++){
					aux_vec[i] = 1./(8.*aux_vec[i]);
				}

				double Fij_1[9] = {6.*aux_vec[0]*Dij[0], 6.*aux_vec[0]*Dij[1], 6.*aux_vec[0]*Dij[2],
						aux_vec[0]*Dij[0],    aux_vec[0]*Dij[1],    aux_vec[0]*Dij[2],
						aux_vec[0]*Dij[0],    aux_vec[0]*Dij[1],    aux_vec[0]*Dij[2]};

				double Fij_2[9] = {   aux_vec[1]*Dij[0],    aux_vec[1]*Dij[1],    aux_vec[1]*Dij[2],
						6.*aux_vec[1]*Dij[0], 6.*aux_vec[1]*Dij[1], 6.*aux_vec[1]*Dij[2],
						aux_vec[1]*Dij[0],    aux_vec[1]*Dij[1],    aux_vec[1]*Dij[2]};

				double Fij_3[9] = {   aux_vec[2]*Dij[0],    aux_vec[2]*Dij[1],    aux_vec[2]*Dij[2],
						aux_vec[2]*Dij[0],    aux_vec[2]*Dij[1],    aux_vec[2]*Dij[2],
						6.*aux_vec[2]*Dij[0], 6.*aux_vec[2]*Dij[1], 6.*aux_vec[2]*Dij[2]};

				pos1 = 3*(pMData->get_AppToPETSc_Ordering(id0) - 1);
				pos2 = 3*(pMData->get_AppToPETSc_Ordering(id1) - 1);
				pos3 = 3*(pMData->get_AppToPETSc_Ordering(id2) - 1);

				int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
				int idxn[3] = {id0-1, id1-1, id2-1};

				ierr = MatSetValues(F,9,idxm,1,&idxn[0],Fij_1,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValues(F,9,idxm,1,&idxn[1],Fij_2,ADD_VALUES);CHKERRQ(ierr);
				ierr = MatSetValues(F,9,idxm,1,&idxn[2],Fij_3,ADD_VALUES);CHKERRQ(ierr);
			}
		}
		return 0;
	}
}
