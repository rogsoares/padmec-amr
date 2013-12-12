#include "EBFV1_elliptic.h"

namespace PRS
{
int EBFV1_elliptic::gradient_F_bdry(pMesh mesh, Mat F, const int &dom){
	switch ( pGCData->getMeshDim() ){
	case 2:	F_bdryEdges(mesh,F,dom); break;
	case 3: F_bdryFaces(mesh,F,dom); break;
	}
	return 0;
}


// For 2-D domains: get boundary edges contribution.
int EBFV1_elliptic::F_bdryEdges(pMesh mesh, Mat F, const int &dom){
	dblarray Dij(3,.0);
	double Fij_column1[4], Fij_column2[4];

	int numBE1 = 0;
	double sumx = .0;
	double sumy = .0;

	EIter eit = M_edgeIter(mesh);
	while (pEntity edge = EIter_next(eit)){
		if (!mesh->getRefinementDepth(edge)){
			if ( pGCData->edgeBelongToDomain(edge,dom) ){
				if (pGCData->belongsToBoundary(edge)){
					numBE1++;
					Dij[0] = .0; Dij[1] = .0;
					pGCData->getDij(edge,dom,Dij);

					// get nodes I and J
					pEntity I = (pVertex)edge->get(0,0);
					pEntity J = (pVertex)edge->get(0,1);
					//cout << "BE: " << numBE1 << "\t" << EN_id(I)-1 << " - " << EN_id(J)-1 << "  Dij: " << Dij[0] << " , " << Dij[1] << endl;

					sumx += Dij[0];
					sumy += Dij[1];

					// get nodes' volumes
					double volumeI = pGCData->getVolume(I,dom);
					double volumeJ = pGCData->getVolume(J,dom);

#ifdef _SEEKFORBUGS_
					if ( fabs(volumeI)<1e-12 || fabs(volumeJ)<1e-12 ){
						char msg[256];
						sprintf(msg,"Volume cannot be null. dom = %d\tvolumeI[%d]: %E, volumeJ[%d] = %E\n",dom,EN_id(I),volumeI,EN_id(J),volumeJ);
						std::cout << "edgeBelongToDomain: " << pGCData->edgeBelongToDomain(edge,dom) << endl;
						throw Exception(__LINE__,__FILE__,msg);
					}
#endif //_SEEKFORBUGS_

					// fill edge matrix
					Fij_column1 [0] = 5.*Dij[0]/(6.*volumeI);
					Fij_column1 [1] = 5.*Dij[1]/(6.*volumeI);
					Fij_column1 [2] = Dij[0]/(6.*volumeJ);
					Fij_column1 [3] = Dij[1]/(6.*volumeJ);
					Fij_column2 [0] = Dij[0]/(6.*volumeI);
					Fij_column2 [1] = Dij[1]/(6.*volumeI);
					Fij_column2 [2] = 5.*Dij[0]/(6.*volumeJ);
					Fij_column2 [3] = 5.*Dij[1]/(6.*volumeJ);

					// where edge matrix must be assembled.
					int id0 = pMData->get_AppToPETSc_Ordering(EN_id(I));
					int id1 = pMData->get_AppToPETSc_Ordering(EN_id(J));
					int pos1 = 2*(id0-1);
					int pos2 = 2*(id1-1);
					int idxm[4] = {pos1,pos1+1,pos2,pos2+1};
					int idxn[2] = {id0-1,id1-1};

					ierr = MatSetValues(F,4,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);
					ierr = MatSetValues(F,4,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);
				}
			}
		}
	}
	EIter_delete(eit);

	cout << setprecision(18) << "Sum: " << sumx << "\t" << sumy << endl;

	cout << "Number of BDRY edges counted: " << numBE1 << endl;
	cout << "Number of BDRY edges counted via FMDB: " << numBE1 << endl;
	return 0;
}

int EBFV1_elliptic::F_bdryFaces(pMesh mesh, Mat F, const int &dom){
	int i;
	//int count = 0;
	//dblarray Dij(3,.0);
	double Dij[3] = {.0,.0,.0};
	std::set<int>::iterator iter;
	FIter fit = M_faceIter(mesh);
	while (pEntity face = FIter_next(fit)){
		cout << "OUT\n";
		if ( pGCData->getDij(face,dom,Dij) ){
			cout << "IN\n";

			// get face's vertices I, J and K
			pEntity I = (pVertex)face->get(0,0);
			pEntity J = (pVertex)face->get(0,1);
			pEntity K = (pVertex)face->get(0,2);

			// get node's volume
			double volumeI = pGCData->getVolume(I,dom);
			double volumeJ = pGCData->getVolume(J,dom);
			double volumeK = pGCData->getVolume(K,dom);
			cout << "Volumes: " << volumeI << "\t" << volumeJ << "\t" << volumeK << endl;

#ifdef _SEEKFORBUGS_
			if ( volumeI<1e-12 || volumeJ<1e-12 || volumeK<1e-12 )
				throw Exception(__LINE__,__FILE__,"Volume cannot be null.\n");
#endif //_SEEKFORBUGS_

			// Fij assembly
			double tmp[3] = {1./(8.*volumeI), 1./(8.*volumeJ), 1./(8.*volumeK)};
			double aux[3][3] = {{6.*tmp[0],tmp[0],tmp[0]},{tmp[1],6.*tmp[1],tmp[1]},{tmp[2],tmp[2],6.*tmp[2]}};
			double Fij_column1[9], Fij_column2[9], Fij_column3[9];
			for (i=0; i<3; i++){
				Fij_column1[3*i] = aux[i][0]*Dij[0];
				Fij_column1[3*i+1] = aux[i][0]*Dij[1];
				Fij_column1[3*i+2] = aux[i][0]*Dij[2];

				Fij_column2[3*i] = aux[i][1]*Dij[0];
				Fij_column2[3*i+1] = aux[i][1]*Dij[1];
				Fij_column2[3*i+2] = aux[i][1]*Dij[2];

				Fij_column3[3*i] = aux[i][2]*Dij[0];
				Fij_column3[3*i+1] = aux[i][2]*Dij[1];
				Fij_column3[3*i+2] = aux[i][2]*Dij[2];
			}

			// index for global Fg
			// where edge matrix must be assembled.
			int id0 = pMData->get_AppToPETSc_Ordering(EN_id(I));
			int id1 = pMData->get_AppToPETSc_Ordering(EN_id(J));
			int id2 = pMData->get_AppToPETSc_Ordering(EN_id(K));

			int pos1 = 3*(id0-1);
			int pos2 = 3*(id1-1);
			int pos3 = 3*(id2-1);

			int idxm[9] = {pos1,pos1+1,pos1+2, pos2,pos2+1,pos2+2, pos3,pos3+1,pos3+2};
			int idxn[3] = {id0-1, id1-1, id2-1};

			ierr = MatSetValues(F,9,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);CHKERRQ(ierr);
			ierr = MatSetValues(F,9,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);CHKERRQ(ierr);
			ierr = MatSetValues(F,9,idxm,1,&idxn[2],Fij_column3,ADD_VALUES);CHKERRQ(ierr);
		}
	}
	FIter_delete(fit); // end of loop over edges
	return 0;
}
}