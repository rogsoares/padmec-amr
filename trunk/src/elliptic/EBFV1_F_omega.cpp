#include "EBFV1_elliptic.h"

namespace PRS
{
	int EBFV1_elliptic::gradient_F_edges(Mat F, pEntity edge, const int &dom, int dim, dblarray &Cij){
		int i;

		// get nodes I and J
		pEntity I = (pVertex)edge->get(0,0);
		pEntity J = (pVertex)edge->get(0,1);
		if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
			std::swap(I,J);
		}

		// Where Fij should be assembled
		int id0 = pMData->get_AppToPETSc_Ordering(EN_id(I));
		int id1 = pMData->get_AppToPETSc_Ordering(EN_id(J));

		//double nrc = pGCData->getNumRemoteCopies(edge) + 1.0;
		//double nrc = (double)pGCData->getNumRC(edge,dom)+1.0;
		double nrc = 1.0;//(double)pGCData->getNumRC(theMesh,edge) + 1.0;

		double volumeI = pGCData->getVolume(I,dom);
		double volumeJ = pGCData->getVolume(J,dom);

		//printf("Vol[%d} - [%d]: %.8f\t[%d]: %.8f\n",dom,EN_id(I),volumeI,EN_id(J),volumeJ);

	#ifdef _SEEKFORBUGS_
		if ( fabs(volumeI)<1e-12 || fabs(volumeJ)<1e-12 ){
			char msg[256];
			sprintf(msg,"Volume cannot be null. dom = %d\tvolumeI[%d]: %E, volumeJ[%d] = %E\n",dom,EN_id(I),volumeI,EN_id(J),volumeJ);
			std::cout << "edgeBelongToDomain: " << pGCData->edgeBelongToDomain(edge,dom) << endl;
			throw Exception(__LINE__,__FILE__,msg);
		}
	#endif //_SEEKFORBUGS_

		double aux1 = (1./(2.*volumeI))/nrc;
		double aux2 = (-1./(2.*volumeJ))/nrc;
		// Assembly Fij
		double Fij_column1[2*dim], Fij_column2[2*dim];
		for (i=0; i<dim; i++){
			Fij_column1[i] = Cij[i]*aux1;
			Fij_column2[i] = Fij_column1[i];
			Fij_column1[i+dim] = Cij[i]*aux2;
			Fij_column2[i+dim] = Fij_column1[i+dim];
		}

		int pos1 = dim*(id0-1);
		int pos2 = dim*(id1-1);
		// says on which rows Fij must be assembled into Fg
		int idxm[2*dim];
		for (i=0; i<dim; i++){
			idxm[i] = pos1+i;
			idxm[dim+i] = pos2+i;
		}
		// says on which columns Fij must be assembled into Fg
		int idxn[2] = {id0-1,id1-1};

//		if (id0==7 && id1==11) {
//
//			static ofstream fid;
//			static bool key=true;
//			if (key){
//				char fname[256]; sprintf(fname,"Fij__%d-of-%d.txt",P_pid(),P_size());
//				fid.open(fname);
//				key=false;
//			}
//			fid << "\nedge: " << id0 << "  " << id1 << " nrc = " << nrc-1.0 << " dom  = " << dom << endl;
//			cout << "edge: " << id0 << "  " << id1 << " nrc = " << nrc-1.0 << endl;
//			for (i=0;i<2*dim;i++) fid << Fij_column1[i] << " " << Fij_column2[i] << endl;
		ierr = MatSetValues(F,2*dim,idxm,1,&idxn[0],Fij_column1,ADD_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(F,2*dim,idxm,1,&idxn[1],Fij_column2,ADD_VALUES);CHKERRQ(ierr);
		return 0;
//		}
	}
}
