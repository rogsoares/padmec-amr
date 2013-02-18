#include "EBFV1_elliptic.h"

// =============================================================================
// Matrix G is equivalent to the matrix of discretization of elliptic terms of
// DFVC method in a dual Voronoi mesh. G is assembled with a loop over all edges.
// In parallel, if an edge is on the partition boundary, Gij must be divided by
// (NRC) Number of Remote Copies. Each partition will add Gij/NRC to the global
// matrix E and the resulting Gij contribution will be the same as it was assem-
// bled in serial.
// =============================================================================

//static int gcount = 0;

namespace PRS
{
	int EBFV1_elliptic::divergence_G(Mat G, pEntity edge, const int &dom, int dim, dblarray &Cij){
		int i, j;
		double edgeLength, sign = 1.0;

		// get nodes I and J
		pEntity I = (pVertex)edge->get(0,0);
		pEntity J = (pVertex)edge->get(0,1);

		// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
		if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
			std::swap(I,J);
			sign = -sign;
		}

		// index for global Fg
		const int id0 = pMData->get_AppToPETSc_Ordering(EN_id(I));
		const int id1 = pMData->get_AppToPETSc_Ordering(EN_id(J));

		// permeability tensor
		const double *Permeability = pSimPar->getPermeability(dom);

		// nodal mobility
		double MobI = pPPData->getTotalMobility(I);
		double MobJ = pPPData->getTotalMobility(J);

		// average mobility
		double MobIJ = 0.5*(MobI + MobJ);

		pGCData->getEdgeVec_Unitary(edge,Lij);
		for (i=0; i<dim; i++) Lij[i] *= sign;
		edgeLength = pGCData->getEdgeLength(edge);

		// ####################################################################
		// what is calculated below: Gij = -(K*Mob_IJ)*Lij/(length)*Cij*|1 -1|
		//											                    |-1 1|
		// ####################################################################
		// KL = Permeability * Lij
		double KL[3] = {.0, .0, .0};

		int pos = 0;
		if ( pSimPar->is_K_Isotropic() ){
			for (i=0; i<dim; i++){
				KL[i] = Permeability[pos]*Lij[i];
				pos += dim+1;
			}
		}
		else{
			for (i=0; i<dim; i++)
				for (j=0; j<dim; j++) KL[i] += Permeability[dim*i+j]*Lij[j];
		}

		double aux = .0;
		for (i=0; i<dim; i++) aux += Cij[i]*KL[i];
		aux /= -edgeLength;
		double nrc = 1.0;
		//			// if an edge is on partition boundary, the result from all processors
		//			//  contribution for this edge must be equal as if was ran in serial
		//			// Gij (serial) = Gij (parallel) if edge is on partition boundary
		aux /= nrc;

		double Gij_I[2] = { MobIJ*aux,-MobIJ*aux};
		double Gij_J[2] = {-MobIJ*aux, MobIJ*aux};
		int rows[2] = {id0-1,id1-1};	// where Gij must be assembled in global G
		int cols[2] = {id0-1,id1-1};	// where Gij must be assembled in global G

		ierr = MatSetValues(G,1,&rows[0],2,cols,Gij_I,ADD_VALUES);CHKERRQ(ierr);
		ierr = MatSetValues(G,1,&rows[1],2,cols,Gij_J,ADD_VALUES);CHKERRQ(ierr);
		return 0;
	}
}
