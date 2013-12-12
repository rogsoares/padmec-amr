#include "EBFV1_elliptic.h"

// =============================================================================
// Matrix E is responsible to project the gradient over an orthogonal plane to
// edge. E is assembled with a loop over all edges. In parallel, if an edge is
// on the partition boundary, Eij must be divided by (NRC) Number of Remote
// Copies. Each partition will add Eij/NRC to the global matrix E and the resul-
// ting Eij contribution will be the same as it was assembled in serial.
// =============================================================================


namespace PRS
{
int EBFV1_elliptic::divergence_E(Mat E, pEntity edge, const int &dom, int dim, dblarray &Cij){
	const double I2D[4] = {1.0,.0,.0,1.0};
	const double I3D[9] = {1.0,.0,.0,.0,1.0,.0,.0,.0,1.0};
	int i, j, k;
	double sign = 1.0;
	// get nodes I and J
	pEntity I = (pVertex)edge->get(0,0);
	pEntity J = (pVertex)edge->get(0,1);

	// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
	if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
		std::swap(I,J);
		sign = -sign;
	}

	// index for global Fg
	int id0 = pMData->get_AppToPETSc_Ordering(EN_id(I));
	int id1 = pMData->get_AppToPETSc_Ordering(EN_id(J));
		
// 	cout << "node : " << EN_id(I) << "\t assembled in: " << id0 << endl; 
// 	cout << "node : " << EN_id(J) << "\t assembled in: " << id1 << endl;
	//STOP();

	// get absolute K tensor: dim x dim
	const double *K = pSimPar->getPermeability(dom);

	// get mobility for nodes I and J
	const double MobI = pPPData->getTotalMobility(I);
	const double MobJ = pPPData->getTotalMobility(J);
	const double MobIJ = 0.5*(MobI + MobJ);

	// ####################################################################
	// what is calculated below: Eij = -(K*Mob_IJ)/2*|I-Lij*Lij I-Lij*Lij|*Cij
	//												 |Lij*Lij-I Lij*Lij-I|
	// ####################################################################

	// product between two vectors: Lij*Lij => matrix 3x3
	k = 0;
	pGCData->getEdgeVec_Unitary(edge,Lij);
//	printf("Simadapt-- versor: %f %f\n",Lij[0],Lij[1]);
	for (i=0; i<dim; i++) Lij[i] *= sign;
	double matLij[dim*dim];
	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++) matLij[k++] = Lij[i]*Lij[j];

	const double *Identity = (dim==2)?I2D:I3D;
	double matSubtrac_ILij[dim*dim];
	for (k=0; k<dim*dim; k++) matSubtrac_ILij[k] = -0.5*(Identity[k] - matLij[k]);

	dblarray EA(dim*dim,.0);
	k = 0;
	int pos = 0;
	if ( pSimPar->is_K_Isotropic() ){
		int pos1 = 0;
		int pos2 = 0;
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++) EA[k++] = K[pos1]*matSubtrac_ILij[pos2+j];
			pos1 += dim+1;
			pos2 += dim;
		}
	}else{
		for (i=0; i<dim; i++){
			for (j=0; j<dim; j++){
				for (k=0; k<dim; k++)
					EA[pos] += K[dim*i+k]*matSubtrac_ILij[dim*k+j];
				pos++;
			}
		}
	}

	/*
	 * Se o meio é isotropico:
	 *
	 *     |K00*matSubtrac_ILij00   K00*matSubtrac_ILij01|
	 * EA =|                                             |
	 *     |K11*matSubtrac_ILij10   K11*matSubtrac_ILij11|
	 *
	 *
	 * */

	double Eij_part1[3] = {.0,.0,.0};
	for (i=0; i<dim; i++){
		for (j=0; j<dim; j++)
			Eij_part1[i] += Cij[j]*EA[dim*j+i];

		// if an edge is on partition boundary, the result from all processors
		//  contribution for this edge must be equal as if was ran in serial
		// Eij (serial) = Eij (parallel) if edge is on partition boundary
		double nrc = 1.0;// (double)pGCData->getNumRC(theMesh,edge) + 1.0;
		Eij_part1[i] *= MobIJ/nrc;
	}

	// fill Eij matrix
	double Eij_row1[2*dim], Eij_row2[2*dim];
	for (i=0; i<dim; i++){
		Eij_row1[i] = Eij_part1[i];
		Eij_row1[i+dim] = Eij_row1[i];
		Eij_row2[i] = -Eij_row1[i];
		Eij_row2[i+dim] = -Eij_row1[i];
	}

	// where Eij must be assembled in global E matrix
	int pos1 = dim*(id0-1);
	int pos2 = dim*(id1-1);
	int idxn[2*dim];
	int idxm[2] = {id0-1,id1-1};
	for (i=0; i<dim; i++){
		idxn[i] = pos1+i;
		idxn[dim+i] = pos2+i;
	}

	ierr = MatSetValues(E,1,&idxm[0],2*dim,idxn,Eij_row1,ADD_VALUES); CHKERRQ(ierr);
	ierr = MatSetValues(E,1,&idxm[1],2*dim,idxn,Eij_row2,ADD_VALUES); CHKERRQ(ierr);
	K = 0;

	return 0;
}
}

