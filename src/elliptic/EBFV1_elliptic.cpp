#include "EBFV1_elliptic.h"
#include <time.h>

namespace PRS           // PRS: Petroleum Reservoir Simulator
{
	EBFV1_elliptic::EBFV1_elliptic(){
	}
	
	EBFV1_elliptic::EBFV1_elliptic(pMesh mesh, PhysicPropData *ppd, SimulatorParameters *sp, GeomData *gcd, MeshData *md){
		matvec_struct = new Data_struct;
		pPPData = ppd;
		pGCData = gcd;
		pSimPar = sp;
		pMData = md;
		theMesh = mesh;
		one_eighth = (double)(1.0/8.0);
	}
	
	EBFV1_elliptic::~EBFV1_elliptic(){
		if (matvec_struct){
			delete[] matvec_struct->rows;
			matvec_struct->rows = 0;
			delete matvec_struct;
			matvec_struct = 0;
		}
	}
	
	// solves system of equation for pressure field
	double EBFV1_elliptic::solver(pMesh theMesh){
		if (pSimPar->userRequiresAdaptation()){
			matvec_struct = new Data_struct;
		}
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tIN\n";
		#endif
		assembly_EFG_RHS(theMesh);
		setMatrixFreeOperation(theMesh);
		updatePressure(theMesh);
		calculatePressureGradient();
		freeMemory();
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: pressure solver\tOUT\n";
		#endif
		return 0;
	}
	
	double EBFV1_elliptic::updatePressure(pMesh theMesh){
		double startt = MPI_Wtime();
		
		PetscScalar *sol, val;
		PetscInt i,m,n,row,col=0;
		PetscInt numGN = pMData->getNum_GNodes();
		Mat mSol,mLSol;
		
		// create a column matrix to receive output vector values (mSol)
		MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,numGN,1,0,PETSC_NULL,0,PETSC_NULL,&mSol);
		
		int nLIDs, *IDs_ptr;
		pMData->getRemoteIDs(nLIDs,&IDs_ptr);
		
		// transference process: from vector to column matrix
		VecGetArray(output,&sol);
		MatSetValues(mSol,matvec_struct->nrows,matvec_struct->rows,1,&col,sol,INSERT_VALUES);
		VecRestoreArray(output,&sol);
		MatAssemblyBegin(mSol,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(mSol,MAT_FINAL_ASSEMBLY);
		//VecView(output,PETSC_VIEWER_STDOUT_WORLD); //throw 1;
		
		// remote values cannot be gotten from remote matrix positions. transfer (via MatGetSubMatrixRaw) necessary remote values 
		// to each process to a second column matrix (mLSol).
		MatGetSubMatrixRaw(mSol,nLIDs,IDs_ptr,1,&col,PETSC_DECIDE,MAT_INITIAL_MATRIX,&mLSol);
		MatDestroy(mSol);
		MatGetOwnershipRange(mLSol,&m,&n);
		
		// loop over IDs_ptr[i]
		row = m;
		for(i=0; i<nLIDs;i++){
			int ID = pMData->get_PETScToApp_Ordering(IDs_ptr[i]+1);
			MatGetValues(mLSol,1,&row,1,&col,&val);
			pPPData->setPressure(ID-1,val);
			row++;
		}
		MatDestroy(mLSol);
		
		static bool key = true;
		if (key){
			int nnodes = M_numVertices(theMesh);
			for(i=0;i<nnodes;i++){
				int ID = pMData->get_AppToPETSc_Ordering(i+1);
				if ( pMData->getDirichletValue(ID,&val) ){
					pPPData->setPressure(i,val);
				}
			}
			key = false;
		}
		return 0;
	}
	
	double EBFV1_elliptic::freeMemory(){
		double startt = MPI_Wtime();
		// free matrices
		if (!pSimPar->useDefectCorrection()){
			MatDestroy(matrix);
		}
		MatDestroy(matvec_struct->G);
		for(int i=0; i<pSimPar->getNumDomains(); i++){
			MatDestroy(matvec_struct->E[i]);
			MatDestroy(matvec_struct->F[i]);
		}
		
		/// free vectors
		VecDestroy(matvec_struct->RHS);
		VecDestroy(matvec_struct->z);
		VecDestroy(output);
		delete[] matvec_struct->F;
		delete[] matvec_struct->E;
		matvec_struct->F = 0;
		matvec_struct->E = 0;
		if (pSimPar->userRequiresAdaptation()){
			delete[] matvec_struct->rows; matvec_struct->rows = 0;
			delete matvec_struct; matvec_struct = 0;
		}
		return MPI_Wtime()-startt;
	}
}
