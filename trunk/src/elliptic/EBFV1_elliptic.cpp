#include "EBFV1_elliptic.h"
#include <time.h>

namespace PRS           // PRS: Petroleum Reservoir Simulator
{
	EBFV1_elliptic::EBFV1_elliptic(){
	}

	EBFV1_elliptic::EBFV1_elliptic(pMesh mesh, PhysicPropData *ppd,
			SimulatorParameters *sp, GeomData *gcd,
			MeshData *md){
		matvec_struct = new Data_struct;
		pPPData = ppd;
		pGCData = gcd;
		pSimPar = sp;
		pMData = md;
		theMesh = mesh;
		Lij.assign(3,.0);
		pVec = new Vectors;
		DF_key = true;
	}

	EBFV1_elliptic::~EBFV1_elliptic(){
		delete[] matvec_struct->rows;
		delete matvec_struct;
	}

	// solves system of equation for pressure field
	double EBFV1_elliptic::solver(){
		if (!P_pid()) std::cout << "Elliptic solver...";
		double cpu_time = assembly_EFG_RHS();
//		cout << "assembly time: " << cpu_time << endl;
		// which scheme should be used to solve pressure field
		if (pSimPar->useDefectCorrection())
			cpu_time += solveIteratively();
		else
			cpu_time += setMatrixFreeOperation();
//		cout << "solver time: " << cpu_time << endl;
		cpu_time += updatePressure();
//		cout << "updatePressure time: " << cpu_time << endl;
#ifdef CRUMPTON_EXAMPLE
		// Output data (VTK)
		//pSimPar->printOutVTK(theMesh,pPPData,exportSolutionToVTK);
		MPI_Barrier(MPI_COMM_WORLD);
		// terminate elliptic equation evaluation
		char msg[256]; sprintf(msg,"CRUMPTON_EXAMPLE:\t\tCPU time elapsed: %f\n",cpu_time);
		throw Exception(__LINE__,__FILE__,msg);
#endif
		cpu_time += pressureGradient();
//		cout << "pressureGradient time: " << cpu_time << endl;
	//	STOP();
		cpu_time += freeMemory();
		if (!P_pid()) std::cout << "done.\n\n";
		return cpu_time;
	}

	double EBFV1_elliptic::updatePressure(){
		double startt = MPI_Wtime();

		PetscScalar *sol, val;
		PetscInt i,m,n,row,col=0;
		PetscInt numGN = pMData->getNum_GNodes();
		Mat mSol,mLSol;

		// create a column matrix to receive output vector values (mSol)
		ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,numGN,1,0,PETSC_NULL,0,PETSC_NULL,&mSol);CHKERRQ(ierr);

		int nLIDs, *IDs_ptr;
		pMData->getRemoteIDs(nLIDs,&IDs_ptr);

		// transference process: from vector to column matrix
		ierr = VecGetArray(output,&sol);CHKERRQ(ierr);
		ierr = MatSetValues(mSol,matvec_struct->nrows,matvec_struct->rows,1,&col,sol,INSERT_VALUES);CHKERRQ(ierr);
		ierr = VecRestoreArray(output,&sol);CHKERRQ(ierr);
		ierr = MatAssemblyBegin(mSol,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(mSol,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		//ierr = VecView(output,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); //throw 1;

		// remote values cannot be gotten from remote matrix positions.
		// transfer (via MatGetSubMatrixRaw) necessary remote values to each process
		// to a second column matrix (mLSol).
		ierr = MatGetSubMatrixRaw(mSol,nLIDs,IDs_ptr,1,&col,PETSC_DECIDE,MAT_INITIAL_MATRIX,&mLSol);CHKERRQ(ierr);
		ierr = MatDestroy(mSol);CHKERRQ(ierr);
		ierr = MatGetOwnershipRange(mLSol,&m,&n);CHKERRQ(ierr);

		// loop over IDs_ptr[i]
		row = m;
		for(i=0; i<nLIDs;i++){
			int ID = pMData->get_PETScToApp_Ordering(IDs_ptr[i]+1);
			pVertex node = theMesh->getVertex(ID);
			if (!node) throw Exception(__LINE__,__FILE__,"Node does not exist.\n");
			ierr = MatGetValues(mLSol,1,&row,1,&col,&val);CHKERRQ(ierr);
			pPPData->setPressure(node,val);
			//printf("p[%d] = %f\n",EN_id(node),val);
			row++;
		}
		//STOP();
		ierr = MatDestroy(mLSol);CHKERRQ(ierr);

		static bool key = true;
		if (key){
			VIter vit = M_vertexIter(theMesh);
			while (pEntity node = VIter_next(vit)){
				int ID = pMData->get_AppToPETSc_Ordering(EN_id(node));
				if ( pMData->getDirichletValue(ID,&val) )
					pPPData->setPressure(node,val);
			}
			VIter_delete(vit);
			key = false;
		}
		return MPI_Wtime()-startt;
	}

	double EBFV1_elliptic::freeMemory(){
		double startt = MPI_Wtime();
		// free matrices
		if (!pSimPar->useDefectCorrection()){
			ierr = MatDestroy(matrix);CHKERRQ(ierr);
		}
		ierr = MatDestroy(matvec_struct->G);CHKERRQ(ierr);
		for(int i=0; i<pSimPar->getNumDomains(); i++){
			ierr = MatDestroy(matvec_struct->E[i]);CHKERRQ(ierr);
			ierr = MatDestroy(matvec_struct->F[i]);CHKERRQ(ierr);
		}
		/// free vectors
		ierr = VecDestroy(matvec_struct->RHS);CHKERRQ(ierr);
		ierr = VecDestroy(matvec_struct->z);CHKERRQ(ierr);
		ierr = VecDestroy(output);CHKERRQ(ierr);
		double endt = MPI_Wtime();
		return endt-startt;
	}
}
