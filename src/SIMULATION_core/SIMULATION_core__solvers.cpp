/*
 * SIMULATION_core__solvers.cpp
 *
 *  Created on: 27/08/2012
 *      Author: rogsoares
 */

#include "SIMULATION_core.h"
#include <sstream>

void readmesh(pMesh m,char* filename);
void setMeshFile();

namespace PRS{

	int SIMULATION_core::solver(){
		switch ( simFlag ){
		case STEADY_STATE:
			this->steadyState();
			break;
		case TRANSIENT:
			this->transient();
			break;
		default:
			throw Exception(__LINE__,__FILE__, Exception::INIT_ERROR );
		}
		return 0;
	}

	/*                       ATENCAO
	 * A chamada do solver de presso faz usa da malha apontada por theMesh, logo esta  quem deve ser adaptada. Isto  feito por meio
	 * da malha pIData->m1 e NUNCA por pIData->m2. Aps a adaptacao, a interpolacao  de pIData->m2 para pIData->m1, o preprocessamento
	 *  sobre pIData->m1. Quando o ciclo de adaptacao se encerrar, o solver, que faz usa de theMesh, ira calcular o novo campo de pressoes
	 * sobre pIData->m1, pois este ponteiro aponta para o mesmo endereco de memoria que theMesh.
	 *
	 * Rogrio, 15/05/2013
	 */
	int SIMULATION_core::steadyState(){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
		bool adapt;
		double tol1 = pSimPar->getToleranceForAllElements();
		double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();
		int numFields = 1;
	
		pPPData->setSimulationState(false);
		pIData->m1 = theMesh;
		PADMEC_mesh *pm = new PADMEC_mesh;
		int count = 0;

		do{
			pElliptic_eq->solver(theMesh);
			adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),numFields);
			pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
			if ( adapt ){
				makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSaturation_Old);
				pMeshAdapt->rodar(pErrorAnalysis,pIData->m1);
				deleteMesh(pIData->m1);pIData->m1 = 0;
				deleteMesh(theMesh); theMesh = 0;
				pIData->m1 = MS_newMesh(0);
				readmesh(pIData->m1,"Final_01.msh");
				theMesh = pIData->m1;
				makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);
				PADMEC_GAMBIARRA(pIData->m1);
				cout<< "Interpolador"<<endl;
				interpolation(pIData,pSimPar->getInterpolationMethod());
				pSimPar->printOutVTK(theMesh,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
				//STOP();
				EBFV1_preprocessor(pIData->m1,pGCData);
	#ifdef __ADAPTATION_DEBUG__
				validate_EBFV1(pGCData,pIData->m1,pSimPar->setOfDomains);
				if (!pSimPar->setOfDomains.size()){
					throw Exception(__LINE__,__FILE__,"Num domains NULL!\n");
				}
	#endif
				updatePointersData(theMesh);
				// After interpolation of pressure field calculate (not interpolate) the gradient pressure for the new mesh
				deleteMesh(pm);
				deleteMesh(pIData->m2); pIData->m2 = 0;
				pIData->m2 = MS_newMesh(0);
			}
			count++;
		}while (adapt);

		PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
		cout<< "Loops :"<< count <<endl;
		return 0;
	}
	
	int SIMULATION_core::transient(){
		PetscPrintf(PETSC_COMM_WORLD,"\n\nStart simulation:\n-----------------------------------------------\n");
		LogFiles(OPENLG,0,0,0,0,pSimPar->getOutputPathName(),pSimPar->useRestart(),pSimPar->getTStepNumber(),pSimPar->getCPU_time());

		double timeStep;
		double time_step_summation = .0;
		double adaptStep = 0;				// performs mesh adaptation every 20 timeSteps
		bool adapt;
		double tol1 = pSimPar->getToleranceForAllElements();
		double tol2 = pSimPar->getToleranceForAllElements_excludingSingularities();
		int numFields = 2;
		pIData->m1 = theMesh;
		PADMEC_mesh *pm = new PADMEC_mesh;

		while ( !pSimPar->finishSimulation() ){
			pElliptic_eq->solver(theMesh);
			//pSimPar->printOutVTK(pIData->m1,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
			pHyperbolic_eq->solver(theMesh,timeStep);
			pSimPar->printOutVTK(pIData->m1,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
			if ( pSimPar->userRequiresAdaptation() ){
				adapt = calculate_ErrorAnalysis(pErrorAnalysis,theMesh,pSimPar,tol1,tol2,pPPData->get_getPFuncArray(),numFields);
				//pSimPar->printOutVTK(pIData->m1,pPPData,pErrorAnalysis,pSimPar,exportSolutionToVTK);
				if (adapt){
					makeMeshCopy2(pIData->m1,pm,pPPData->getPressure,pPPData->getSaturation_Old);
					pMeshAdapt->rodar(pErrorAnalysis,pIData->m1);
					deleteMesh(pIData->m1);pIData->m1 = 0;
					deleteMesh(theMesh); theMesh = 0;
					pIData->m1 = MS_newMesh(0);
					readmesh(pIData->m1,"Final_01.msh");
					theMesh = pIData->m1;
					makeMeshCopy2(pm,pIData->m2,pPPData->setPressure,pPPData->setSaturation);
					PADMEC_GAMBIARRA(pIData->m1);
					cout<< "Interpolador"<<endl;
					interpolation(pIData,pSimPar->getInterpolationMethod());
					//pInterpolateData(pIData);
					EBFV1_preprocessor(pIData->m1,pGCData);

					#ifdef __ADAPTATION_DEBUG__
						validate_EBFV1(pGCData,pIData->m1,pSimPar->setOfDomains);
						if (!pSimPar->setOfDomains.size()){
							throw Exception(__LINE__,__FILE__,"Num domains NULL!\n");
						}
					#endif

					updatePointersData(theMesh);
					deleteMesh(pm);
					deleteMesh(pIData->m2); pIData->m2 = 0;
					pIData->m2 = MS_newMesh(0);
				}
			}
		}

		PetscPrintf(PETSC_COMM_WORLD,"\n\nEnd of simulation:\n-----------------------------------------------\n");
		return 0;
	}
}

void setMeshFile(){

	pMesh m;
	m = MS_newMesh(0);

	readmesh(m,"Final_01.msh");
	PADMEC_GAMBIARRA(m);

	double coords[3];
	ofstream fid;
	fid.open("AdaptedMesh.msh");
	fid << "$NOD\n";
	fid << M_numVertices(m) << endl;
	VIter vit = M_vertexIter(m);
	pEntity node;
	while( (node = VIter_next(vit)) ){
		V_coord(node,coords);
		fid << EN_id(node) << " " << coords[0] << " " << coords[1]<< " " << coords[2] << endl;
	}
	VIter_delete(vit);
	fid << "$ENDNOD\n$ELM\n";

	// count elements flagged
	std::list<pEntity> flaggedVertices;
	vit = M_vertexIter(m);
	while( (node = VIter_next(vit)) ){
		if ( node->getClassification() ){
			int flag = GEN_tag(node->getClassification());
			cout << flag << endl;
			if ( flag==51 || flag==10 || flag==1100 ){
				flaggedVertices.push_back(node);
			}
		}
	}
	VIter_delete(vit);
	std::list<pEntity> flaggedEdges;
	pEntity edge;
	EIter eit = M_edgeIter(m);
	while( (edge = EIter_next(eit)) ){
		if ( edge->getClassification() ){
			int flag = GEN_tag(edge->getClassification());
			if ( flag==2000 ){
				flaggedEdges.push_back(edge);
			}
		}
	}
	EIter_delete(eit);

	fid << (int)flaggedVertices.size() + (int)flaggedEdges.size() + M_numFaces(m) << endl;
	int k = 0;
	std::list<pEntity>::iterator iter;
	for(iter=flaggedVertices.begin(); iter != flaggedVertices.end(); iter++){
		fid << ++k << " 15 " << GEN_tag((*iter)->getClassification()) << " 1 1 " << EN_id(*iter) << endl;
	}
	for(iter=flaggedEdges.begin(); iter != flaggedEdges.end(); iter++){
		fid << ++k << " 1 2000 1 2 " << EN_id((*iter)->get(0,0)) << " " << EN_id((*iter)->get(0,1)) << endl;
	}
	pEntity face;
	FIter fit = M_faceIter(m);
	while( (face = FIter_next(fit)) ){
		fid << ++k << " 2 3300 1 3 " << EN_id(face->get(0,0)) << " " << EN_id(face->get(0,1)) << " " << EN_id(face->get(0,2)) << endl;
	}
	FIter_delete(fit);
	fid << "$ENDELM\n";
	fid.close();
}

void readmesh(pMesh m,char* filename){
	cout << "Lendo malha... ";
	ifstream fid;
	fid.open(filename);
	int NbNod;
	int iNod;
	double x,y,z;
	char line[256];
	fid.getline (line,256);
	fid >> NbNod;
	cout << "\nNodes: " << NbNod << endl;
	for(int i=0;i<NbNod;i++){
		fid >> iNod >> x >> y >> z;
		m->createVertex(iNod,x,y,z,0);
	}
	fid.getline (line,256);
	fid.getline (line,256);
	fid.getline (line,256);
	fid >> NbNod;
	cout << "\nElements: " << NbNod << endl;
	int face_id = 0;
	for (int i=0; i<NbNod; i++){
		int iNbNod,iTyp,iGrp,iElm,iNbSub,id;
		mVertex *nod[100];
		fid >> iElm >> iTyp >> iGrp >> iNbSub >> iNbNod;
		for(int i=0;i<iNbNod;i++){
			fid >> id;
			nod[i] = m->getVertex(id);
		}

		mEntity *theEntity = 0;
		switch(iTyp){
		case 2 :
			theEntity = m->createFaceWithVertices(nod[0],nod[1],nod[2],m->getGEntity(iGrp,2));
			EN_setID((pEntity)theEntity,++face_id);
			break;
		case 1 :
			theEntity = m->createEdge(nod[0],nod[1],m->getGEntity(iGrp,1));
			break;
		case 15 :
			mVertex *v = nod[0];{
				v->classify(m->getGEntity(iGrp,0));
			}
		}
	}
	cout << "OK!\n";
}

