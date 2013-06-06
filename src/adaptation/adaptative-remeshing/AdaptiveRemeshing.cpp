/*
 * AdaptiveRemeshing.cpp
 *
 *  Created on: 26/04/2013
 *      Author: Saulo Menezes
 */
#include "AdaptiveRemeshing.h"
#include "ErrorEstimator.h"
#include "MRE.h"
#include "FileManager.h"
#include "Rebuilder3.h"

void AdaptiveRemeshing::rodar(ErrorAnalysis  *pErrorAnalysis, pMesh & theMesh){
	if (!theMesh){
		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
	}

	//AOMD::AOMD_Util::Instance()->exportGmshFile("SAIDAMSHENTRADA.msh",theMesh);

	std::list<pEntity> elementList;
	std::set<pEntity> elementSet;
	pErrorAnalysis->getRefUnrefElementsList(theMesh,elementList,elementSet);
	//	STOP();
	cout << "Number of elements to be (un)refined: " << elementList.size() << endl;

	if (!elementList.size()){
		throw Exception(__LINE__,__FILE__,"No elements to be (un)refined");
	}
	else{
		theMesh->modifyState(3,2,1);
		theMesh->modifyState(3,1,1);
		theMesh->modifyState(3,0);

		theMesh->modifyState(2,1);
		theMesh->modifyState(2,3);
		theMesh->modifyState(2,0);

		theMesh->modifyState(1,3);
		theMesh->modifyState(1,2);
		theMesh->modifyState(1,0);

		theMesh->modifyState(0,2);
		theMesh->modifyState(0,1);
		theMesh->modifyState(0,3);

		ErrorEstimator *E = new ErrorEstimator(theMesh);
		E->LeitorDeLista(pErrorAnalysis, theMesh);
		E->VertexList (theMesh);

		AOMD::AOMD_Util::Instance()->exportGmshFile("Malha_0.msh",theMesh);

		FileManager *Reader = new FileManager();
		Reader->FileReader("ListadeElementos"); // Esta funcao está lendo o arquivo que contem a lista de ids dos elementos a ser removidos e salvando em um set ListofElements (lista de ids), este set será usado pelas funcoes BoundaryElements3D e BoundaryElements2D

		//int I = Reader->ReturnIdent(); // retorna o identificador, se os elementos a ser removidos sao faces ou tetras.
		//exportSolutionToVTK(theMesh,0,0,0,"vtks/theMesh_0.vtk");

		int I=0;
		if (M_numRegions(theMesh)!=0){ // verificando se a malha é 3D ou 2D...
			I=3; // malha 3D
		}
		else{
			I=2; // malha 2D
		}

		cout << "Num tetras: " << M_numRegions(theMesh) << endl;
		cout << "Num vertices: " << M_numVertices(theMesh) << endl;
		cout << "Num faces: " << M_numFaces(theMesh) << endl;

		theMesh->modifyState(3,2,1);
		theMesh->modifyState(3,1,1);
		theMesh->modifyState(3,0);

		theMesh->modifyState(2,1);
		theMesh->modifyState(2,3);
		theMesh->modifyState(2,0);

		theMesh->modifyState(1,3);
		theMesh->modifyState(1,2);
		theMesh->modifyState(1,0);

		theMesh->modifyState(0,2);
		theMesh->modifyState(0,1);
		theMesh->modifyState(0,3);

		int VerCL = 0;

		double Cl = 0.5;

		cout << "o Cl inserido é: " << Cl << endl;

		Remover *p = new Remover(theMesh);
		p->ClReader("ListaCls"); // Esta funcao lê o arquivo com a lista de pontos com seus respectivos cls
		int GrtID = p->Iterador(Reader->ReturnListofElements(), I, theMesh); // PAULO DIZ QUE O GMSH JA DA OS ELEMENTOS DE FRONTEIRA, VERIFICAR ISSO.


		// Atribuindo CL = 0 a todos os nodes
		VIter vit = M_vertexIter(theMesh);
		while (pVertex Vertex = VIter_next(vit)){
			EN_attachDataDbl (Vertex, MD_lookupMeshDataId("CL"), 0);
		}
		VIter_delete(vit);

		if (I==2){
			/*
			 * recebe o set de faces a remover set<pFace> BoundaryFaces e a partir dele preenche os boundaryedges e boundarynodes,
			 * inicialmente recebe todos os faces (nao apenas os de contorno), mas depois das funcoes de remoção só sobram os boundary mesmo, daí o nome do set.
			 */
			p->BoundaryElements2D(Reader->ReturnListofElements(), Cl, VerCL,theMesh);

			// Pegando o valor da maior aresta ANTES de rodar a funcao Boundaryelements2D novamente, quando ela rodar de novo, esse valor vai ser substituido, por isso tenho que guardá-lo agora.
			float GE = p->GreatestEdge;
			p->SaveBGMView1(Cl,GE); // Salvando o arquivo da View da malha de background.
			VerCL=1;

			//   AGORA TENHO QUE CHAMAR AS FUNCOES IDENTIFICADORAS.
			int PSCounter = 0;  // Contador de plane surfaces

			//PSCounter=p->identifysurfaces_2D(theMesh); // Esta funcao identifica os Plane Surfaces antes da remocao dos elementos, assim se quando a funcao RemoveStgElements remover mais algum, será necessario identificar os novos elementos de contorno também.

			p->removetetras(theMesh); // removendo os tetraedros indicados na lista
			p->removeinternalfaces(theMesh); // removendo as faces indicadas na lista

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			p->removestgNodes(theMesh);  // Não está removendo os elementos da malha... Só está removendo os do contorno...

			theMesh->modifyState(1,2,0);
			theMesh->modifyState(1,2,1);
			theMesh->modifyState(2,1,0);
			theMesh->modifyState(2,1,1);

			PSCounter=p->Identify_and_Remove_Edges_2D(theMesh);
			p->removeinternaledges(theMesh); // removendo as arestas das faces que foram removidas

			theMesh->modifyState(0,1,0);
			theMesh->modifyState(0,1);

			theMesh->modifyState(0,1,0);
			theMesh->modifyState(0,1,1);
			theMesh->modifyState(0,2,0);
			theMesh->modifyState(0,2,1);


			// Vou deixar as funcoes que removem elementos da fronteira desligadas por enquanto.
			// Esta funcao remove elementos "estranhos" (dentes nos contornos, elementos que estejam causando uniao de lineloops, faces 		soltas no interior do buraco...)
			// cout << "Vai comecar SaveBGMView2" << endl;
			p->SaveBGMView2(GE); // Salvando o arquivo da View da malha de background.
			p->removeinternalnodes(theMesh);

			// Agora salvando os arquivos:
			FileManager *fm;
			fm = new FileManager(theMesh, p->ReturnBoundaryEdges(), p->ReturnBoundaryNodes(),PSCounter);
			fm->SaveFileGeometry_2D(GE, theMesh,PSCounter);  // Neste momento o cl dos elementos é o comprimento da maior aresta
			fm->ResetFileManager();

			p->removeexternaledges(theMesh);
			p->removeexternalnodes(theMesh);
			cout << "Procurando o SF" << endl;

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			// Agora a chamando as funcoes MakeMesh e RefineMesh para recriar a malha com a geometria recém salva em fm->SaveFileGeometry_2D();
			Rebuilder3 *Q;
			Q = new Rebuilder3();
			Q->MakeMesh_2D("Geometria2D");

			// A partir daqui, começa a geração da malha de adaptação

			pMesh theMesh2;
			theMesh2 = MS_newMesh(0);
			cout << "Lendo Geometria2D.msh" << endl;
			M_load(theMesh2,"Geometria2D.msh");
			cout << "Lido o arquivo de malha gerado" << endl;
			cout << "Num edges: " << M_numEdges(theMesh2) << endl;


			theMesh2->modifyState(3,2,1);
			theMesh2->modifyState(3,1,1);
			theMesh2->modifyState(3,0);

			theMesh2->modifyState(2,1);
			theMesh2->modifyState(2,3);
			theMesh2->modifyState(2,0);

			theMesh2->modifyState(1,3);
			theMesh2->modifyState(1,2);
			theMesh2->modifyState(1,0);

			theMesh2->modifyState(0,2);
			theMesh2->modifyState(0,1);
			theMesh2->modifyState(0,3);

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);
			double GID = (double)GrtID;
			fm->SaveFileMshRenumber(theMesh, "theMesh.msh",0); // Este terceiro parametro é o ponto inicial de contagem dos ids
			fm->SaveFileMshRenumber2(theMesh2, "theMesh2.msh",GID);

			Q->MakeMerge_2D("theMesh.msh", "theMesh2.msh","Final_01.msh");

			cout<<" ERRO "<<endl;
			theMesh = MS_newMesh(0);
			deleteMesh(theMesh);
			M_load(theMesh,"Final_01.msh");

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			//  AOMD::AOMD_Util::Instance()->exportGmshFile("SAIDAFINALLOOP.msh",theMesh);
			cout << "A MAIOR ARESTA É " << GE << endl;
			delete E; E = 0;
			delete p; p = 0;
			delete Q; Q = 0;
			delete Reader; Reader = 0;
			delete fm; fm = 0;

		}
	}

}
