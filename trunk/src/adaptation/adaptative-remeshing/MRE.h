#include "includes.h"
#include <string>

class Remover
{
public:

	// Contructor and Destructor
	Remover(pMesh theMesh);
	~Remover();


	set<pEdge> PlaneSurfaceEdges;
	set<pVertex> GeometryBoundNodes; // acho que este está sem uso... verificar...
	set<pVertex> GeomBoundNodes;
	set<pEdge> GeomBoundEdges;
	set<pEdge> MeshBoundEdges;
	set<pEdge> MeshBoundNodes;
	set<pVertex> AdaptNodes;
	set<pEdge> AdaptEdges;
	set<pVertex> BoundaryNodes;
	set<pEdge> BoundaryEdges;
	set<pEdge> RemovefromBoundEdges; // tive que colocar pra nao dar erro na contagem.
	set<pVertex> RemovefromBoundNodes;
	set<pFace> RemovefromBoundFaces;
	set<pFace> BoundaryFaces; // transformei o antigo set Faces em set BoundaryFaces
	set<pFace> GeomBoundFaces; // criado o GeomBoundFaces

	set<pRegion> Tetras; // Originalmente so este set é pRegion e o resto é pFace.
	vector <int> argi;
	map <string,pVertex> CommonVertex;
	float GreatestEdge;

	set <pFace> triangles;
	map <int,float> VertexClsAux; // Colocar esse como uma variavel local da funcao...
	map <int,float> VertexCls;

	typedef struct eleMatriz{
		int ID;
 		float x1,x2,x3,y1,y2,y3,z1,z2,z3;
 		float Cl1, Cl2, Cl3;
		pVertex Pt1, Pt2, Pt3;



	} EleMatriz;

	class operadorEleMatriz {
	public:
		bool operator()(const EleMatriz a, const EleMatriz b) {
			if(a.ID < b.ID) {
				return true;
			}
			return false;
		}
	};

	

	void ClReader (char *CLS);
	void SaveBGMView1(double Cl,float GE); // Rotina para salvar view da malha de background
	void SaveBGMView2(float GE); // Rotina para salvar view da malha de background
	void Merge_Edges(pMesh theMesh2);
	void ResetMRE();
	void removeinternalfaces (pMesh theMesh);
	void removetetras (pMesh theMesh);
	int Iterador (set <int> ListofElements, int I, pMesh theMesh);
	void BoundaryElements3D (set <int> ListofElements, pMesh theMesh);
	void BoundaryElements2D (set <int> ListofElements, double Cl,int VerCL, pMesh theMesh);
	void removeinternaledges (pMesh theMesh);
	void removeinternalnodes (pMesh theMesh);
	void removeexternaledges(pMesh theMesh);
	void removeexternalnodes(pMesh theMesh);
	void RemoveStgElements_2D(pMesh theMesh, float GE);
	void RemoveStgElements_3D(pMesh theMesh);
	void removestgNodes(pMesh theMesh);
	
	void CriadordeLista(pMesh theMesh); // lembrar de retirar esse daqui depois
	void CriadordeLista3D(pMesh theMesh); // lembrar de retirar esse daqui depois
	int identifysurfaces_2D(pMesh Mesh);
	int Identify_and_Remove_Edges_2D(pMesh theMesh);
	int identifyvolumes_3D(pMesh Mesh);

set <EleMatriz, operadorEleMatriz> STriangles;
	eleMatriz * MtxTriangles; // As coordenadas e cls dos triangulos que formam a malha de background
	set <int> HoleNeighbours(pMesh theMesh2);
	set <pFace> ReturnBoundaryEdges();
	set <pFace> ReturnBoundaryNodes();
	set <pFace> ReturnBoundaryFaces();
	


private:
	///char meshFilename[];
	
	
};


