/*
 * mesh.h
 *
 *  Created on: 01/01/2015
 *      Author: rogerio
 */

#ifndef MESH_H_
#define MESH_H_

#include "utilities.h"

typedef double Coords;

struct VertexInfo{
	Coords* coords;
	int physical;
	int geom;
	int numRCopies;
};

struct EdgeInfo{
	int MPV_id;		// Middle Point Vertex ID: for refinement purposes
	int physical;	// physical flag: dirichlet, neumann, etc
	int geom;		// geometry identification
	int numRCopies;
};

struct TriInfo{
	int id0;
	int id1;
	int id2;
	int geom;
	int physical;
	int numRCopies;
};

struct QuadInfo{
	int id0;
	int id1;
	int id2;
	int id3;
	int geom;
	int physical;
};

struct TetraInfo{
	int id0;
	int id1;
	int id2;
	int id3;
	int geom;
	int physical;
};

enum ELEM_TYPE{TRI, QUAD, TETRA};
enum REF_MOMENT{BEFORE, AFTER};

namespace MeshDB{

	class Mesh{
	public:
		Mesh();
		~Mesh();

		void createVertex(int ID, VertexInfo* vinfo);
		void createVertex(int ID, double x, double y, double z, int physical=0, int geom=0);
		void getVertex(int ID,VertexInfo** vinfo);
		void createEdge(int id0, int id1, int physical, int geom);
		void createEdgeDataStructure();
		void deleteEdgeDataStructure();
		void createTriangle(int id0, int id1, int id2, int physical, int geom);
		void setVertex(int ID, int physical, int geom);
		void createTetrahedron(int id0, int id1, int id2, int id3, int physical, int geom);

		int getNumVertices(REF_MOMENT rm) const;
		int getNumEdges(REF_MOMENT rm) const;
		int getNumTriangles(REF_MOMENT rm) const;
		int getNumQuad() const { return 0; }
		int getNumTetras(REF_MOMENT rm) const;

		void refine_mesh(int refLevel);

		void getEdge(int id0, int id1, EdgeInfo**);
		void findEdge(int id0, int id1, bool& found);
		void findEdge(int id0, int id1, bool& found0, bool& found1, std::map<int, std::map<int, EdgeInfo*> >::iterator& iter);
		void printMeshStatistic(const char* filename) const;

		int getDim() const;
		void setDim(int d);

		ELEM_TYPE getElemType() const;
		string getElementType() const;
		void setElemType(ELEM_TYPE et);

		void setChacteristics();

		void read(const char* filename);
		void write(const char* filename);

	private:

		int dim;
		ELEM_TYPE elem_type;

		int numVertices_before;
		int numVertices_after;
		int numEdges_before;
		int numTriangles_before;
		int numQuad_before;
		int numTetra_before;

		int _refLevel;

		// where vertices, edges and triangles are stored
		std::map<int, VertexInfo*> VertexDB;
		std::map<int, std::map<int, EdgeInfo*> >  EdgeDB;
		std::list<TriInfo*> TriangleDB;
		std::list<TetraInfo*> TetraDB;

		void printVertexList(ofstream& fid) const;
		void printEdgeList(ofstream& fid) const;
		void printTriangleList(ofstream& fid) const;

		void refine_TRI();
		void refine_QUAD();
		void refine_TETRA();

		void getTetraVerticesCoords(TetraInfo* tinfo, Coords** p1, Coords** p2, Coords** p3, Coords** p4);
		void getVertexID_EdgeTetra(int id0, int id1, int& max_vertex_ID, int& ID);
		void calculate_volume();

		// PARALLEL STUFFES
		// ---------------------------------------------------------------------------
		int rank;
		int get_rank() const;

		int nproc;
		int get_nproc() const;

		// serial mesh partition into nproc parts
		void mesh_partition();

		// migrate mesh entities from rank 0 to all other processes
		void mesh_distribution();

		// create consistent global vertex ID over partition's boundaries.
		void bdry_linkSetup();
	};

}

void readMesh(const char* filename, MeshDB::Mesh* mesh);
void middlePoint(const Coords* p1, const Coords* p2, Coords* mp);
void middlePoint(const Coords* p1, const Coords* p2, const Coords* p3, const Coords* p4, Coords* mp);

/*
 * MeshDB: uma malha de elementos finitos pode ser armazenada atraves da estrutura Mesh (Dentro do namespace MeshDB).
 * O comando
 *
 * 1) read: ler on nós (ID, x, y, z) e conectiviadades dos elementos (apenas triangulos, por enquanto) assim como
 *    os flags definidos sobre a geometria e que ficam associados as entidade da malha apos sua gereção.
 *
 * 2) refine_mesh: subdivide os triangulos usando os pontos médio das suas arestas. Assim, um triangulo gera quatro novos
 *    triangulos. Para que este refinamento ocorra, é preciso criar uma estrutura de arestas (lista de arestas) sobre onde
 *    serão criados os novos nós da malha refinada. Partindo de uma malha M_k, o processo de refinamento é feito elemento-
 *    a-elemento e se excerra quando o número de elementos da malha M_k é atingido quando então, surge uma nova malha M_k+1.
 *    Ao término, a malha M_k é eliminada da memória da seguinte forma: para cada 4 novos elementos filhos gerados, o elemento
 *    pai é automaticamente removida da estrutura. As arestas de M_k são removidas no final do processo.
 *
 * Estrutura de arestas:
 *
 * A criação da estrutura de arestas é criada usando uma combinação do comando std::map e é baseada no ID global dos nós
 * que definem a aresta. O primeiro map define um tipo chamado EdgeDB que será usado para acessar de fato as arestas da malha.
 * Este map é do tipo:
 *
 * 		std::map<edge_id0, std::map<edge_id1,EdgeInfo*>>
 *
 * onde:
 *
 *      edge_id0                    : refere-se a um numero inteiro e que representa o menor ID dos dois ID dos nós da aresta.
 *      std::map<edge_id1,EdgeInfo*>: map que armazena o segundo ID da aresta (o maior) e associa a esta aresta alguma informação
 *                                    a seu respeito.
 *
 *  Considere a malha abaixo:
 *  1-3-7
 *  7-4-3
 *  4-6-2
 *  2-5-4
 *  4-3-5
 *
 *  Suas arestas são:
 *  1-7
 *  1-3
 *  2-4
 *  2-5
 *  2-6
 *  3-4
 *  3-5
 *  3-7
 *  4-5
 *  4-6
 *  4-7
 *
 * as quais podem ser agrupadas da seguinte forma:
 * 1 - 3,7
 * 2 - 4,5,6
 * 3 - 4,5,7
 * 4 - 5,6,7
 *
 * Ou seja, as arestas podem ser agrupadas obedecendo a seguinte ordem:
 *
 * ID_I - ID_j1, ID_j2, ... , ID_jN
 */




#endif /* MESH_H_ */
