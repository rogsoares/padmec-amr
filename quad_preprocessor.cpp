/*
 * refine_TRI.cpp
 *
 *  Created on: 02/01/2015
 *      Author: rogerio
 */

#include "mesh.h"
/*
namespace MeshDB{
void  Mesh::quad_preprocessor(){

	// inicializar Cij, Dij atraves de um loop sobre as arestas
	// ---------------------------------------------------------------------
	int id0, k = 0;
	std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
	std::map<int, EdgeInfo*>::const_iterator iter2;
	fid << "Edges:\n";
	for(iter1 = EdgeDB.begin();iter1!=EdgeDB.end();iter1++){
		id0 = iter1->first;
		iter2 = iter1->second.begin();
		for(;iter2!=iter1->second.end();iter2++){
			EdgeInfo* einfo = iter2->second;
			einfo->Cij = new double[3];
			einfo->Cij[0] = einfo->Cij[1] = einfo->Cij[2] = 0;

			if (einfo->physical == 2000){
				einfo->Dij = new double[3];
				einfo->Dij[0] = einfo->Dij[1] = einfo->Dij[2] = 0;
			}

		}
	}
	// fim
	// --------------------------------------------------------------------------





	VertexInfo* V1=0;
	VertexInfo* V2=0;

	EdgeInfo* einfo;
	QuadInfo* tinfo=0;
	QuadInfo* tinfo_new=0;

	int I,J,K,L;

	std::map<int, VertexInfo*>::const_iterator VIter = VertexDB.end();
	int max_vertex_ID = VIter->first;

	std::list<QuadInfo*>::iterator iter1 = QuadDB.begin();
	int elem_counter = 0;

	// loop sobre os elementos da malha
	while (elem_counter < getNumQuad(BEFORE)){
		tinfo = *iter1;
		I = tinfo->id0;
		J = tinfo->id1;
		K = tinfo->id2;
		L = tinfo->id3;

		// I-J
		// J-K
		// K-L
		// L-I

		// get triangle edges:
		getEdge(I,J,&einfo);



	}
}

}

*/
