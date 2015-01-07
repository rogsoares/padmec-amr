/*
 * refine_TETRA.cpp
 *
 *  Created on: 02/01/2015
 *      Author: rogerio
 */

#include "mesh.h"

namespace MeshDB{

	void Mesh::refine_TETRA(){
		TetraInfo* tinfo=0;
		int I,J,K,L,R,S,T,U,V,X;

		std::map<int, VertexInfo*>::const_iterator VIter = VertexDB.end();
		int max_vertex_ID = VIter->first;

		std::list<TetraInfo*>::iterator iter1 = TetraDB.begin();
		int nelem = getNumTetras(AFTER);
		int elem_counter = 0;
		while (elem_counter < nelem){
			tinfo = *iter1;
			I = tinfo->id0;
			J = tinfo->id1;
			K = tinfo->id2;
			L = tinfo->id3;

			// get tetrahedron edges:
			getVertexID_EdgeTetra(I,J,max_vertex_ID,S);
			getVertexID_EdgeTetra(I,K,max_vertex_ID,R);
			getVertexID_EdgeTetra(J,K,max_vertex_ID,V);
			getVertexID_EdgeTetra(I,L,max_vertex_ID,T);
			getVertexID_EdgeTetra(J,L,max_vertex_ID,X);
			getVertexID_EdgeTetra(L,K,max_vertex_ID,U);

			// create element: I-R-T-S | K-R-S-V | L-T-X-U | J-V-X-S | R-S-U-V | S-T-X-U | R-T-U-S | V-U-X-S
			createTetrahedron(I,R,T,S,tinfo->physical,tinfo->geom);
			createTetrahedron(K,R,U,V,tinfo->physical,tinfo->geom);
			createTetrahedron(L,T,U,X,tinfo->physical,tinfo->geom);
			createTetrahedron(J,V,X,S,tinfo->physical,tinfo->geom);
			createTetrahedron(R,S,U,V,tinfo->physical,tinfo->geom);
			createTetrahedron(S,T,U,X,tinfo->physical,tinfo->geom);
			createTetrahedron(R,T,U,S,tinfo->physical,tinfo->geom);
			createTetrahedron(V,U,X,S,tinfo->physical,tinfo->geom);

			iter1 = TetraDB.erase(iter1);
			elem_counter++;
		}
	}

	void Mesh::getVertexID_EdgeTetra(int id0, int id1, int& max_vertex_ID, int& ID){
		VertexInfo* V1=0;
		VertexInfo* V2=0;
		VertexInfo* V_new=0;
		EdgeInfo* einfo;

		getEdge(id0,id1,&einfo);
		if (einfo->MPV_id < 0){
			++max_vertex_ID;
			getVertex(id0,&V1);
			getVertex(id1,&V2);
			V_new = new VertexInfo;
			V_new->coords = new Coords[3];
			middlePoint(V1->coords,V2->coords,V_new->coords);
			V_new->geom = einfo->geom;
			V_new->physical = einfo->physical;
			createVertex(max_vertex_ID,V_new);
			einfo->MPV_id = max_vertex_ID;
			V1 = V2 = 0;
		}
		ID = einfo->MPV_id;
	}
}


