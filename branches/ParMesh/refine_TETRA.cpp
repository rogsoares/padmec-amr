/*
 * refine_TETRA.cpp
 *
 *  Created on: 02/01/2015
 *      Author: rogerio
 */

#include "mesh.h"

namespace MeshDB{

	void Mesh::refine_TETRA(){
		VertexInfo* V1=0;
		VertexInfo* V2=0;
		VertexInfo* V_new=0;
		EdgeInfo* einfo;

		TetraInfo* tinfo=0;
		TetraInfo* tinfo_new=0;

		int I,J,K,L,R,S,T,U,V,X;

		std::map<int, VertexInfo*>::const_iterator VIter = VertexDB.end();
		int max_vertex_ID = VIter->first;

		std::list<TetraInfo*>::iterator iter1 = TetraDB.begin();
		int nelem = getNumTetras();
		int elem_counter = 0;
		while (elem_counter < nelem){
			tinfo = *iter1;
			I = tinfo->id0;
			J = tinfo->id1;
			K = tinfo->id2;
			L = tinfo->id3;

			// get tetrahedron edges:
			getEdge(I,J,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(I,&V1);
				this->getVertex(J,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			S = einfo->MPV_id;

			getEdge(I,K,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(I,&V1);
				this->getVertex(K,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			R = einfo->MPV_id;

			getEdge(J,K,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(J,&V1);
				this->getVertex(K,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			V = einfo->MPV_id;

			getEdge(I,L,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(I,&V1);
				this->getVertex(L,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			T = einfo->MPV_id;

			getEdge(J,L,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(J,&V1);
				this->getVertex(L,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			X = einfo->MPV_id;

			getEdge(L,K,&einfo);
			if (einfo->MPV_id < 0){
				++max_vertex_ID;
				this->getVertex(L,&V1);
				this->getVertex(K,&V2);
				V_new = new VertexInfo;
				V_new->coords = new Coords[3];
				middlePoint(V1->coords,V2->coords,V_new->coords);
				V_new->geom = einfo->geom;
				V_new->physical = einfo->physical;
				this->createVertex(max_vertex_ID,V_new);
				einfo->MPV_id = max_vertex_ID;
				V1 = V2 = 0;
			}
			U = einfo->MPV_id;

			// create element: I-R-T-S | K-R-S-V | L-T-X-U | J-V-X-S | R-S-U-V | S-T-X-U | R-T-U-S | V-U-X-S
			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = I;
			tinfo_new->id1 = R;
			tinfo_new->id2 = T;
			tinfo_new->id3 = S;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = K;
			tinfo_new->id1 = R;
			tinfo_new->id2 = U;
			tinfo_new->id3 = V;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = L;
			tinfo_new->id1 = T;
			tinfo_new->id2 = U;
			tinfo_new->id3 = X;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = J;
			tinfo_new->id1 = V;
			tinfo_new->id2 = X;
			tinfo_new->id3 = S;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = R;
			tinfo_new->id1 = S;
			tinfo_new->id2 = U;
			tinfo_new->id3 = V;
			TetraDB.push_back(tinfo_new);

//			// create element: I-R-T-S | K-R-S-V | L-T-X-U | J-V-X-S | R-S-U-V | S-T-X-U | R-T-U-S | V-U-X-S
			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = S;
			tinfo_new->id1 = T;
			tinfo_new->id2 = X;
			tinfo_new->id3 = U;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = R;
			tinfo_new->id1 = T;
			tinfo_new->id2 = U;
			tinfo_new->id3 = S;
			TetraDB.push_back(tinfo_new);

			tinfo_new = new TetraInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = V;
			tinfo_new->id1 = U;
			tinfo_new->id2 = X;
			tinfo_new->id3 = S;
			TetraDB.push_back(tinfo_new);

			//tinfo=0;
			iter1 = TetraDB.erase(iter1);
			elem_counter++;
		}
	}
}


