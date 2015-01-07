/*
 * refine_TRI.cpp
 *
 *  Created on: 02/01/2015
 *      Author: rogerio
 */

#include "mesh.h"

namespace MeshDB{
void  Mesh::refine_TRI(){
		VertexInfo* V1=0;
		VertexInfo* V2=0;
		VertexInfo* V_new=0;
		EdgeInfo* einfo;
		TriInfo* tinfo=0;
		TriInfo* tinfo_new=0;

		int I,J,K,R,S,T;

		std::map<int, VertexInfo*>::const_iterator VIter = VertexDB.end();
		int max_vertex_ID = VIter->first;

		std::list<TriInfo*>::iterator iter1 = TriangleDB.begin();
		int nelem = getNumTriangles(AFTER);
		int elem_counter = 0;
		while (elem_counter < nelem){
			tinfo = *iter1;
			I = tinfo->id0;
			J = tinfo->id1;
			K = tinfo->id2;

			// get triangle edges:
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
			}
			R = einfo->MPV_id;

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
			}
			T = einfo->MPV_id;

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
			}
			S = einfo->MPV_id;


			// create element: I-R-T, S-R-T, R-S-J, S-T-K
			tinfo_new = new TriInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = I;
			tinfo_new->id1 = R;
			tinfo_new->id2 = T;
			TriangleDB.push_back(tinfo_new);

			tinfo_new = new TriInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = S;
			tinfo_new->id1 = R;
			tinfo_new->id2 = T;
			TriangleDB.push_back(tinfo_new);

			tinfo_new = new TriInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = R;
			tinfo_new->id1 = S;
			tinfo_new->id2 = J;
			TriangleDB.push_back(tinfo_new);

			tinfo_new = new TriInfo;
			tinfo_new->geom = tinfo->geom;
			tinfo_new->physical = tinfo->physical;
			tinfo_new->id0 = S;
			tinfo_new->id1 = T;
			tinfo_new->id2 = K;
			TriangleDB.push_back(tinfo_new);

			tinfo=0;
			iter1 = TriangleDB.erase(iter1);
			elem_counter++;
		}
	}
}
