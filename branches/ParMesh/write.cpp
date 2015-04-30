/*
 * write.cpp
 *
 *  Created on: 03/01/2015
 *      Author: rogerio
 */

#include "mesh.h"

namespace MeshDB{
	void Mesh::write(const char* filename){

		char fname[512];
		sprintf(fname,"%s-%d-of-%d.msh",filename,this->get_rank(),this->get_nproc());

		ofstream fid;
		fid.open(fname);

		int v_count = 0;
		int e_count = 0;
		int f_count = getNumTriangles(AFTER);
		int q_count = getNumQuad(AFTER);
		int t_count = getNumTetras(AFTER);


		// Vertices: ID - x - y - z
		// ------------------------------------------------------------------------
		fid << "$NOD\n" << getNumVertices(AFTER) << std::endl;
		VertexInfo* vinfo=0;
		std::map<int, VertexInfo*>::iterator vit = VertexDB.begin();
		std::map<int, VertexInfo*> VertexDB_aux;
		for(;vit!=VertexDB.end(); vit++){
			int ID = vit->first;
			vinfo = vit->second;
			fid << ID << " " << vinfo->coords[0] << " " << vinfo->coords[1] << " " << vinfo->coords[2] << endl;
			if (vinfo->geom != -1){
				VertexDB_aux[ID] = vinfo;
				v_count++;
			}
		}
		fid << "$ENDNOD\n";
		vinfo=0;

		// count how many edges are on geometry (boundary)
		std::list<int*> bedges;
		int* bedge_info=0;
		std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
		std::map<int, EdgeInfo*>::const_iterator iter2;
		for(iter1 = EdgeDB.begin();iter1!=EdgeDB.end();iter1++){
			int id0 = iter1->first;
			iter2 = iter1->second.begin();
			for(;iter2!=iter1->second.end();iter2++){
				int id1 = iter2->first;
				EdgeInfo* einfo = iter2->second;
				if (einfo->geom != -1){
					bedge_info = new int[4];
					bedge_info[0] = id0;
					bedge_info[1] = id1;
					bedge_info[2] = einfo->physical;
					bedge_info[3] = einfo->geom;
					bedges.push_back(bedge_info);
					e_count++;
				}
			}
		}

		// number of element to be printed:
		fid << "$ELM\n";
		fid << v_count + e_count + q_count + f_count + t_count << endl;

		int count = 0;

		// print vertices flagged over geometry
		vit = VertexDB_aux.begin();
		for(;vit!=VertexDB_aux.end(); vit++){
			int ID = vit->first;
			vinfo = vit->second;
			// 1 15 10 1 1 1
			fid << ++count << " 15 " << vinfo->physical << " " << vinfo->geom << " 1 " << ID << endl;
		}

		// print edges flagged over geometry
		std::list<int*>::iterator eit = bedges.begin();
		for(;eit!=bedges.end();eit++){
			bedge_info = *eit;
			//fid << ++count << " 1 " << bedge_info[2] << " " << bedge_info[3] << " 2 " << bedge_info[0] << " " << bedge_info[1] << endl;
		}

		// print triangles over geometry (if 3-D)
		TriInfo* finfo=0;
		std::list<TriInfo*>::iterator fit = TriangleDB.begin();
		for(;fit!=TriangleDB.end();fit++){
			finfo = *fit;
			fid << ++count << " 2 " << finfo->physical << " " << finfo->geom << " 3 " << finfo->id0 << " " << finfo->id1 << " " << finfo->id2 << endl;
		}

		// print tetrahedra
		TetraInfo* tinfo=0;
		std::list<TetraInfo*>::iterator tit = TetraDB.begin();
		for(;tit!=TetraDB.end();tit++){
			tinfo = *tit;
			fid << ++count << " 4 " << tinfo->physical << " " << tinfo->geom << " 4 " << tinfo->id0 << " " << tinfo->id1 << " " << tinfo->id2 << " " << tinfo->id3 << endl;
		}
		fid << "$ENDELM\n";
		fid.close();
	}

	void Mesh::printMeshStatistic(const char* filename) const{
		char fname[512];
		sprintf(fname,"MeshStatistic__nproc-%d",get_nproc());



		ofstream fid;
		fid.open(fname);

		fid << "Mesh statistics:\n";
		fid << "-----------------------------------------------------------\n\n";
		fid << "Number of refinement levels: " << _refLevel << endl;
		fid << "Mesh dimension             : " << this->getDim() << endl;
		fid << "Mesh element type          : " << getElementType() << endl;
		fid << "\n                               Before    After\n";
		fid << "Number of vertices           : " << getNumVertices(BEFORE) << "      " << getNumVertices(AFTER) << endl;
		fid << "Number of edges over geometry: " << getNumEdges(BEFORE) << "      " << getNumEdges(AFTER) << endl;
		fid << "Number of triangles          : " << getNumTriangles(BEFORE) << "      " << getNumTriangles(AFTER) << endl;
		fid << "Number of quad               : " << getNumQuad(BEFORE) << "      " << getNumQuad(AFTER) << endl;
		fid << "Number of tetrahedra         : " << getNumTetras(BEFORE) << "      " << getNumTetras(AFTER) << endl;

	}

	void Mesh::printVertexList(ofstream& fid) const{
		fid << "Vertices:\n";
		VertexInfo *vinfo;
		std::map<int, VertexInfo*>::const_iterator iter = VertexDB.begin();
		for(;iter!=VertexDB.end();iter++){
			vinfo = iter->second;
			fid << iter->first << " " << "coords: " << vinfo->coords[0] << " " << vinfo->coords[1] << " " << vinfo->coords[2] << " G:" << vinfo->geom << " P:" << vinfo->physical << "\n";
		}
	}

	void Mesh::printEdgeList(ofstream& fid) const{
		int id0, k = 0;
		std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
		std::map<int, EdgeInfo*>::const_iterator iter2;
		fid << "Edges:\n";
		for(iter1 = EdgeDB.begin();iter1!=EdgeDB.end();iter1++){
			id0 = iter1->first;
			iter2 = iter1->second.begin();
			for(;iter2!=iter1->second.end();iter2++){
				EdgeInfo* einfo = iter2->second;
				fid << ++k << ": " << id0 << " - " << iter2->first << " G:" << einfo->geom << " P:" << einfo->physical << endl;
			}
		}
	}

	void Mesh::printTriangleList(ofstream& fid) const{
		int k = 0;
		fid << "Triangles:\n";
		std::list<TriInfo*>::const_iterator iter = TriangleDB.begin();
		TriInfo* tinfo;
		for(;iter!=TriangleDB.end();iter++){
			tinfo = *iter;
			fid << ++k << ": " << tinfo->id0 << " " << tinfo->id1 << " " << tinfo->id2 << " " << tinfo->geom << " " << tinfo->physical << endl;
		}
	}
}
