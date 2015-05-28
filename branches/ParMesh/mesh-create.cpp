#include "mesh.h"

namespace MeshDB{

	void Mesh::createVertex(int ID, VertexInfo* vinfo){
		VertexDB[ID] = vinfo;
	}

	void Mesh::createVertex(int ID, double x, double y, double z, int physical, int geom){
		VertexInfo *vinfo = new VertexInfo;
		vinfo->coords = new Coords[3];
		vinfo->coords[0] = x;
		vinfo->coords[1] = y;
		vinfo->coords[2] = z;
		vinfo->geom = geom;
		vinfo->physical = physical;
		VertexDB[ID] = vinfo;
	}

	void Mesh::createEdge(int id0, int id1, int physical, int geom){
		EdgeInfo* einfo = new EdgeInfo;
		einfo->geom = geom;
		einfo->physical = physical;
		einfo->MPV_id = -1;
		std::map<int,EdgeInfo*> secondPart;
		std::map<int, std::map<int, EdgeInfo*> >::iterator iter;
		bool found0, found1;

		if (id0 > id1){
			swap(id0,id1);
		}
		findEdge(id0,id1,found0,found1,iter);
		if (found0 && !found1){
			iter->second.insert( std::pair<int,EdgeInfo*>(id1,einfo) );
		}
		else if (!found0 && !found1){
			secondPart.insert( std::pair<int,EdgeInfo*>(id1,einfo) );
			EdgeDB[id0] = secondPart;
		}
	}

	void Mesh::createTriangle(int id0, int id1, int id2, int physical, int geom){
		TriInfo* tinfo = new TriInfo;
		tinfo->id0 = id0;
		tinfo->id1 = id1;
		tinfo->id2 = id2;
		tinfo->geom = geom;
		tinfo->physical = physical;
		TriangleDB.push_back(tinfo);
	}

	void Mesh::createQuad(int id0, int id1, int id2, int id3,int physical, int geom){
		QuadInfo* qinfo = new QuadInfo;
		qinfo->id0 = id0;
		qinfo->id1 = id1;
		qinfo->id2 = id2;
		qinfo->id3 = id3;
		qinfo->geom = geom;
		qinfo->physical = physical;
		QuadDB.push_back(qinfo);
	}

	void Mesh::createTetrahedron(int id0, int id1, int id2, int id3,int physical, int geom){
		TetraInfo* tinfo = new TetraInfo;
		tinfo->id0 = id0;
		tinfo->id1 = id1;
		tinfo->id2 = id2;
		tinfo->id3 = id3;
		tinfo->geom = geom;
		tinfo->physical = physical;
		TetraDB.push_back(tinfo);
	}

	void Mesh::createEdgeDataStructure(){
		cout << "Creating edge data structure:\n";
		if (getDim()==2){
			switch ( getElemType() ){
			case TRI:{
				std::list<TriInfo*>::iterator iter = TriangleDB.begin();
				TriInfo* tinfo;
				for(;iter!=TriangleDB.end();iter++){
					tinfo = *iter;
					createEdge(tinfo->id0,tinfo->id1,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id0,tinfo->id2,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id1,tinfo->id2,tinfo->physical,tinfo->geom);
				}
				break;
			}
			case QUAD:{
				std::list<QuadInfo*>::iterator iter = QuadDB.begin();
				QuadInfo* qinfo;
				for(;iter!=QuadDB.end();iter++){
					qinfo = *iter;
					createEdge(qinfo->id0,qinfo->id1,qinfo->physical,qinfo->geom);
					createEdge(qinfo->id1,qinfo->id2,qinfo->physical,qinfo->geom);
					createEdge(qinfo->id2,qinfo->id3,qinfo->physical,qinfo->geom);
					createEdge(qinfo->id3,qinfo->id0,qinfo->physical,qinfo->geom);
				}
				break;
			}
			default:
				cout << "Error! Exiting....\n";
				exit(1);
			}
		}
		else{
			if (this->getElemType()==TETRA){
				std::list<TetraInfo*>::iterator iter = TetraDB.begin();
				TetraInfo* tinfo;
				for(;iter!=TetraDB.end();iter++){
					tinfo = *iter;
					createEdge(tinfo->id0,tinfo->id1,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id0,tinfo->id2,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id0,tinfo->id3,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id1,tinfo->id2,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id1,tinfo->id3,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id2,tinfo->id3,tinfo->physical,tinfo->geom);
				}
			}
		}

		if (!numEdges_before){
			std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
			for(iter1=EdgeDB.begin(); iter1!=EdgeDB.end(); iter1++){
				numEdges_before += (int)iter1->second.size();
			}
		}
	}

	void Mesh::deleteEdgeDataStructure(){
		std::map<int, std::map<int, EdgeInfo*> >::iterator iter1;
		std::map<int, EdgeInfo*>::iterator iter2;
		for(iter1 = EdgeDB.begin();iter1!=EdgeDB.end();iter1++){
			iter2 = iter1->second.begin();
			for(;iter2!=iter1->second.end();iter2++){
				delete iter2->second;
			}
			iter1->second.clear();
		}
		EdgeDB.clear();
		//cout << "Edge Data Structure Deleted:\n";
	}

	void Mesh::calculate_volume(){
		Coords* p1=0;
		Coords* p2=0;
		Coords* p3=0;
		Coords* p4=0;
		std::list<TetraInfo*>::iterator iter = TetraDB.begin();
		double v = 0;
		TetraInfo* tinfo=0;
		for(;iter!=TetraDB.end();iter++){
			tinfo = *iter;
			getTetraVerticesCoords(tinfo,&p1,&p2,&p3,&p4);
			v += R_Volume(p1,p2,p3,p4);
			//cout << "tetra volume [" << tinfo->id0 << " " << tinfo->id1 << " " << tinfo->id2 << " " << tinfo->id3 << "] = " << v << endl;
		}
		tinfo=0;
		//cout << "Tetrahedral mesh volume: " << v <<  endl;
	}
}
