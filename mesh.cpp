#include "mesh.h"

namespace MeshDB{

	Mesh::Mesh(){}
	Mesh::~Mesh(){}

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

	void Mesh::getVertex(int ID, VertexInfo** vinfo){
		std::map<int, VertexInfo*>::iterator iter = VertexDB.find(ID);
		if (iter==VertexDB.end()){
			cout << "WARNING: Vertex does not exist.\n";
			vinfo = 0;
		}
		*vinfo = iter->second;
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

	void Mesh::setVertex(int ID, int physical, int geom){
		std::map<int, VertexInfo*>::iterator iter = VertexDB.find(ID);
		if (iter==VertexDB.end()){
			cout << "WARNING: vertex not found!\n";
		}
		else{
			VertexInfo *vinfo = iter->second;
			vinfo->geom = geom;
			vinfo->physical = physical;
		}
	}

	int Mesh::getNumVertices() const{
		return  (int)VertexDB.size();
	}

	int Mesh::getNumEdges() const{
		int k = 0;
		std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
		for(iter1=EdgeDB.begin(); iter1!=EdgeDB.end(); iter1++){
			k += (int)iter1->second.size();
		}
		return k;
	}

	int Mesh::getNumTriangles() const{
		return (int)TriangleDB.size();
	}

	int Mesh::getNumTetras() const{
		return (int)TetraDB.size();
	}

	void Mesh::getEdge(int id0, int id1, EdgeInfo** einfo){
		if (id0>id1){
			swap(id0,id1);
		}
		std::map<int, std::map<int, EdgeInfo*> >::iterator iter1;
		iter1 = EdgeDB.find(id0);
		if (iter1==EdgeDB.end()){
			cout << "Edge not found!\n";
		}

		std::map<int,EdgeInfo*> map2 = iter1->second;
		std::map<int,EdgeInfo*>::iterator iter2 = map2.find(id1);
		if (iter2 == map2.end()){
			cout << "Edge not found!\n";
		}
		*einfo = iter2->second;
	}

	void Mesh::findEdge(int id0, int id1, bool& found){
		found = false;
		if (id0>id1){
			swap(id0,id1);
		}
		std::map<int, std::map<int, EdgeInfo*> >::iterator iter1;
		iter1 = EdgeDB.find(id0);			// find first edge vertex ID
		if (iter1 != EdgeDB.end()){		// if found, find second edge vertex ID
			std::map<int,EdgeInfo*>::iterator iter2 = iter1->second.find(id1);
			if (iter2!=iter1->second.end()){
				found = true;
			}
		}
	}

	void Mesh::findEdge(int id0, int id1, bool& found0, bool& found1, std::map<int, std::map<int, EdgeInfo*> >::iterator& iter_out){
		iter_out = EdgeDB.find(id0);			// find first edge vertex ID
		if (iter_out != EdgeDB.end()){			// if found, find second edge vertex ID
			found0 = true;
			std::map<int,EdgeInfo*>::iterator iter = iter_out->second.find(id1);
			found1 = (iter!=iter_out->second.end())?true:false;
		}
		else{
			found0 = false;
			found1 = false;
		}
	}

	void Mesh::printMeshStatistic() const{
		ofstream fid;
		fid.open("MeshStatistic.txt");
		fid << "Mesh statistics:\n";
		fid << "Number of vertices           : " << getNumVertices() << endl;
		fid << "Number of edges over geometry: " << getNumEdges() << endl;
		fid << "Number of triangles          : " << getNumTriangles() << endl;
		fid << "Number of tetrahedra         : " << this->getNumTetras() << endl;
//		printVertexList(fid);
//		printEdgeList(fid);
//		printTriangleList(fid);
//		cout << "\n";
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

	void Mesh::createEdgeDataStructure(){
		//cout << "Creaing edge data structure:\n";
		if (getDim()==2){
			if (this->getElemType()==TRI){
				std::list<TriInfo*>::iterator iter = TriangleDB.begin();
				TriInfo* tinfo;
				for(;iter!=TriangleDB.end();iter++){
					tinfo = *iter;
					createEdge(tinfo->id0,tinfo->id1,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id0,tinfo->id2,tinfo->physical,tinfo->geom);
					createEdge(tinfo->id1,tinfo->id2,tinfo->physical,tinfo->geom);
				}
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

	int Mesh::getDim() const{
		return dim;
	}

	void Mesh::setDim(int d){
		dim = d;
	}

	ELEM_TYPE Mesh::getElemType() const{
		return elem_type;
	}

	void Mesh::setElemType(ELEM_TYPE et){
		elem_type = et;
	}

	void Mesh::setChacteristics(){
		setDim(2);
		setElemType(TRI);
		if (TetraDB.size()){
			setDim(3);
			setElemType(TETRA);
		}
	}

	void Mesh::getTetraVerticesCoords(TetraInfo* tinfo, Coords** p1, Coords** p2, Coords** p3, Coords** p4){
		VertexInfo* V1=0;
		VertexInfo* V2=0;
		VertexInfo* V3=0;
		VertexInfo* V4=0;
		getVertex(tinfo->id0,&V1);
		getVertex(tinfo->id1,&V2);
		getVertex(tinfo->id2,&V3);
		getVertex(tinfo->id3,&V4);
		*p1 = V1->coords;
		*p2 = V2->coords;
		*p3 = V3->coords;
		*p4 = V4->coords;
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
			cout << "tetra volume [" << tinfo->id0 << " " << tinfo->id1 << " " << tinfo->id2 << " " << tinfo->id3 << "] = " << v << endl;
		}
		tinfo=0;
		cout << "Tetrahedral mesh volume: " << v <<  endl;
	}
}



















