#include "mesh.h"

namespace MeshDB{

	Mesh::Mesh(){
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);	/* get current process id */
		MPI_Comm_size(MPI_COMM_WORLD, &nproc);	/* get number of processes */
	}
	Mesh::~Mesh(){
	}

	void Mesh::getVertex(int ID, VertexInfo** vinfo){
		std::map<int, VertexInfo*>::iterator iter = VertexDB.find(ID);
		if (iter==VertexDB.end()){
			cout << "WARNING: Vertex does not exist.\n";
			vinfo = 0;
		}
		*vinfo = iter->second;
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

	int Mesh::getNumVertices(REF_MOMENT rm) const{
		int n;
		switch(rm){
		case BEFORE:
			n = numVertices_before;
			break;
		case AFTER:
			n = (int)VertexDB.size();
			break;
		}
		return n;
	}

	int Mesh::getNumEdges(REF_MOMENT rm) const{
		int n;
		switch(rm){
		case BEFORE:
			n = numEdges_before;
			break;
		case AFTER:
			n = 0;
			std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;
			for(iter1=EdgeDB.begin(); iter1!=EdgeDB.end(); iter1++){
				n += (int)iter1->second.size();
			}
			break;
		}
		return n;
	}

	int Mesh::getNumTriangles(REF_MOMENT rm) const{
		int n;
		switch(rm){
		case BEFORE:
			n = numTriangles_before;
			break;
		case AFTER:
			n = (int)TriangleDB.size();
			break;
		}
		return n;
	}

	int Mesh::getNumQuad(REF_MOMENT rm) const{
		int n;
		switch(rm){
		case BEFORE:
			n = numQuad_before;
			break;
		case AFTER:
			n = (int)QuadDB.size();
			break;
		}
		return n;
	}

	int Mesh::getNumTetras(REF_MOMENT rm) const{
		int n;
		switch(rm){
		case BEFORE:
			n = numTetra_before;
			break;
		case AFTER:
			n = (int)TetraDB.size();
			break;
		}
		return n;
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

	int Mesh::getDim() const{
		return dim;
	}

	void Mesh::setDim(int eType){
		dim = 2;
		if (eType>3){
			dim = 3;
		}
	}

	ELEM_TYPE Mesh::getElemType() const{
		return elem_type;
	}

	string Mesh::getElementType() const{
		string etype;
		switch(elem_type){
		case TRI:
			etype = "Triangles";
			break;
		case QUAD:
			etype = "Quad";
			break;
		case TETRA:
			etype = "Tetra";
		}
		return etype;
	}

	void Mesh::setElemType(int eType){
		switch (eType){
		case 2:
			elem_type = TRI;
			break;
		case 3:
			elem_type = QUAD;
			break;
		case 4:
			elem_type = TETRA;
			break;
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
}
