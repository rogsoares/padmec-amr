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

	void Mesh::verify_coeff_calculation(){

	/*
	crie uma funcao verify_coeff_calculation() que deve ser chamada apenas ao fim do pre-processador.
	 So depois que voce calcular tudo é que voce verifica se tudo está OK.

Verificar:

1) Cij
2) Dij
3) volumes

                               J1________
                               |        |    Ci-> element centroid...
                           C2..|....C1  |    IJi-> conection edge
                      J2___:___I___:____J4   Here, I is a internal node
                           :   |   :    |
                           C3..|...:C4  |
                               |        |
                               J3_______|

A soma dos Cij (e Dij se o nó for de contorno) em torno de um VOLUME DE CONTROLE deve ser zero. Ou todos os vetores apontam para dentro ou todos
apontam para fora do volume de controle. Se há 100 nós, haverá 100 somatórios.

Esse somatorio deve ser feito POR ARESTA. Inicialize um somador para cada nó que receberá o Cij para o nó I da aresta e
-Cij para o nó J desta aresta. Ao final do LOOP NAS ARESTAS o somatório para cada nó terá sido feito. Para um nó entre dois dominio
haverá dois somadores (não esquecer de incluir od Dij!)

A soma de todos os volume de controle deve ser igual a ao “volume” dos dominos.

Exemplo:  se o dominio é um quadrado de lado 1 entao a soma dos volumes de controle deve ser 1 também.

Exemplo: se o dominio é um quadrado dividivo ao meio, entao a soma da será 1/2 para ambos os domínios.
	*/
        cout<<"verifying coefficient calculation...\n\n";
        double CijSum[3]={0,0,0},DijSum[3]={0,0,0};
        int id0, k = 0;

        std::map<int, std::map<int, EdgeInfo*> >::const_iterator iter1;  //associado em Mesh : std::map<int, std::map<int, EdgeInfo*> >  EdgeDB;
        std::map<int, EdgeInfo*>::const_iterator iter2;  //iterador para o valor mapeado de iter1, percorre arestas...

        //nesse loop, percorre-se o banco de dados de arestas
        for(iter1 = EdgeDB.begin();iter1!=EdgeDB.end();iter1++){ //iterator aponta para o primeiro elemento de EdgeDB, e percorre cada elemento, até que chegue ao ultimo
            id0 = iter1->first; //iter1->first : valor-chave de EdgeDB, apontado por iter1. Nesse caso, o valor de indentificacao de um nó(inteiro) de aresta do elemento de EdgeDB ...
            iter2 = iter1->second.begin(); // aqui, iter2 aponta para a primeira aresta do elemento de indice iter1->first de EdgeDB
            for(;iter2!=iter1->second.end();iter2++){  //nesse loop, percorre-se as arestas de algum elemento...

                EdgeInfo* einfo = iter2->second;
                cout<<"\nCij ("<<einfo->physical<<") : "<<einfo->Cij[einfo->physical][0]<<","<<endl;


               // einfo->Cij[einfo->physical]=new double[3];
               // einfo->Cij[einfo->physical][0] = einfo->Cij[einfo->physical][1] = einfo->Cij[einfo->physical][2] = 0; //map <flag,vetor>

                 //   if (einfo->physical != 3300){ //dominio de contorno,etc...
                   //     einfo->Dij = new double[3];
                    //    einfo->Dij[0] = einfo->Dij[1] = einfo->Dij[2] = 0;
                   // }

            }
	}
	// fim
	// --------------------------------------------------------------------------




        /* calculate volume of control volume and associate it to elements nodes
           o volume dos volumes de controle correspondem as areas dos elementos
           atribui-se a cada nó de um elemento 1/4 de sua área((b1*h1)/2+(b2*h2)/2). Isso é feito de forma
           acumulativa a fim de que ao término do "loop" sobre os elementos da malha,
           todos os nós tenham seus devidos volumes calculados */

        ///+y
        /// ^
        /// |
        /// L_____K                  B   *   H          V[0]->I
        /// | A2 /|           A1=(||J-I||*||K-J||)/2    V[1]->J
        /// |  /A1|           A2=(||K-L||*||L-I||)/2    V[2]->K
        /// I/____J ---->+x                             V[3]->L

        double B,H,A1,A2,A;
       /* B=sqrt(pow((V[1]->coords[0]-V[0]->coords[0]),2)+pow((V[1]->coords[1]-V[0]->coords[1]),2));//b1
        H=sqrt(pow((V[2]->coords[0]-V[1]->coords[0]),2)+pow((V[2]->coords[1]-V[1]->coords[1]),2));//h1
        A1=(B*H)/2;
        B=sqrt(pow((V[2]->coords[0]-V[3]->coords[0]),2)+pow((V[2]->coords[1]-V[3]->coords[1]),2));//b2
        H=sqrt(pow((V[3]->coords[0]-V[0]->coords[0]),2)+pow((V[3]->coords[1]-V[0]->coords[1]),2));//h2
        A2=(B*H)/2;
        A=(A1+A2)/4.0;     // element area
        //cout << "area = " <<4*A<< endl;

        for (int k=0; k<4; k++){
            V[k]->volume += A;
        }*/


        cout<<"\nCij sum -> ("<<CijSum[0]<<","<<CijSum[1]<<","<<CijSum[2]<<")\n";
        cout<<"Dij sum -> ("<<DijSum[0]<<","<<DijSum[1]<<","<<DijSum[2]<<")\n";

	}
}
