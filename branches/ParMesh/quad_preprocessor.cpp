/*
 * quad_preprocessor.cpp
 *
 *  Created on: 28/04/2015
 *      Author: Rogerio Soares
 */

#include "mesh.h"
#include "math.h"


/* Gera arquivo de saida contendo os coeficientes geometricos Cij e Dij */

namespace MeshDB{
void  Mesh::quad_preprocessor(){

    double** p=new double*[getNumQuad(BEFORE)*4];//stores Cij's values

	VertexInfo* V[4]; //vertices coordinates
	EdgeInfo* einfo;
	QuadInfo* qinfo=0;
	QuadInfo* qinfo_new=0;
	double EdgeCG[3]={0,0,0};
    double CijSum[3]={0,0,0},DijSum[3]={0,0,0};

	int I,J,K,L;

	std::map<int, VertexInfo*>::const_iterator VIter = VertexDB.begin();
	int max_vertex_ID = VIter->first;

	std::list<QuadInfo*>::iterator iter3 = QuadDB.begin();
	int elem_counter,edge_counter = 0;

	// loop sobre os elementos da malha
	cout<<"\nNumber of mesh elements(quadrilaterals) : "<<getNumQuad(BEFORE);
	cout<<"\nLooping over mesh elements...\n";

	//while (elem_counter < getNumQuad(BEFORE)){
	for(;iter3!=QuadDB.end();iter3++){

		double*ElemCG=new double[3];
		double Cij[3],Dij[3];
		EdgeInfo* QuadEdgeAdress[4];
		int currentedge=0; //used to index pointer p

        qinfo = *iter3; //select current quad
        int dom=qinfo->physical;


		I = qinfo->id0;
		J = qinfo->id1;
		K = qinfo->id2;
		L = qinfo->id3;

		//get quadrilateral edges
		getEdge(I,J,&einfo);QuadEdgeAdress[0]=einfo;
        getEdge(J,K,&einfo);QuadEdgeAdress[1]=einfo;
        getEdge(K,L,&einfo);QuadEdgeAdress[2]=einfo;
        getEdge(L,I,&einfo);QuadEdgeAdress[3]=einfo;

        //get vertex coords
        getVertex(I,&V[0]); getVertex(J,&V[1]);
        getVertex(K,&V[2]); getVertex(L,&V[3]);


        //get element gravity center
         middlePoint(V[0]->coords, V[1]->coords, V[2]->coords, V[3]->coords, ElemCG);


        //cout<<"\nQUAD CG: ("<<ElemCG[0]<<","<<ElemCG[1]<<","<<ElemCG[2]<<")\n\n";

        int numBE = 0;
        //loop over quadrialateral element edges...
        for(int i=0;i<getNumQuad(BEFORE);i++){

            int j=0;
            int id0,id1=0;
            //adjust indexes
            if((i+1)<getNumQuad(BEFORE)){
                j=i+1;
            }
                else{
                    i=3;j=0;
                }

           switch(i){
                case 0:id0=I;break;
                case 1:id0=J;break;
                case 2:id0=K;break;
                case 3:id0=L;break;
            }
            switch(j){
                case 0:id1=I;break;
                case 1:id1=J;break;
                case 2:id1=K;break;
                case 3:id1=L;break;


            }

            //get edge centroid
            middlePoint(V[i]->coords,V[j]->coords,EdgeCG);
            //cout<<"\nEDGE CG: ("<<EdgeCG[0]<<","<<EdgeCG[1]<<","<<EdgeCG[2]<<")\n\n";

            // edge vector, used as a reference vector
            double IJ[2] = {V[i]->coords[0]-V[j]->coords[0], V[i]->coords[1]-V[j]->coords[1]};
            // cout<<"\nid0 : "<<id0<<"\nid1 : "<<id1<<endl;

            // vector IJ must point from the smaller vertex ID to the greater
            if ( id0 > id1 ){
                for (int jj=0; jj<2; jj++){
                    IJ[jj] = -IJ[jj];
                }
            }
            // vector: from element center to edge middle point, used as a reference vector
            double v[2] = {EdgeCG[0]-ElemCG[0],EdgeCG[1]-ElemCG[1]};

            // Cij is orthogonal to v
            for (int jj=0; jj<2; jj++){
                Cij[jj] = .0;
            }


            // Cij must point as if I inside the CV and J outside
            double innerprod = v[1]*IJ[0] + (-v[0])*IJ[1];
            if ( innerprod <= .0 ){
                for (int jj=0; jj<2; jj++){
                    v[jj] = -v[jj];
                }
            }

            // associate Cij coefficient to edge
            Cij[0] += v[1];
            Cij[1] += -v[0];

            // apos calculado o vetor Cij, este é "guardado" na aresta da malha
            p[currentedge]=new double[3];
            p[currentedge][0]=p[currentedge][1]=0;
            p[currentedge][0]=Cij[0];
            p[currentedge][1]=Cij[1];
            QuadEdgeAdress[i]->Cij[dom]=p[currentedge];
            CijSum[0]+=QuadEdgeAdress[i]->Cij[dom][0];CijSum[1]+=QuadEdgeAdress[i]->Cij[dom][1];
            currentedge++;


           //cout << "Cij = \t" << QuadEdgeAdress[i]->Cij[QuadEdgeAdress[i]->physical][0] << "\t" << QuadEdgeAdress[i]->Cij[QuadEdgeAdress[i]->physical][1] << "\t" << QuadEdgeAdress[i]->Cij[QuadEdgeAdress[i]->physical][2];
           // cout<<" :: domain flag : "<<dom<<endl;
            //Check if Net Cij is Zero

            /// Calculate Dij coefficient only for boundary edges.

            if (QuadEdgeAdress[i]->physical!=dom){
            numBE++;//cout<<"\nboundary edge found :\t"<<id0<<"-"<<id1<<"\n";
            QuadEdgeAdress[i]->Dij=new double[3];
            QuadEdgeAdress[i]->Dij[0]=QuadEdgeAdress[i]->Dij[1]=QuadEdgeAdress[i]->Dij[2]=0;
            // Dij vector is orthogonal to edge (it's unknown Dij orientation)
            Dij[0] = -(V[i]->coords[1]-V[j]->coords[1])/2.0;
            Dij[1] =  (V[i]->coords[0]-V[j]->coords[0])/2.0;

             // make Dij points to outside domain. First, take the face that uses edge and its flag


            // vector: from element center to edge middle point, used as a reference vector
            v[0] = EdgeCG[0]-ElemCG[0];
            v[1] = EdgeCG[1]-ElemCG[1];

            ///if(Flag do Quad Analisado no loop > Flag do Quad que contém a mesma aresta){
            //Dij points from the lowest flag to the highest flag domain

            ///}

            // Dij must point to outside element
            double innerprod = Dij[0]*v[0] + Dij[1]*v[1];
            if (  innerprod <= .0 ){
                for (j=0; j<2; j++){
                    Dij[j] = -Dij[j];
                }
            }

                QuadEdgeAdress[i]->Dij[0]=Dij[0];
                QuadEdgeAdress[i]->Dij[1]=Dij[1];
                cout << "Dij = \t" << QuadEdgeAdress[i]->Dij[0] << "\t" << QuadEdgeAdress[i]->Dij[1] << "\t" << QuadEdgeAdress[i]->Dij[2]<<endl;
                DijSum[0]+=QuadEdgeAdress[i]->Dij[0]; DijSum[1]+=QuadEdgeAdress[i]->Dij[1];

            }

        }//end of edge loops

	}//end of Quad loops
	//verify_coeff_calculation() ;
    cout<<"\nCij sum -> ("<<CijSum[0]<<","<<CijSum[1]<<","<<CijSum[2]<<")\n";
    cout<<"Dij sum -> ("<<DijSum[0]<<","<<DijSum[1]<<","<<DijSum[2]<<")\n";










}
}
