/*******************************************************************************************
 * This is the pre-processor code for the EBFV1 (Edge based finite volume 1) formulation.
 * It's based on the FVM discretization presented in Carvalho (2005) thesis. This pre-pro-
 * cessor leads only with triagular elements for serial and paralllel simulations.
 *
 * Developed by: Rogerio S. da Silva     2005-2010-2012
 ********************************************************************************************/

#include "EBFV1__pre-processors.h"


// calculates: Cij, Dij and nodal volume
// ---------------------------------------------------------------------
int EBFV1_preprocessor_2D(pMesh theMesh, void *pData, int &ndom){
	cout << "calculateCoefficients2D_old\n";

	GeomData *pGCData = (GeomData*)pData;
	pGCData->setMeshDim(theMesh->getDim());

	theMesh->modifyState(2,1);	// create edge data structure
	theMesh->modifyState(0,2);	// create adjacency around nodes
	theMesh->modifyState(1,2);	// create faces around edge

	cout << "Mesh dimension: 2D\n";
	int i,j,flag,dom;
	double edgeCenter[3], I[3], J[3];

	std::vector<double> Cij(3), Dij(3);
	pEntity edge, face;

	/// initialize coefficients
	initializeCoefficients(theMesh,pGCData);

	double H = pGCData->getReservoirHeight();
	std::set<int> setOfDomain;

	// for each face, calculate all data needed and attach them to the its edges
	FIter fit = M_faceIter(theMesh);
	while ( (face = FIter_next(fit)) ){
		// look up only for leave elements (without children)
		if ( !theMesh->getRefinementDepth(face) ){
			dom = EN_getFlag(face);

			setOfDomain.insert(dom);		// store every new domain
			double faceCenter[3] = {.0, .0, .0};
			getFCenter(face,faceCenter);

			// loop over all three face's edges to calculate Cij and Dij
			for (i=0; i<3; i++){
				edge = F_edge(face, i);
				edgeCenter[0] = edgeCenter[1] = .0;
				E_center(edge,edgeCenter);				// edge centroid
				V_coord(edge->get(0,0),I);
				V_coord(edge->get(0,1),J);

				int id0 = EN_id(edge->get(0,0));
				int id1 = EN_id(edge->get(0,1));

				// edge vector, used as a reference vector
				double IJ[2] = {J[0]-I[0], J[1]-I[1]};

				// vector IJ must point from the smaller vertex ID to the greater
				if ( id0 > id1 ){
					for (j=0; j<2; j++){
						IJ[j] = -IJ[j];
					}
				}

				// vector: from element center to edge middle point, used as a reference vector
				double v[2] = {edgeCenter[0]-faceCenter[0],edgeCenter[1]-faceCenter[1]};

				// Cij is orthogonal to v
				for (j=0; j<2; j++){
					Cij[j] = .0;
				}

				// Cij must point as if I inside the CV and J outside
				double innerprod = v[1]*IJ[0] + (-v[0])*IJ[1];
				if ( innerprod <= .0 ){
					for (j=0; j<2; j++){
						v[j] = -v[j];
					}
				}

				// associate Cij coefficient to edge
				pGCData->getCij(edge,dom,Cij);
				Cij[0] += v[1];
				Cij[1] += -v[0];

				for (j=0; j<2; j++){
					Cij[j] = H*Cij[j];
				}
				pGCData->setCij(edge,dom,Cij);
			}

			// calculate volume of control volume and associate it to elements nodes
			const double porosity = .0;
			double A = F_area(face)/3.0;	// element area
			for (j=0; j<3; j++){
				double v1 = pGCData->getVolume(face->get(0,j),dom);
				double v2 = pGCData->getWeightedVolume(face->get(0,j));
				v1 += A;
				v2 += A*porosity;
				pGCData->setVolume(face->get(0,j),dom,v1);
				pGCData->setWeightedVolume(face->get(0,j),v2);
			}
		}
	}
	FIter_delete(fit);

	// Calculate Dij coefficient only for boundary edges.
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		if ( !theMesh->getRefinementDepth(edge) ){
			flag = EN_getFlag(edge);
			std::set<int>::iterator iter = setOfDomain.find( flag );
			if ( iter == setOfDomain.end() ){
				double I[3] = {.0,.0,.0}, J[3] = {.0,.0,.0};
				V_coord(edge->get(0,0),I);
				V_coord(edge->get(0,1),J);

				// Dij vector is orthogonal to edge (it's unknown Dij orientation)
				Dij[0] = -(J[1]-I[1])/2.0;
				Dij[1] =  (J[0]-I[0])/2.0;

				// make Dij points to outside domain. First, take face that uses edge and its flag
				int domains[2] = {0,0};
				for (i=0; i<E_numFaces(edge); i++){
					face = E_face(edge,i);
					if (!face){
						throw Exception(__LINE__,__FILE__,"Null face!\n");
					}
					// TODO: change this for heterogeneous media!!!!
					domains[i] = 3300; //EN_getFlag(face);
				}
				// that the reference face to make Dij points to outside
				face = E_face(edge, 0);

				if (!face){
					throw Exception(__LINE__,__FILE__,"Null face!\n");
				}

				// now, get face's center and edge's center
				double faceCenter[3] = {.0, .0, .0}, edgeCenter[3] = {.0, .0, .0};
				getFCenter(face,faceCenter);	// element centroid
				E_center(edge,edgeCenter);		// edge centroid

				// vector: from element center to edge middle point, used as a reference vector
				double v[2] = {edgeCenter[0]-faceCenter[0],edgeCenter[1]-faceCenter[1]};

				// Dij must point to outside element
				double innerprod = Dij[0]*v[0] + Dij[1]*v[1];
				if (  innerprod <= .0 ){
					for (j=0; j<2; j++){
						Dij[j] = -Dij[j];
					}
				}

				// associate to edge domains flags to wjich it belongs
				EN_attachDataInt(edge,MD_lookupMeshDataId("dom1"),domains[0]);
				EN_attachDataInt(edge,MD_lookupMeshDataId("dom1"),domains[1]);

				// associate to edge Dij vector
				// REMEMBER THAT: Dij points domains[0] -> domains[1] ALWAYS!!!
				for (j=0; j<2; j++){
					Dij[j] *= H;
				}
				pGCData->setDij(edge,domains[0],domains[1],Dij);
			}
		}
	}
	EIter_delete(eit);

	/*
	 * Parallel step:
	 *
	 * unifyCijAmongProcessors: for any edge IJ on partition boundary, all ranks
	 * 		that share it calculate vector CIJ locally from local tetrahedrals.
	 * 		What is made here consist to sum (unify) all local Cij and retrieve
	 * 		the result to all processes that share IJ. This procedure is equivalent
	 * 		to the MPI function "MPI_Allreduce" where all ranks send a specific
	 * 		value and receive the sum.
	 *
	 * unifyVolumesAmongProcessors: the same as unifyCijAmongProcessors but with
	 * 		scalars.
	 *
	 * NOTE: All partitions must be aware of the total number of domains.
	 *
	 * */
//	AllgatherDomains(setOfDomain);
//	unifyCijAmongProcessors(theMesh,setOfDomain,pGCData);
//	unifyVolumesAmongProcessors(theMesh,setOfDomain,"volume",pGCData);


	// TODO: mudar argumentos de entrada para validacao
	//validate_EBFV1(pPPP);

	calculateEdgeLength(theMesh,pGCData);
	calculateCijNorm(theMesh,pGCData,setOfDomain);


	/*
	 * For 2-D domains, multiply volume by reservoir height (H) for 2D/3D simulations (physics occurs only on 2-D but
	 * reservoir volume is considered)
	 * =======================================================================================
	 */
	double vt = .0;
	double vol,wvol;
	std::set<int>::iterator iter = setOfDomain.begin();
	for(;iter!=setOfDomain.end();iter++){
		VIter vit = M_vertexIter(theMesh);
		while (pEntity node = VIter_next(vit)){
			// do not divide vol by number of remote copies!
			vol = H*pGCData->getVolume(node,*iter);
			vt += vol;
			wvol = H*0.2*vol;
			if (vol == .0){
				throw Exception(__LINE__,__FILE__,"Volume is NULL!\n");
			}
			pGCData->setWeightedVolume(node,wvol);
		}
		VIter_delete(vit);
	}
	pGCData->setTotalReservoirVolume(vt);
	return 0;
}


void initializeCoefficients(pMesh theMesh, GeomData *pGCData){
	VIter vit = M_vertexIter(theMesh);
	while (pEntity v = VIter_next(vit)){
		pGCData->setWeightedVolume(v,0.0);
		pGCData->setVolume(v,3300,.0);
	}
	VIter_delete(vit);

	std::vector<double> vec(3,.0);
	EIter eit = M_edgeIter(theMesh);
	while ( pEdge edge = EIter_next(eit) ){
		pGCData->setCij(edge,3300,vec);
		if (EN_getFlag(edge)==2000){
			pGCData->setDij(edge,2000,0,vec);
		}
	}
	EIter_delete(eit);
}
