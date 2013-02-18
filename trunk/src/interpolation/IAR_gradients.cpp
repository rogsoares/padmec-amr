///*
// * IAR_gradients.cpp
// *
// * IAR_gradients (Interpolation for Adaptative Remeshing)
// *
// *  Created on: 15/02/2013
// *      Author: rogsoares
// */
//
//
//#include "KSP_solver.h"
//
//double calculate_Gradients(InterpolationData* pIntpData){
//	int i;
//	int NumOfNodes= M_numVertices(theMesh);
//	double *dx, *dy;
//
//	//PETSc variables
//	Vec Fx, Fy, gradx, grady;
//	Mat M;
//	PetscErrorCode ierr;
//
//	//create vectors Fx, Fy, gradx and grady
//	ierr = VecCreate(PETSC_COMM_WORLD,&Fx);CHKERRQ(ierr);
//	ierr = VecSetSizes(Fx,PETSC_DECIDE,NumOfNodes);CHKERRQ(ierr);
//	ierr = VecSetFromOptions(Fx);CHKERRQ(ierr);
//	ierr = VecDuplicate(Fx,&Fy);CHKERRQ(ierr);
//	ierr = VecDuplicate(Fx,&gradx);CHKERRQ(ierr);
//	ierr = VecDuplicate(Fx,&grady);CHKERRQ(ierr);
//
//	//create an sparse matrix
//	ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,NumOfNodes,NumOfNodes,80,PETSC_NULL,80,PETSC_NULL,&M); CHKERRQ(ierr);
//
//	//Assembly
//	Assembly_Mat_Vec(theMesh,M,Fx,Fy);
//
//	//ierr = VecView(Fx,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
//
//	//Solve the linear system
//	KSP_solver(M, M, Fx, gradx, KSPCG, PCJACOBI);
//	KSP_solver(M, M, Fy, grady, KSPCG, PCJACOBI);
//
//	//store first order derivatives
//	store_FirstOrderDerivatives(theMesh,gradx,grady);
//
//	//Second order derivatives
//	calculate_SecondOrderDerivatives(theMesh);
//
//	//destroy vectors, matrix and KSP's
//	ierr = VecDestroy(Fx);CHKERRQ(ierr);
//	ierr = VecDestroy(Fy);CHKERRQ(ierr);
//	ierr = VecDestroy(gradx);CHKERRQ(ierr);
//	ierr = VecDestroy(grady);CHKERRQ(ierr);
//	ierr = MatDestroy(M);CHKERRQ(ierr);
//}
//
//
//double Assembly_Mat_Vec(pMesh theMesh,Mat M, Vec Fx, Vec Fy) {
//	pEntity vertice_1,vertice_2,vertice_3,vertice_4;
//	double coord1[3],coord2[3],coord3[3], value, value1, value2, value3, Fx_values[3],Fy_values[3], M_values[9];
//	PetscErrorCode ierr;
//
//	// to reduce de number of flop
//	const double byTwelve = 1./12.;
//	const double bySix = 1./6.;
//
//	//Assembly
//	FIter fit = M_faceIter(theMesh);
//	while (pEntity entity = FIter_next(fit)){
//
//		//get the coordinates and the id's from the element's vertices
//		vertice_1 = entity->get(0,0);
//		vertice_2 = entity->get(0,1);
//		vertice_3 = entity->get(0,2);
//		V_coord(vertice_1,coord1);
//		V_coord(vertice_2,coord2);
//		V_coord(vertice_3,coord3);
//		int position[3] = {EN_id(vertice_1)-1, EN_id(vertice_2)-1, EN_id(vertice_3)-1};
//
//		//calculate area
//		double At = F_area(coord1,coord2,coord3);
//
//		//get node values
//		value1 = getNodeValue(EN_id(vertice_1));
//		value2 = getNodeValue(EN_id(vertice_2));
//		value3 = getNodeValue(EN_id(vertice_3));
//
//		//calculate the contributions
//		for(int i=0;i<3;i++){
//			Fx_values[i] = (coord2[1]-coord3[1])*value1*bySix + (coord3[1]-coord1[1])*value2*bySix + (coord1[1]-coord2[1])*value3*bySix;
//			Fy_values[i] = (coord3[0]-coord2[0])*value1*bySix + (coord1[0]-coord3[0])*value2*bySix + (coord2[0]-coord1[0])*value3*bySix;
//			for(int j=0;j<3;j++){
//				if (i==j){
//					M_values[3*i + j] = At*bySix; //At/6;
//				}
//				else{
//					M_values[3*i + j] = At*byTwelve; //At/12;
//				}
//			}
//		}
//		VecSetValues(Fx,3,position,Fx_values, ADD_VALUES);
//		VecSetValues(Fy,3,position,Fy_values, ADD_VALUES);
//		MatSetValues(M,3,position,3,position,M_values, ADD_VALUES);
//	}
//	FIter_delete(fit);
//	ierr = MatAssemblyBegin(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//	ierr = MatAssemblyEnd(M,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//	ierr = VecAssemblyBegin(Fx);CHKERRQ(ierr);
//	ierr = VecAssemblyEnd(Fx);CHKERRQ(ierr);
//	ierr = VecAssemblyBegin(Fy);CHKERRQ(ierr);
//	ierr = VecAssemblyEnd(Fy);CHKERRQ(ierr);
//}
//
//double store_FirstOrderDerivatives(pMesh theMesh,Vec gradx, Vec grady){
//	PetscErrorCode ierr;
//	double *dx, *dy;
//	ierr = VecGetArray(gradx,&dx);CHKERRQ(ierr);
//	ierr = VecGetArray(grady,&dy);CHKERRQ(ierr);
//
//	//loop on vertices to calculate the secord order derivatives
//	VIter vit = M_vertexIter(theMesh);
//	while ( pEntity vertex = VIter_next(vit) ){
//		double coord[3];
//		int id = EN_id(vertex);
//		V_coord(vertex,coord);
//
//		//store the first order derivatives
//		setGrad(id,1,dx[id - 1]);
//		setGrad(id,2,dy[id - 1]);
//	}
//	VIter_delete(vit);
//	ierr = VecRestoreArray(gradx,&dx);CHKERRQ(ierr);
//	ierr = VecRestoreArray(grady,&dy);CHKERRQ(ierr);
//}
//
//double calculate_SecondOrderDerivatives(pMesh theMesh) {
//	//loop on vertices to calculate the secord order derivatives
//	VIter vit = M_vertexIter(theMesh);
//	while ( pEntity vertex = VIter_next(vit) ){
//		//find out number of elements made with vertex
//		int NumOfNeighbours = V_numFaces(vertex);
//		if (!NumOfNeighbours){
//			throw Exception(__LINE__,__FILE__,"Any neighbors found!\n");
//		}
//		pEntity vertice1, vertice2, vertice3;
//		int id = EN_id(vertex);
//		double coord[3];
//		V_coord(vertex,coord);
//
//		//calculate second order derivatives on all elements around a certain vertex and compute a media
//		double dx2 = 0;
//		double dy2 = 0;
//		double dxdy = 0;
//
//		//find the elements that have the vertex and calculate their contribution
//		for(int i=0; i<NumOfNeighbours; i++){
//			double ptn1[3],ptn2[3],ptn3[3];
//
//			//get the coordinates and id's
//			pEntity face = vertex->get(2,i);
//			if (!face){
//				throw Exception(__LINE__,__FILE__,"Null face.\n");
//			}
//
//			//get coordinates and id from vertice 1
//			vertice1 = face->get(0,0);
//			V_coord(vertice1,ptn1);
//			int id1 = EN_id(vertice1);
//			double Dx1 = getGrad(id1,1);
//			double Dy1 = getGrad(id1,2);
//
//			//get coordinates and id from vertice 2
//			vertice2 = face->get(0,1);
//			V_coord(vertice2,ptn2);
//			int id2 = EN_id(vertice2);
//			double Dx2 = getGrad(id2,1);
//			double Dy2 = getGrad(id2,2);
//
//			//get coordinates and id from vertice 3
//			vertice3 = face->get(0,2);
//			V_coord(vertice3,ptn3);
//			int id3 = EN_id(vertice3);
//			double Dx3 = getGrad(id3,1);
//			double Dy3 = getGrad(id3,2);
//
//			//calculate area
//			double Area = F_area(ptn1,ptn2,ptn3);
//
//			//calculate the contributions
//			dx2 = dx2 + 0.5*(ptn2[1]-ptn3[1])*Dx1/Area + 0.5*(ptn3[1]-ptn1[1])*Dx2/Area + 0.5*(ptn1[1]-ptn2[1])*Dx3/Area;
//			dy2 = dy2 + 0.5*(ptn3[0]-ptn2[0])*Dy1/Area + 0.5*(ptn1[0]-ptn3[0])*Dy2/Area + 0.5*(ptn2[0]-ptn1[0])*Dy3/Area;
//			dxdy = dxdy + 0.5*(ptn2[1]-ptn3[1])*Dy1/Area + 0.5*(ptn3[1]-ptn1[1])*Dy2/Area + 0.5*(ptn1[1]-ptn2[1])*Dy3/Area;
//		}
//
//		//calculate the second order derivatives as a media
//		dx2 = dx2/NumOfNeighbours;
//		setGrad(id,3,dx2);
//		dy2 = dy2/NumOfNeighbours;
//		setGrad(id,4,dy2);
//		dxdy = dxdy/NumOfNeighbours;
//		setGrad(id,5,dxdy);
//	}
//	VIter_delete(vit);
//}
