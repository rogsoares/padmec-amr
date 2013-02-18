/*
 * interpolateBetweenMeshes.cpp
 *
 *  Created on: 11/02/2013
 *      Author: rogsoares
 */


#include "interpolation.h"

void calculate_GeometricCoefficients(InterpolationDataStruct* pIData, int dim){
	pEntity vertice_1,vertice_2,vertice_3,vertice_4; //element vertices
	double ptn1[3],ptn2[3],ptn3[3],ptn4[3];

	if(dim==2){
		VIter vit = M_vertexIter(pIData->m2);
		while ( pEntity vertex = VIter_next(vit) ){

			// get coordinate of a point from the mesh to be interpolated
			double xyz[3];
			V_coord(vertex,xyz);

			// search entity in which the point above are
			pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);
			if (!element){
				cout << "Entity not found! Exiting...\n";
				exit(1);
			}

			//get coordinates
			vertice_1 = element->get(0,0);
			vertice_2 = element->get(0,1);
			vertice_3 = element->get(0,2);
			V_coord(vertice_1,ptn1);
			V_coord(vertice_2,ptn2);
			V_coord(vertice_3,ptn3);

			//calculate areas
			double a1 = F_area(xyz,ptn2,ptn3);
			double a2 = F_area(xyz,ptn1,ptn3);
			double a3 = F_area(xyz,ptn1,ptn2);
			double a = a1 + a2 + a3;
			a1 = a1/a;
			a2 = a2/a;
			a3 = a3/a;

			//calculate coefficients
			EN_attachDataDbl(element, MD_lookupMeshDataId("GC_a1"),a1);
			EN_attachDataDbl(element, MD_lookupMeshDataId("GC_a2"),a2);
			EN_attachDataDbl(element, MD_lookupMeshDataId("GC_a3"),a3);
		}
		VIter_delete(vit);
	}
//	else {
//		VIter vit = M_vertexIter(pIData->m1);
//		while ( pEntity vertex = VIter_next(vit) ){
//
//			// get coordinate of a point from the mesh to be interpolated
//			double xyz[3];
//			V_coord(vertex,xyz);
//			int row = EN_id(vertex);
//			// search entity in which the point above are
//
//			pEntity element = (pEntity)Octree_Search(xyz,theOctree);
//			if (!element){
//				cout << "Entity not found! Exiting...\n";
//				exit(1);
//			}
//
//			//get coordinates
//			vertice_1= element->get(0,0); V_coord(vertice_1,ptn1);
//			vertice_2= element->get(0,1); V_coord(vertice_2,ptn2);
//			vertice_3= element->get(0,2); V_coord(vertice_3,ptn3);
//			vertice_4= element->get(0,3); V_coord(vertice_4,ptn4);
//
//			//calculate areas
//			double v1 = R_Volume(xyz,ptn2,ptn3,ptn4);
//			double v2 = R_Volume(xyz,ptn1,ptn3,ptn4);
//			double v3 = R_Volume(xyz,ptn1,ptn2,ptn4);
//			double v4 = R_Volume(xyz,ptn1,ptn2,ptn3);
//			double v = v1 + v2 + v3 + v4;
//			v1 = v1/v; v2= v2/v; v3= v3/v; v4= v4/v;
//
//			//calculate coefficients
//			setGeomCoeff(row,0,v1);
//			setGeomCoeff(row,1,v2);
//			setGeomCoeff(row,2,v3);
//			setGeomCoeff(row,3,v4);
//		}
//
//		VIter_delete(vit);
//	}
}


void calculate_LinearInterpolation(InterpolationDataStruct* pIData, int dim){
	string geocoeffstr[4] = {"GC_a1","GC_a2","GC_a3","GC_a4"};
	double xyz[3], scalar, val;

	//Loop over vertices (from the mesh to where the data is being transfered)
	VIter vit = M_vertexIter(pIData->m1);
	while ( pEntity vertex = VIter_next(vit) ){
		//get the node coordinate

		V_coord(vertex,xyz);

		//search entity in which the node is
		pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);

		// Interpolate all fields assigned to mesh node
		int size = dim + 1;
		double geo_coeff = .0;
		for (int field=0; field<pIData->numFields; field++){
			// get element's nodes data per field
			for(int i=0; i<size; i++){
				scalar = pIData->pGetDblFunctions[field]( element->get(0,i) );
				EN_getDataDbl(element, MD_lookupMeshDataId(geocoeffstr[i].c_str()), &geo_coeff);
				val += scalar*geo_coeff;
			}
			pIData->pSetDblFunctions[field](vertex,val);
		}
	}
	VIter_delete(vit);
}

//double calculate_QuadraticInterpolation(pMesh theMesh, int dim){
//	double NodeValue, GeomCoeff, residual, coord[3], linear, quadratic;
//	pEntity v;
//	double xyz[3];
//	//Loop over vertices (from the mesh to where the data is being transfered)
//	VIter vit = M_vertexIter(theMesh);
//	while ( pEntity vertex = VIter_next(vit) ){
//		V_coord(vertex,xyz);
//		int row = EN_id(vertex);
//		//search entity in which the node is
//		pEntity element = (pEntity)Octree_Search(xyz,theOctree);
//		if (!element){
//			cout << "WARNING:  Entity not found!\n";
//		}
//
//		//calculate residual
//		residual = 0;
//		for( int i = 0; i < dim+1; i++){
//			v = element->get(0,i);
//			int id= EN_id(v);
//			V_coord(v,coord);
//			GeomCoeff= getGeomCoeff(row,i);
//			double dx = getGrad(id,1);
//			double dy = getGrad(id,2);
//			double dx2 = getGrad(id,3);
//			double dy2 = getGrad(id,4);
//			double dxdy = getGrad(id,5);
//
//			double Res = (xyz[0]-coord[0])*dx + (xyz[1]-coord[1])*dy;
//			Res += (xyz[0]-coord[0])*(xyz[0]-coord[0])*dx2/2;
//			Res += (xyz[1]-coord[1])*(xyz[1]-coord[1])*dy2/2;
//			Res += (xyz[0]-coord[0])*(xyz[1]-coord[1])*dxdy;
//			residual += GeomCoeff*Res;
//		}
//
//		//calculate quadratic interpolation
//		linear = getInterpolatedValues(row);
//		quadratic = linear + residual;
//		setInterpolatedValues(row,quadratic);
//	}
//	VIter_delete(vit);
//}
//
//
//void calculate_DerivativesError(pMesh theMesh){
//	double summ_dx = 0; double summ2_dx = 0; double MaxError_dx = 0;
//	double summ_dy = 0; double summ2_dy = 0; double MaxError_dy = 0;
//	double summ_dx2 = 0; double summ2_dx2 = 0; double MaxError_dx2 = 0;
//	double summ_dy2 = 0; double summ2_dy2 = 0; double MaxError_dy2 = 0;
//	double summ_dxdy = 0; double summ2_dxdy = 0; double MaxError_dxdy = 0;
//
//	//loop on the mesh vertices
//	VIter vit = M_vertexIter(theMesh);
//	while (pEntity vertex = VIter_next(vit)){
//
//		int id = EN_id(vertex);
//		//get coordinates
//		double xyz[3] = {0,0,0};
//		V_coord(vertex, xyz);
//
//		//ATTENTION: THE CORRECT GRADIENTS VALUES MUST BE PUT HERE BEFORE USE THIS FUNCTION
//		double dx = 0;
//		double dy = 2*xyz[1];
//		double dx2 = 0;
//		double dy2 = 2;
//		double dxdy = 0;
//
//		//calculate the global error (2-norm) related to the gradients
//		double error_dx = getGrad(id,1) - dx;
//		double error_dy = getGrad(id,2) - dy;
//		double error_dx2 = getGrad(id,3) - dx2;
//		double error_dy2 = getGrad(id,4) - dy2;
//		double error_dxdy = getGrad(id,5) - dxdy;
//
//		summ_dx = summ_dx + pow(error_dx,2);
//		summ2_dx = summ2_dx + pow(dx,2);
//
//		summ_dy = summ_dy + pow(error_dy,2);
//		summ2_dy = summ2_dy + pow(dy,2);
//
//		summ_dx2 = summ_dx2 + pow(error_dx2,2);
//		summ2_dx2 = summ2_dx2 + pow(dx2,2);
//
//		summ_dy2 = summ_dy2 + pow(error_dy2,2);
//		summ2_dy2 = summ2_dy2 + pow(dy2,2);
//
//		summ_dxdy = summ_dxdy + pow(error_dxdy,2);
//		summ2_dxdy = summ2_dxdy + pow(dxdy,2);
//
//		//calculate the maximmun node error
//		error_dx = fabs(error_dx);
//		if ( error_dx > MaxError_dx ) MaxError_dx = error_dx;
//
//		error_dy = fabs(error_dy);
//		if ( error_dy > MaxError_dy ) MaxError_dy = error_dy;
//
//		error_dx2 = fabs(error_dx2);
//		if ( error_dx2 > MaxError_dx2 ) MaxError_dx2 = error_dx2;
//
//		error_dy2 = fabs(error_dy2);
//		if ( error_dy2 > MaxError_dy2 ) MaxError_dy2 = error_dy2;
//
//		error_dxdy = fabs(error_dxdy);
//		if ( error_dxdy > MaxError_dxdy ) MaxError_dxdy = error_dxdy;
//
//	}
//}
