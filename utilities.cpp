/*
 * Utilities.cpp
 *
 *  Created on: 01/01/2015
 *      Author: rogerio
 */

#include "utilities.h"

void middlePoint(const double* p1, const double* p2, double* mp){
	mp[0] = 0.5*(p1[0]+p2[0]);
	mp[1] = 0.5*(p1[1]+p2[1]);
	mp[2] = 0.5*(p1[2]+p2[2]);
}

void middlePoint(const double* p1, const double* p2, const double* p3, const double* p4, double* mp){
	for (int i=0;i<3;i++){
		mp[i] = 0.25*(p1[i]+p2[i]+p3[i]+p4[i]);
	}
}

double R_Volume(const double* p1, const double* p2, const double* p3, const double* p4){
	double a[3] = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] } ;
	double b[3] = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] } ;
	double c[3] = { p4[0]-p1[0], p4[1]-p1[1], p4[2]-p1[2] } ;
	double normal[3];
	computeCrossProduct( a, b, normal ) ;
	return computeDotProduct( normal, c )/6.0;
}

void computeCrossProduct(const double* a, const double* b, double* n){
	n[0] = a[1]*b[2] - b[1]*a[2] ;
	n[1] = a[2]*b[0] - a[0]*b[2] ;
	n[2] = a[0]*b[1] - a[1]*b[0] ;
}

double computeDotProduct(const double* a, const double* b){
	double sum = 0.0 ;
	for( int i = 0 ; i < 3 ; ++i ){
		sum += a[i]*b[i];
	}
	return sum ;
}

