/*
 * GeomData.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{
	void GeomData::initializeElementMatrix(int n){
		numElem = n;
		geoElementMat.allocateMemory(numElem,18);
	}

	void GeomData::setElementMatrices(int row, const double* Aij, const double* Ajk, const double* Aik){
		for (int i=0; i<6; i++){
			geoElementMat.setValue(row,i,Aij[i]);
			geoElementMat.setValue(row,i+6,Ajk[i]);
			geoElementMat.setValue(row,i+12,Aik[i]);
		}
	}

	void GeomData::getElementMatrices(int row, double* Aij, double* Ajk, double* Aik){
		for (int i=0; i<6; i++){
			Aij[i] = geoElementMat.getValue(row,i);
			Ajk[i] = geoElementMat.getValue(row,i+6);
			Aik[i] = geoElementMat.getValue(row,i+12);
		}
	}

	void GeomData::getVolume_MEBFV(int id, double &Vi){
		Vi = volume_MEBFV[id-1];
	}
}
