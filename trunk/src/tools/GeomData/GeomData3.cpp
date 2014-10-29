/*
 * GeomData.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: rogerio
 */

#include "GeomData.h"

namespace PRS{
	void GeomData::initializeElementMatrix(int n){
		this->numElem = n;
		this->geoElementMat.allocateMemory(numElem,18);
		idxElemMat.allocateMemory(numElem,3);
	}

	void GeomData::setElementMatrices(int row, const double* Aij, const double* Ajk, const double* Aik){
		for (int i=0; i<6; i++){
			this->geoElementMat.setValue(row,i,Aij[i]);
			this->geoElementMat.setValue(row,i+6,Ajk[i]);
			this->geoElementMat.setValue(row,i+12,Aik[i]);
		}
	}

	void GeomData::getElementMatrices(int row, double* Aij, double* Ajk, double* Aik){
		for (int i=0; i<6; i++){
			Aij[i] = this->geoElementMat.getValue(row,i);
			Ajk[i] = this->geoElementMat.getValue(row,i+6);
			Aik[i] = this->geoElementMat.getValue(row,i+12);
		}
	}

	void GeomData::getIJK_IDs(int row, int& id0,int& id1,int& id2){
		id0 = idxElemMat.getValue(row,0);
		id1 = idxElemMat.getValue(row,1);
		id2 = idxElemMat.getValue(row,2);
	}

	void GeomData::setIJK_IDs(int row, int id0, int id1, int id2){
		idxElemMat.setValue(row,0,id0);
		idxElemMat.setValue(row,1,id1);
		idxElemMat.setValue(row,2,id1);
	}

	void GeomData::getVolume_MEBFV(int id, double &Vi){
		Vi = volume_MEBFV[id-1];
	}
}
