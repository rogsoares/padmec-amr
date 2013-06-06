/*
 * MeshAdaptaion__Remeshing.h
 *
 *  Created on: 24/01/2013
 *      Author: rogsoares
 */

#ifndef MESHADAPTAION__REMESHING_H_
#define MESHADAPTAION__REMESHING_H_

#include "AMR.h"

class AdaptiveRemeshing : public AMR{
public:

	AdaptiveRemeshing(){}
	~AdaptiveRemeshing(){}

	//void rodar(ErrorAnalysis*,pMesh);
	void rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh);
};


#endif /* MESHADAPTAION__REMESHING_H_ */
