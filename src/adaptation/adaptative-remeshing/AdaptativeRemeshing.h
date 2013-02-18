/*
 * MeshAdaptaion__Remeshing.h
 *
 *  Created on: 24/01/2013
 *      Author: rogsoares
 */

#ifndef MESHADAPTAION__REMESHING_H_
#define MESHADAPTAION__REMESHING_H_

#include "AMR.h"

class AdaptativeRemeshing : public AMR{
public:

	AdaptativeRemeshing(){}
	~AdaptativeRemeshing(){}

	void run(ErrorAnalysis*,pMesh);
};


#endif /* MESHADAPTAION__REMESHING_H_ */
