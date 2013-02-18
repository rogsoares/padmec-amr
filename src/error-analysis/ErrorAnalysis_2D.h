/*
 * ErrorAnalysis-2D.h
 *
 *  Created on: 24/04/2012
 *      Author: rogsoares
 */

#ifndef ERRORANALYSIS_2D_H_
#define ERRORANALYSIS_2D_H_

#include "ErrorAnalysis.h"


/*! \class ErrorAnalysis_2D
 *  \brief Used to perform error analysis in 2-D unstructured (triangles) meshes. General analysis functions are implemented in base class.
 */
class ErrorAnalysis_2D : public ErrorAnalysis{

public:

	ErrorAnalysis_2D(){}
	~ErrorAnalysis_2D(){}

private:
	/*! \brief: Realize
	 * \param theMesh FMDB mesh object
	 * \param GetPFuncGrad function pointer to get nodal gradient
	 */
	void calculate_ElementsError(pMesh theMesh, SimulatorParameters* pSimPar, GetPFuncGrad);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_SmoothedGradientNorm(pMesh, SimulatorParameters *pSimPar, GetPFuncGrad);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_SmoothedGradientNorm_excludingSingularities(pMesh, SimulatorParameters*, GetPFuncGrad);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_CharacteristicDimensionLength(pMesh);

	/*! \brief: Realize the data transfer from the fine to the coarse mesh.
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	void calculate_DegreeOfRefinement(pMesh, int, int, bool);
};



#endif /* ERRORANALYSIS_2D_H_ */
