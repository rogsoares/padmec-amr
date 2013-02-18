/*
 * ErrosAnalysis.h
 *
 *  Created on: 20/02/2012
 *      Author: rogsoares
 */

#ifndef ERRORANALYSIS_H_
#define ERRORANALYSIS_H_

/*
 * Get estimated error from mesh elements.
 */

#include "SimulatorParameters.h"
using namespace PRS;

class ErrorAnalysis;

/*! \brief: Main function in the error analysis process.
 * \param ErrorAnalysis* polymorphic object pointer to ErrorAnalysis_2D or ErrorAnalysis_3D class
 * \param pMesh FMDB mesh object
 * \param double error tolerance for all mesh elements
 * \param double error tolerance for all mesh elements excluding singularities
 * \param GetPFuncGrad* pointer to an array of function pointers to get nodal gradients (0 - p_grad, 1 - Sw_grad)
 * \param int size of GetPFuncGrad array
 */

bool calculate_ErrorAnalysis(ErrorAnalysis*, pMesh, SimulatorParameters*, double, double, GetPFuncGrad*, int);

/*! \class ErrorAnalysis
 *  \brief This class does some basic operations to perform error analysis and serves as base class for further implementation for 2-D and 3-D
 *  unstructured meshes. A polymorphic object must be allocated for 2-D (ErrorAnalysis_2D) or 3-D (ErrorAnalysis_3D) to perform the desired
 *  calculations.
 */
class ErrorAnalysis{

public:

	ErrorAnalysis(){
		elem_id = MD_lookupMeshDataId( "elem_error" );
		levelRef_id = MD_lookupMeshDataId( "elem_ref" );
		sing_id = MD_lookupMeshDataId( "elem_sing" );
		cdl_id = MD_lookupMeshDataId( "elem_cdl" );
	}
	~ErrorAnalysis(){}

	/*
	 * That the main function in the error analysis class. Once it has been called, user can use getLevelOfRefinement function to drive the
	 * adaptations process.
	 */
	static double getElementError(pEntity elem){
		double error = .0;
		EN_getDataDbl(elem,MD_lookupMeshDataId( "elem_error" ),&error);
		return error;
	}

	// set/get level of refinement
	static int getLevelOfRefinement(pEntity elem){
		int level;
		EN_getDataInt(elem, MD_lookupMeshDataId( "elem_ref" ), &level);
		return level;
	}

	//calculate functions
	virtual void calculate_ElementsError(pMesh theMesh, SimulatorParameters *pSimPar, GetPFuncGrad) = 0;
	virtual void calculate_SmoothedGradientNorm(pMesh, SimulatorParameters *pSimPar, GetPFuncGrad) = 0;
	virtual void calculate_SmoothedGradientNorm_excludingSingularities(pMesh, SimulatorParameters *pSimPar, GetPFuncGrad) = 0;
	virtual void calculate_CharacteristicDimensionLength(pMesh) = 0;
	virtual void calculate_DegreeOfRefinement(pMesh, int, int, bool) = 0;

	void calculate_GlobalError(pMesh, GetPFuncGrad);
	void calculate_AvgError(pMesh,double,bool);
	void initializeParameters(pMesh);

	/*! \brief: Find max element depth and max and min refinement level.
	 * \param pMesh FMDB mesh object
	 */
	void calculate_MaxMinErrorData(pMesh);

	// Inline functions
	void setMaxDepth(int depth) { maxDepth = depth; }
	void setMaxRefinementFlag(int level) { maxRefFlag = level; }
	void setMinRefinementFlag(int level) { minRefFlag = level; }

	int getMaxDepth() const { return maxDepth; }
	int getMaxRefinementFlag() const { return maxRefFlag; }
	int getMinRefinementFlag() const { return minRefFlag; }

	//!brief write on file error analysis data
	void monitoring(pMesh theMesh);

//protected:

	/*! \brief: count total number of elements with depth=0 (without children) or the same but excluding singularities. Two variables are set here
	 *          setNumFaces_excludingSingularities(numfaces_exSing) and setNumFaces(numfaces);
	 * \param pMesh FMDB mesh object
	 * \param bool says if must count all elements or excluding those on singularities regions
	 */
	void countElements(pMesh theMesh, bool excludingSingularities);

	/*! \brief:
	 * \param fine fine mesh
	 * \param coarse coarse mesh.
	 */
	double calculate_ErrorSum(pMesh);

	// set/get functions
	void setNumElements_excludingSingularities(int n) { numElements_excludingSingularities = n; }
	int getNumElements_excludingSingularities() const { return numElements_excludingSingularities; }

	void setSmoothedGradNorm(double n) { SGN = n; }
	double getSmoothedGradNorm() const { return SGN; }

	void setSmoothedGradNorm_singular(double n) { SGN_sing = n; }
	double getSmoothedGradNorm_singular() const { return SGN_sing; }

	void setNumElements(int n) { numElements = n; }
	int getNumElements() const { return numElements; }

	int isSingular(pEntity elem) const{
		int sing;
		EN_getDataInt(elem,sing_id,&sing);
		return sing;
	}
	void setElementAsSingular(pEntity elem){
		EN_attachDataInt(elem,sing_id,1);
	}

	void resetElementAsSingular(pEntity elem){
		EN_attachDataInt(elem,sing_id,0);
	}

	void setGlobalError(double ge){ globalError = ge; }
	double getGlobalError() const { return globalError; }

	void setAverageError(double ae){ averageError = ae; }
	double getAverageError() const { return averageError; }

	void setAverageError_excludingSingularities(double ae){ averageError_excludingSingularities = ae; }
	double getAverageError_excludingSingularities() const { return averageError_excludingSingularities; }

	void setElementError(pEntity elem, double error){
		EN_attachDataDbl(elem,elem_id,error);
	}

	void setLevelOfRefinement(pEntity elem, int level){
		EN_attachDataInt(elem,levelRef_id,level);
	}

	// set characteristic Dimension Length
	void setElement_CDL(pEntity elem, double cdl) {
		EN_attachDataDbl(elem,cdl_id,cdl);
	}

	// get characteristic Dimension Length
	double getElement_CDL(pEntity elem) const {
		double CDL = .0;
		EN_getDataDbl(elem,cdl_id,&CDL);
		return CDL;
	}

	bool checkMaximumNumberOfSubdivision(pMesh theMesh, const int &maxNumberOfSubdivision);

private:

	int singular;
	int numElements;
	int numElements_excludingSingularities;
	double SGN;
	double SGN_sing;

	double globalError;
	double averageError;
	double averageError_excludingSingularities;

	pMeshDataId elem_id;
	pMeshDataId levelRef_id;
	pMeshDataId sing_id;
	pMeshDataId cdl_id;

	int maxDepth;		// says how many times an element has been subdivided
	int maxRefFlag;		// says the highest element's level flag for refinement
	int minRefFlag;		// says the lowest element's level flag for refinement

	ofstream fid;		// error analysis output
};


#endif /* ERROSAnalysis_H_ */
