/*
 *  interpolateScalarFields.h
 *
 *
 *  Created by Rogerio Soares on 05/02/12.
 *  Copyright 2012 Universidade Federal de Pernambuco. All rights reserved.
 *
 *
 *  Interpolates data after h-type refinement.
 */

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "auxiliar.h"
#include "OctreeCreate2.h"

/*! \brief: Set of objects needed for interpolation functions. For any other interpolation function added which needs more information,
 * a new pointer, for example, it must be included into the InterpolationDataStruct struct
 *
 */
typedef double(*GetDblFunction)(pEntity);
typedef void(*SetDblFunction)(pEntity,double);
struct InterpolationDataStruct{
	pMesh m1;							// From:   sends interpolated data to m2
	pMesh m2;							// To:     receives interpolated data from m1
	Octree *theOctree;					// Octree_ structure for finding elements

	// Array size equal number fo fields
	GetDblFunction* pGetDblFunctions; 	// pointer to an array of function pointer. Type: get
	SetDblFunction* pSetDblFunctions;	// pointer to an array of function pointer. Type: set

	int numFields;
	bool (*isElementSpecial)(pEntity);
	int (*getLevelOfRefinement)(pEntity);
};

/*
 * List of interpolation methods
 * ----------------------------------------------------------------------------------------------------------------------------------
 */
void hRefinement(InterpolationDataStruct*);
void Linear(InterpolationDataStruct*);
void Quadratic(InterpolationDataStruct*);
void Adaptative(InterpolationDataStruct*);
void Conservative(InterpolationDataStruct*);
void PureInjection(InterpolationDataStruct*);
void HalfWeighting(InterpolationDataStruct*);
void FullWighting(InterpolationDataStruct*);


/*
 * Auxiliary function for interpolation between meshes
 */
void calculate_GeometricCoefficients(InterpolationDataStruct*, int dim);
void calculate_LinearInterpolation(InterpolationDataStruct*, int dim);
double calculate_Gradients(InterpolationDataStruct*);
void calculate_DerivativesError(InterpolationDataStruct*);
double calculate_QuadraticInterpolation(InterpolationDataStruct*);


//void interpolate(pMesh theMesh, std::list<mEntity*> &leaves, int parentDepth, GetPFuncScalar* pArray_GetFunctions, SetPFuncScalar* pArray_SetFunctions, int numFields);
//
//// after interpolation has been done, all edges flagged containing a node interpolated must be unflagged for the next simulation
//void unflagEdgesAsInterpolated(pMesh);

/*! \brief: Perform a linear interpolation for h-type adaptation.
 * \param theMesh FMDB mesh object
 * \param GetPFuncScalar pointer for an array of function pointer to get nodal scalar
 * \param SetPFuncScalar pointer for an array of function pointer to set nodal scalar
 * \param numField array length
 */
//void interpolateScalarField(pMesh theMesh, ErrorAnalysis *, AMR*, GetPFuncScalar*, SetPFuncScalar*, int);



typedef mMeshEntityContainer::iter iterall;

//enum RestrictionOperator {PURE_INJECTION, HALF_WEIGHTING, FULL_WEIGHTING};
//enum InterpolationOperator {LINEAR, QUADRATIC, ADAPTATIVE, CONSERVATIVE};

/*! \class DataTransfer
 *  \brief This class contains all the variables and procedures related to the multigrid data transfer.
 */

//class DataTransfer {
//
//  public:
//
//  /*! \brief: Constructor is called when the class is created. It creates the GeomCoeff, pLinear, NodeValues and Octree_.
//   * \param mesh1 GeomCoeff, NodeValues and Octree_ is constructed based on this mesh.
//   * \param mesh2 pLinear is constructed based on this mesh.
//   */
//  DataTransfer(pMesh mesh1, pMesh mesh2);
//  ~DataTransfer();
//
//  /*! \brief: Realize the data transfer from the fine to the coarse mesh.
//   * \param fine fine mesh
//   * \param coarse coarse mesh.
//   */
//
//  void fine_to_coarse(pMesh fine, pMesh coarse);
//
//    /*! \brief: Realize the data transfer from the coarse to the fine mesh.
//   * \param coarse coarse mesh.
//   * \param fine fine mesh
//   */
//
//  void coarse_to_fine(pMesh coarse, pMesh fine);
//
// /*! \brief: Sets the restriction method.
//   * \param RestrictionOperator type of restriction mesthod (pure injection, half weighting or full weighting).
//   */
//
//  void setRestrictionOperator(RestrictionOperator rOp){ restrictionOp = rOp; }
//
//
//   /*! \brief: Sets the interpolation method.
//   * \param InterpolationOperator type of interpolation mesthod (linear, quadratic, adaptative or conservative).
//   */
//
//  void setInterpolationOperator(InterpolationOperator iOp){ interpolationOp = iOp; }
//
//  /*! \brief: Returns the restriction method.
//   */
//  RestrictionOperator getRestrictionOperator() const{ return restrictionOp; }
//
//
//  const double* getInterpolatedData() const { return pInterpolatedVals; }
//
//   /*! \brief: Allocate the nodes value in a array.
//   * \param id It is the node id from the mesh from where the data is being transfered.
//   * \param value Value to be stored.
//   */
//  void setNodeValue(int id, double value) {pNodeValue[id] = value;}
//
//   /*! \brief: Get the desired node value.
//   * \param id It is the node id from the mesh from where the data is being transfered.
//   */
//
//  double getNodeValue(int id) {return pNodeValue[id];}
//
//    /*! \brief: Get the linear interpolation result.
//   * \param id It is the node id from the mesh to where the data is being transfered.
//   */
//
//  double getInterpolatedValues(int id) {return pInterpolatedVals[id];}
//
//  // debuggin functions
//  void showScalarField() const;
//
//  const double* getScalarField_baseMesh() const { return pNodeValue; }
//
//  bool isHalfWeighting();
//
//  private:
//
////variables
//
//    RestrictionOperator restrictionOp; 		//Restriction method
//    InterpolationOperator interpolationOp; 	//Interpolation method
//    double **pGeomCoeff;			//Geometric Coefficients matrix
//    double *pNodeValue;				//Node Values array
//
//    double *pInterpolatedVals;		//Interpolation result
//    double **pGrad;				//Gradients
//
////functions
//
//
//  /*! \brief: Assembles the Matrix and the Vectors needed on the Quadratic Interpolation.
//   * \param theMesh The mesh on which the gradients must be calculated.
//   * \param M The M-matrix (finite elements method)
//   * \param Fx The RHS (Right Hand Side) vector used to calculate du/dx.
//   * \param Fy The RHS (Right Hand Side) vector used to calculate du/dy.
//   */
//    double Assembly_Mat_Vec(pMesh theMesh,Mat M, Vec Fx, Vec Fy);
//
//  /*! \brief: Once the first order derivatives are computed, the second order derivatives are calculated as a media.
//   * \param theMesh The mesh on which the gradients must be calculated.
//   */
//
//    double calculate_SecondOrderDerivatives(pMesh theMesh);
//
//  /*! \brief: Assembles the Matrix and the Vectors needed on the Quadratic Interpolation.
//   * \param theMesh The mesh on which the gradients must be calculated.
//   * \param gradx The vector on which the first order derivatives (du/dx) are stored.
//   * \param grady The vector on which the first order derivatives (du/dy) are stored.
//   */
//
//    double store_FirstOrderDerivatives(pMesh theMesh,Vec gradx, Vec grady);
//
//  /*! \brief: Compute the global error related to the derivatives using the 2-norm.
//   * \param theMesh The mesh on which the gradients must be calculated.
//   */
//

//
//  /*! \brief: Once the interpolation is done, the neighborhood is used to weight the interpolation value.
//   * \param mesh The mesh on which the gradients must be calculated.
//   */
//
//    void addNeighborsWeighting(pMesh mesh);
//
//    /*! \brief: Returns the interpolation method.
//   */
//
//  InterpolationOperator getInterpolationOperator() const{
//	  return interpolationOp;
//  }
//
//  /*! \brief: Allocate the calculated geometric coefficients in a matrix.
//   * \param row It is the node id from the mesh to where the data is being transfered.
//   * \param column It can be 3 or 4, dependind on the dimension. It is related to the nodes of the element found with the Octree_.
//   * \param value Value to be stored.
//   */
//
//  void setGeomCoeff(int row, int column, double value) {
//	  pGeomCoeff[row][column]= value;
//  }
//
//  /*! \brief: Get the calculated geometric coefficient.
//   * \param row It is the node id from the mesh to where the data is being transfered.
//   * \param column It can be 3 or 4, dependind on the dimension. It is related to the nodes of the element found with the Octree_.
//   */
//
//  double getGeomCoeff(int row, int column) {
//	  return pGeomCoeff[row][column];
//  }
//
//    /*! \brief: Store the linear interpolation in a array.
//   * \param id It is the node id from the mesh to where the data is being transfered.
//   */
//
//  void setInterpolatedValues(int id, double value){
//	  pInterpolatedVals[id]= value;
//  }
//
//  /*! \brief: Calculate the Geometric Coefficients (Area or Volume coordinates).
//   * \param mesh mesh from where the data is being transfered.
//   * \param dim mesh dimension.
//   */
//  void calculate_GeomCoeff(pMesh mesh, int dim);
//
//  /*! \brief: Calculate linear interpolation.
//   * \param mesh mesh from where the data is being transfered.
//   * \param dim mesh dimension.
//   */
//  void calculate_LinearInterpolation(pMesh mesh,pMesh mesh1, int dim);
//
//  /*! \brief: Calculate quadratic interpolation.
//   * \param mesh mesh from where the data is being transfered.
//   * \param dim mesh dimension.
//   */

//
//  /*! \brief: Calculate gradients (used for quadratic interpolation).
//   * \param theMesh mesh from where the data is being transfered.
//   */
//

//
//  /*! \brief: Allocate the gradients in a array.
//   * \param id It is the node id from the mesh from where the data is being transfered.
//   * \param value Value to be stored.
//   */
//
//  void setGrad(int id, int column,double value) {
//	  pGrad[id][column]= value;
//  }
//
//  /*! \brief: Get the desired node gradient.
//   * \param id It is the node id from the mesh from where the data is being transfered.
//   */
//
//  double getGrad(int id, int column) {
//	  return pGrad[id][column];
//  }
//
//};
#endif // INTERPOLATESCALARFIELD_H_
