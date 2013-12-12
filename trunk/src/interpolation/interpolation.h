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

#include "SimulatorParameters.h"
#include "Matrix.h"
#include "OctreeCreate2.h"

/*! \brief: Set of objects needed for interpolation functions. For any other interpolation function added which needs more information,
 * a new pointer, for example, it must be included into the InterpolationDataStruct struct
 */
typedef double(*GetDblFunction)(pEntity);
typedef void(*SetDblFunction)(pEntity,double);
typedef mMeshEntityContainer::iter iterall;

struct InterpolationDataStruct{
	pMesh m1;							// To: receives interpolated data from m1
	pMesh m2;							// From: sends interpolated data to m2
	Octree *theOctree;					// Octree_ structure for finding elements

	// Array size equal number fo fields
	GetDblFunction* pGetDblFunctions; 	// pointer to an array of function pointer. Type: get
	SetDblFunction* pSetDblFunctions;	// pointer to an array of function pointer. Type: set

	int numFields;
	bool (*isElementSpecial)(pEntity);
	int (*getLevelOfRefinement)(pEntity);

	Matrix<double> pInterpolatedVals;			// Interpolation result
	Matrix<double> pGeomCoeff;		  // geometric coefficients are calculated in mesh from (m2) and stored into mesh to (m1)
	Matrix<double> pGrad;						// Gradients
	Matrix<double> pNodeValue;					// Node Values array
};

void interpolation(InterpolationDataStruct*, INTERPOLATION_OPTIONS);

/*
 * List of interpolation methods
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
double calculate_Gradients(InterpolationDataStruct*, int);
void calculate_DerivativesError(InterpolationDataStruct*);
double calculate_QuadraticInterpolation(InterpolationDataStruct*, int);
void finalize(InterpolationDataStruct* pIData);
void initialize(InterpolationDataStruct* pIData);

/*! \brief: Assembles the Matrix and the Vectors needed on the Quadratic Interpolation.
 * \param theMesh The mesh on which the gradients must be calculated.
 * \param M The M-matrix (finite elements method)
 * \param Fx The RHS (Right Hand Side) vector used to calculate du/dx.
 * \param Fy The RHS (Right Hand Side) vector used to calculate du/dy.
 */
double Assembly_Mat_Vec(InterpolationDataStruct* pIData, int field, Mat M, Vec Fx, Vec Fy);

// solve system of equations: Ax=y
double KSP_solver(Mat A, Mat pcMatrix, Vec y, Vec x);

/*! \brief: Once the first order derivatives are computed, the second order derivatives are calculated as a media.
 * \param theMesh The mesh on which the gradients must be calculated.
 */
double calculate_SecondOrderDerivatives(InterpolationDataStruct* pIData);

/*! \brief: Assembles the Matrix and the Vectors needed on the Quadratic Interpolation.
 * \param theMesh The mesh on which the gradients must be calculated.
 * \param gradx The vector on which the first order derivatives (du/dx) are stored.
 * \param grady The vector on which the first order derivatives (du/dy) are stored.
 */
double store_FirstOrderDerivatives(InterpolationDataStruct* pIData, Vec gradx, Vec grady);

#endif
