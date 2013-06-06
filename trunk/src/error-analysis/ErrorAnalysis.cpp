#include "ErrorAnalysis.h"

/*
 Steps for error analysis (taken from Felipe Araujo Dissertation, p47):

 step 1: Gradient Recovery (skip this step because gradients have already been calculated.)
 step 2: Calculate L2 norms of all nodal gradients and elements errors.
		 The error analysis works always favorable to security, it means that when evaluating element error related to two or more variables (e.g pressure and saturation)
		 the highest error will be considered.


 step 3: Calculate global error for solution Neta

         if (global_error > tolerance)

 step 4: calculate mean error based on Neta_1 for whole mesh
 step 5: calculate d1 (new element dimension) parameter for all mesh elements. Defines the degree of refinement for all elements.


 step 6: if (d1 < d_min) d1 = d_min
 step 7: recalculate mean error based on the new Neta_2 for whole domain excluding the singularities

 Singularities must be understood as regions where gradients are too high which leads the error analysis to set those elements with a very high refinement level

 step 8: calculate d2 (new element dimension) parameter for all mesh elements excluding elements with singularities.
 step 9: if (d2 < d_min) d2 = d_min
 step 10: d_new = min(d1,d2)
          In fact, step 10 is inside step 9, because as long the new element degree
 */

bool calculate_ErrorAnalysis(ErrorAnalysis *pEA, pMesh theMesh, SimulatorParameters *pSimPar, double tol1, double tol2, GetPFuncGrad* pGetFuncArray, int numScalarFields){
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tIN\n";
#endif

	// for remeshing use
	pEA->_pSimPar = pSimPar;
	/*
	 * Says if it will be necessary perform a mesh adaptation
	 */
	bool adapt = false;

	/*
	 * Every new error analysis, error and level of refinement of all mesh elements are set to avoid mistakes from previous analysis.
	 */
	pEA->initializeParameters(theMesh);

	/*
	 * Before perform the error analysis, check if some original element number of subdivisions reached the maximum number of subdivision allowed.
	 * If true, return 0
	 */
	if ( pEA->checkMaximumNumberOfSubdivision(theMesh,pSimPar->getMax2D()) ){
#ifdef __ERROR_ANALYSIS_DEBUG__
	cout << "Mesh has reached maximum number of element refinement allowed by user.";
#endif
		return adapt;
	}

	/*
	 *  Perform an error estimation for each scalar field. At the end of each evaluation, the new element
	 *  dimension size will be determined for that can lead to a better safety.
	 */
	cout << setprecision(8) << scientific << fixed;
	for (int field=0; field<numScalarFields; field++){
		pEA->calculate_ElementsError(theMesh,pSimPar,pGetFuncArray[field]);										// step 2:
		pEA->calculate_SmoothedGradientNorm(theMesh,pSimPar,pGetFuncArray[field]);								// step 3:
		pEA->calculate_GlobalError(theMesh,pGetFuncArray[field]);
		pEA->countElements(theMesh,false);// step 4:
		cout << "  G. Error: " << pEA->getGlobalError() << endl;
		cout << "      tol1: " << tol1 << endl;
		if ( pEA->getGlobalError() > tol1 ){
			pEA->calculate_AvgError(theMesh,tol1,false);														// step 5:
			pEA->calculate_CharacteristicDimensionLength(theMesh);
			pEA->calculate_DegreeOfRefinement(theMesh,pSimPar->getMax2D(),pSimPar->getNumSubdivision_perStep(),false);								// step 6:
			/*
			 * Repeat calculations for all mesh elements excluding those not belonging to high gradient regions
			 */

			pEA->calculate_SmoothedGradientNorm_excludingSingularities(theMesh,pSimPar,pGetFuncArray[field]);	// step 7:
			pEA->calculate_GlobalError(theMesh,pGetFuncArray[field]);											// step 8:
			pEA->countElements(theMesh,true);
			cout << "  G. Error: " << pEA->getGlobalError() << endl;
			cout << "      tol2: " << tol2 << endl;
			pEA->calculate_AvgError(theMesh,tol2,true);															// step 9:
			pEA->calculate_DegreeOfRefinement(theMesh,pSimPar->getMax2D(),pSimPar->getNumSubdivision_perStep(),true);													// step 10:
			adapt = true;
		}
		else{
			cout << "Adaptation not necessary!\n--------------------------\n";
		}
		cout << "Num elements: " << pEA->getNumElements() << endl;
	}
	pEA->monitoring(theMesh);
#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: Error Analysis\tOUT\n";
#endif
	return adapt;
}
