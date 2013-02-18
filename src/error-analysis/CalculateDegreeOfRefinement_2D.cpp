/*
 * DefineDegreeOfRefinement.cpp
 *
 *  Created on: 01/03/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis_2D.h"

int getLevelToRefine(const double& ratio);
int getLevelToUnefine(const double& ratio);

/*
 * Define degree of refinement for all elements
 */
void ErrorAnalysis_2D::calculate_DegreeOfRefinement(pMesh theMesh, int maxNumSubdivision, int numSubdivision_perStep, bool excludingSingularities){
#ifdef __ERROR_ANALYSIS_DEBUG__
	bool EADebugging = false;
#endif

	double avgError = this->getAverageError();
	cout << setprecision(8) << scientific;
	int refunref_level;
	FIter fit = M_faceIter (theMesh);
	while (pEntity face = FIter_next(fit)){
		if (!theMesh->getRefinementDepth(face)){
			// get element error
			double error = this->getElementError(face);

			// element characteristic length
			double h_old = this->getElement_CDL(face);

			// define a new characteristic element dimension (h_new1 = d1 na tese)
			//eq. 4.27 e 3.2.2 (p‡g. 32) - tese de Filipe
			double h_new = h_old*(avgError/error);

			//cout << error << "  " << avgError << "  " << h_old << "  " << h_new << endl;

			/*
			 * NOTE:
			 * The following code lines were taken from the error analysis written by Felipe Araujo and update by Bruno Luna in (mar/2006).
			 * The code is in Fortran.
			 */
			double ratio = h_new/h_old;
			// Unrefine procedure
			if(ratio >= 1.01){
				refunref_level = getLevelToUnefine(ratio);
			}
			else{											// Refine procedure
				refunref_level = getLevelToRefine(ratio);
				// if level is the maximum, then it's on a singular region
				if (refunref_level == 4){
					this->setElementAsSingular(face);
				}
			}
			//cout << "refunref_level: " << refunref_level << endl;
			/* Do not allow root element to be refined more than the maximum number of refinement set by user. If root element's depth + leave's
			 * level to refine is greater than maximum number of subdivisions, then correct refunref_level. If it happens, error analysis MUST
			 * NOT be performed once more, because tolerance will not be satisfied due the presence of the singularities.
			 */
			int rootDepth = theMesh->getRefinementDepth(face->root());
			if ( rootDepth + refunref_level > maxNumSubdivision ){
				///cout << "refunref_level: " << refunref_level << endl;
				refunref_level = maxNumSubdivision - rootDepth;
			}
			else{
				// Before set the new element level or refinement, check if the previous one is less or greater the new level.
				// Greater level MUST always prevail
				refunref_level = std::max(refunref_level,this->getLevelOfRefinement(face));
			}

			// stops element to be refined completely saving memory.
			if (refunref_level > numSubdivision_perStep){
				refunref_level = numSubdivision_perStep;
			}
			this->setLevelOfRefinement(face,refunref_level);

#ifdef __ERROR_ANALYSIS_DEBUG__
			if (refunref_level){
				EADebugging = true;
			}
#endif
		}
	}
	FIter_delete(fit);

#ifdef __ERROR_ANALYSIS_DEBUG__
	if (!EADebugging){
		//throw Exception(__LINE__,__FILE__,"Global error is greater than tolerance but any element has been flagged to be refined.\n");
	}
#endif
}

void ErrorAnalysis_2D::calculate_CharacteristicDimensionLength(pMesh theMesh){
	mVertex* v;
	double edIJ[3][3];
	Trellis_Util::mPoint pts[3];

	FIter fit = M_faceIter (theMesh);
	while (pEntity face = FIter_next(fit)){
		if (!theMesh->getRefinementDepth(face)){
			for (int i=0; i<3; i++){
				v = (mVertex*)(face->get(0,i));
				pts[i] = v->point();
			}
			for (int i=0; i<3; i++) {
				edIJ[0][i] = pts[0](i)-pts[1](i);
				edIJ[1][i] = pts[0](i)-pts[2](i);
				edIJ[2][i] = pts[1](i)-pts[2](i);
			}
			double dimLen = .0;
			for (int i=0; i<3; i++) dimLen += sqrt(pow(edIJ[i][0],2)+pow(edIJ[i][1],2));
			double CDL = (double)(dimLen/3.0);
			this->setElement_CDL(face,CDL);
		}
	}
	FIter_delete(fit);
}

int getLevelToRefine(const double& ratio){
	if ( ratio >= 1 && ratio < 1.2 ){
		return 0;
	}
	else if (ratio<1 && ratio>=0.5){
		return 1;
	}
	else if (ratio<0.5 && ratio>=0.25){
		return 2;
	}
	else if (ratio<0.25 && ratio>=0.125){
		return 3;
	}
	else{
		return 4;
	}
}

int getLevelToUnefine(const double& ratio){
	if(ratio > 4){
		return -3;
	}
	else if (ratio <= 4){
		return -2;
	}
	else if (ratio > 2.5){
		return -1;
	}
}
