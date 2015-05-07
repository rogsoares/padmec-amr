///*
// * ErrorAnalysis-2D.h
// *
// *  Created on: 24/04/2012
// *      Author: rogsoares
// */
//
//#ifndef ERRORANALYSIS_2D_H_
//#define ERRORANALYSIS_2D_H_
//
//#include "ErrorAnalysis.h"
//
//
///*! \class ErrorAnalysis_2D
// *  \brief
// */
//class ErrorAnalysis_2D : public ErrorAnalysis{
//public:
//
//	ErrorAnalysis_2D(){}
//	~ErrorAnalysis_2D(){}
//
//private:
//	/*! \brief:
//	 * \param
//	 * \param
//	 */
//	void calculate_ElementsError(SimulatorParameters* pSimPar, GeomData*, void(*)(FIELD,int,int,int,double*),FIELD);
//
//	/*! \brief:
//	 * \param
//	 * \param
//	 */
//	void calculate_SmoothedGradientNorm(SimulatorParameters *pSimPar, GeomData*, void(*)(FIELD,int,int,int,double*), FIELD);
//
//	/*! \brief:
//	 * \param
//	 * \param
//	 */
//	void calculate_SmoothedGradientNorm_Singularity(SimulatorParameters*, GeomData*, void(*)(FIELD,int,int,int,double*), FIELD);
//
//	/*! \brief:
//	 * \param
//	 * \param
//	 */
//	void calculate_DegreeOfRefinement(GeomData*, RefinementStrategies, bool, FIELD);
//	void calculateElemHeightRatio(GeomData*, double);
//	void weightNewElementHeight(bool,std::list<pEntity> &);
//};
//
//#endif /* ERRORANALYSIS_2D_H_ */
