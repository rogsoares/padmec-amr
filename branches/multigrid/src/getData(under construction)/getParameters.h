/*
 * GET PARAMETERS CLASS
 *
 * getParameters.h
 *
 *  Created on: 18/12/2014
 *      Author: Julio Cezar
 */
#ifndef _GETPRMTS_
#define _GETPRMTS_
#include "includes.h"

class getParameters{
	private :
		static int maxgrids; // max number of levels in a chain
		static int RestricionOperator; //0-injection 1-Full-Weighting
		static int Pre;   //number of pre-smoothing iterations
		static int Post;   //number of post-smoothing iterations
		static double D0;
		static double DL;
	public :
        getParameters(const char* action);
        int getData(const char* parameterName);

};

#endif
