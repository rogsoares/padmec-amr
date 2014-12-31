/*
 * SMOOTHER(PRE+POST) CLASS
 *
 * smoother.h
 *
 *  Created on: 18/12/2014
 *      Author: Julio Cezar
 */

#ifndef _SMOOTHER_
#define _SMOOTHER_

class smoother{
    private :
        int PreIts ;//number of pre iterations
        int PostIts ;//number of post iterations
    public :
        smoother();
        int PreSmoothing(Mat A, Vec y, Vec u,bool* firstuse);
        int PostSmoothing(Mat A, Vec y, Vec u);
        ///read about what to do when arrive at the last grid...post or pre? or maybe none(direct solver)...


};

#endif // _SMOOTHER_


