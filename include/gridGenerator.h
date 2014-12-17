#ifndef _GGEN_H
#define _GGEN_H
#include "grid1D.h"

class gridGenerator{
    private:
        int maxgrids;
    public:
        gridGenerator(int getpoints);//constructor
        grid1D* generateGrids(int n,int *MAXGRIDS);//return an array of grids
        void getCoarsePoints(grid1D* pGrid);

};
#endif

