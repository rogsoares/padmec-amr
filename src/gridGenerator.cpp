#include "includes.h"
#include "libincludes.h"

gridGenerator::gridGenerator(int getpoints){

    if (getpoints){
            ifstream opt("parameters/options.dat", ios::binary);

                if (opt) {

                    char buffer[11] ; //offset 11 columns("MAXGRIDS = ") in the file,then get the value

                            opt.read (buffer,11);// read data as a block
                            opt>>maxgrids;

                    opt.close();

                }
    }
}


grid1D* gridGenerator::generateGrids(int n,int*MAXGRIDS){

    int N=n+1;//Total number of points...
    grid1D* newgrids;
    double L=1; //reference...1 is 100% of the domain...
    double h=L/N;//constant //distance between points...
    int level=2; //we begin in the first coarse grid,but we need adjustments so we can use the original grid(finest),so level=2
    double multiplier=0;
    //n_current is a vector that print the number of inner points of each level (from the highest to the lowest) in a file
    int* n_current= new int[maxgrids]; //current number of inner points in the level

    filebuf buffer;
    ostream output(&buffer);
    istream input(&buffer);
    buffer.open ("parameters/coarseData_innerP.dat", ios::in | ios::out | ios::trunc);

    n_current[level-2]=n; //just to be accepted in the next line(check  if(n_current[level-2]==n) )

        while(level<=maxgrids && n_current[level-2]>4 && n_current[level-2]%2!=0){

                if(n_current[level-2]==n)   // the  assumptions in the beginning(n=n_finegrid and level=2) will now be adjusted
                    level=1;//go to the first coarse grid...

            n_current[level-1]=0;//reset number of coarse grid inner points
            multiplier=(pow(2,level))*h;//current distance multiplier  2^level*h

            n_current[level-1]=(L/multiplier)-1; //point determination n=(L/h)-1
            output<<"level "<<level<<" = "<<n_current[level-1]<<endl;// write current level's points in the file

            ++level;// go to the next level in the grid chain
        }

    input.clear();           // clear  eofbit and  failbit
    buffer.close();

    newgrids=new grid1D[level-1];
    *MAXGRIDS=level-1;

    delete[] n_current;

    return newgrids;

}

void gridGenerator::getCoarsePoints(grid1D* pGrid){

    ifstream opt("parameters/coarseData_innerP.dat", ios::binary);
    int j=10; //determines buffer lenght (begins with 10,because ("level j = ") contains 10 elements before the integer the program wants
    char* buffer;
            if (opt) {

                    for (int i=0;i<maxgrids;i++){

                    //there's a limit of 999 grids read... check the switch condition
                        switch(i){
                            case 0 :
                                buffer=new char[j]; //offset 10+j columns("level j+k = ") in the file,then get the value
                                break;
                            case 10 :
                                delete[] buffer;
                                ++j;
                                buffer=new char[j];
                                break;
                            case 100:
                                ++j;
                                delete[] buffer;
                                buffer=new char[j];
                                break;

                         }
                                opt.read (buffer,j);// read data as a block
                                opt>>pGrid[i].n_coarse;
                    }

                opt.close();

            }

    delete[] buffer;

}
