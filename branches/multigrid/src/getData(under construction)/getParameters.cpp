/*
 * getParameters.cpp
 *
 *  Created on: 18/12/2014
 *      Author: Julio Cezar
 */

#include "includes.h"
#include "libincludes.h"
#include <map>




getParameters::getParameters(const char* action){



    if(action=="GET"){
        //start
        cout<<"\ngetting data..\n\n";

        /// options.dat
        ifstream file1("MG_parameters/options.dat", ios::binary);
        char buffer[256]; //offset columns

            if (file1) {
                file1.read (buffer,11);
                    file1>>getParameters::maxgrids;cout<<"\nmaxgrids : "<<maxgrids;
                file1.read (buffer,23);
                    file1>>SetOperator::Rop;cout<<"\nRO : "<<RestricionOperator;
                file1.read (buffer,9);
                    file1>>Pre;cout<<"\nPre : "<<Pre;
                file1.read (buffer,10);
                    file1>>Post;cout<<"\nPost : "<<Post;
                file1.close();
            }

        cout<<"\n\noptions.dat ...ok\n";

        /// boundaries.dat
         ifstream file2("MG_parameters/boundaries.dat", ios::binary);

        if (file2) {
                file2.read (buffer,13);
                    file2>>D0;
                file2.read (buffer,13);
                    file2>>DL;
                file2.close();
            }

        cout<<"\n\nboundaries.dat ...ok\n";
    }
    else{
        getData(action);
    }
}


int getParameters::getData(const char* parameterName){

std::map <string, int>parameter;

parameter["MAXGRIDS"]=1;

/*
    switch(parameter){
    ///options.dat parameters
        case 1 :
            return maxgrids;
        break;

        case 'RESTRICTIONONOPERATOR' :
            return RestricionOperator;
        break;

        case 'PRE' :
            return Pre;
        break;

        case 'POST':
            return Post;
        break;

     ///boundaries.dat parameters
        case 'BOUNDARY0' :
            return D0;
        break;

        case 'BOUNDARYL' :
            return DL;
        break;

    ///invalid parameter
        default :
            cout<<"\n\ninvalid parameter : "<<parameterName<<endl<<endl;
            return 0;
    }*/
}
