/*
 * CPU_Profiling.h
 *
 *  Created on: Sep 19, 2014
 *      Author: rogerio
 */

#ifndef CPU_PROFILE_H_
#define CPU_PROFILE_H_

#include "auxiliar.h"

/* This class provides a simple way to get CPU statistics related to a specific code block.
 *
 * See example below to understand how this class works.
 *
 * CPU_Profile::Start();	// Start clock counting
 * foo_1(); 				// function to be profiled
 * CPU_Profile::End();		// Stop clock time and give a name to identify it
 *
 * CPU_Profile::Start();	// Start clock counting
 * foo_2(); 				// function to be profiled
 * CPU_Profile::End();		// Stop clock time and give a name to identify it
 *
 * CPU_Profile::Start();	// Start clock counting
 * foo_N(); 				// function to be profiled
 * CPU_Profile::End();		// Stop clock time and give a name to identify it
 *
 * // Stop profile procedure. Give a file name to where statistic will be printed
 * CPU_Profile::StatisticsOutPut("Filename.txt");
 *
 * Output:
 *
 * CPU profile: Displays statistics about CPU time
 * ===========================================================================================
 * Function                    CPU[%]                       h/m/s
 *
 * foo_1                        18                          0 53 21
 * foo_2                        38                          1 12 11
 * ...
 * foo_N                         2                          0 1 32
 *
 * User can user a counter to get an average CPU time after a set of function are called several times.
 */

class CPU_Profile{
public:

	CPU_Profile(){
	}
	~CPU_Profile(){}

	static void Start(){
		tic = MPI_Wtime();
	}

	static void End(char* whatprofiling){
		CPUtimeList cpulist = cpumap[whatprofiling];
		cpulist.push_back(MPI_Wtime() - tic);
		cpumap[whatprofiling] = cpulist;
	}

	static void StatisticOutput(const char* filename){
		int sizemap = (int)cpumap.size();
		if (!sizemap){
			cout << "CPU profiling statistics could not be printed.\n";
			return;
		}

		double cumulated, average, h, m, s, total;
		int i;
		double arrayprofile[sizemap];

		ofstream fid;
		fid.open(filename);

		i = 0;
		total = .0;
		for (CPUtimeMapIter iter1 = cpumap.begin(); iter1 != cpumap.end(); iter1++){
			CPUtimeList cpulist = (*iter1).second;
			cumulated = .0;
			//cout << (*iter1).first << ":\n";
			for (CPUtimeListIter iter2 = cpulist.begin(); iter2 != cpulist.end(); iter2++){
				cumulated += *iter2;
				//cout << cumulated << "\n";
			}
			int sizelist = (int)cpulist.size();
			//cout << "total: " << total << "\tsizelist: " << sizelist << endl;
			arrayprofile[i] = (double)cumulated/sizelist; // take average time
			total += (double)arrayprofile[i];
			//cout << "avg: " << arrayprofile[i] << endl;
			i++;
		}

		fid << "CPU profile: Displays statistics about CPU time\n";
		fid << "===========================================================================================\n";
		fid << "Function                    CPU[%]                       h/m/s\n\n";

		const int WIDTH = 31;
		int linewidth;
		i = 0;
		string funcname;
		for (CPUtimeMapIter iter1 = cpumap.begin(); iter1 != cpumap.end(); iter1++){
			//cout << "arrayprofile[i]: " << arrayprofile[i] << endl;
			convertSecToTime(arrayprofile[i],&h,&m,&s);
			funcname = (*iter1).first;
			linewidth = WIDTH - (int)funcname.length();	// (*iter1).first.size() = size of function name string
			fid << (*iter1).first <<  setw(linewidth) << setprecision(0) << fixed << 100.0*arrayprofile[i]/total << "%\t\t\t";
			fid << setprecision(0) << h << " " << m << " " << s << endl;
			i++;
			funcname.clear();
		}
		//cout << "Statistics were cleaned from memory and saved in file.\n";
		cpumap.clear();
	}

	typedef std::list<double> CPUtimeList;
	typedef std::list<double>::iterator CPUtimeListIter;
	typedef std::map<char*,CPUtimeList> CPUtimeMap;
	typedef std::map<char*,CPUtimeList>::iterator CPUtimeMapIter;

	static CPUtimeMap cpumap;						// store all CPU time for every time _whatprofiling is called
	static double tic;								// initial counting
};


#endif /* CPU_PROFILING_H_ */
