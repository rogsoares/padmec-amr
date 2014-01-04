///*
// * Restart.cpp
// *
// *  Created on: 23/07/2009
// *      Author: rogerio
// */
//
//#include "Restart.h"
//
//namespace PRS{
//
//	void restartSimulation(SimulatorParameters *pSimPar, PhysicPropData *pPPData, pMesh theMesh){
//
//		ifstream fid;
//		string filename(pSimPar->getRestartFilename());
//
//		// open vtk file and check if operation has succeeded
//		fid.open(filename.c_str());
//		if ( !fid.is_open() ){
//			char msg[256]; sprintf(msg,"Restart file '%s' required but it could not be opened or do not exist.\n ",pSimPar->getRestartFilename().c_str());
//			throw Exception(__LINE__,__FILE__,msg);
//		}
//
//		/*
//		 * Get file step
//		 * VTK file names obey the following rule:
//		 *
//		 * case_name__rank-of-size__step-n.vtk
//		 *
//		 * where:
//		 * 	case_name - string name given by user for case problem identification
//		 *  rank-of-size - rank is the processor ID and size the number of processors used
//		 *  step-n - n is an integer number for the step (here, 'step' is not the time-step)
//		 *
//		 *  Ex.: FiveSpot-3D-anisotropic-heterogeneous__0-of-1__step-12.vtk
//		 *  case_name = "FiveSpot-3D-anisotropic-heterogeneous"
//		 *  rank-of-size = 0-of-1 (It means that only one processor was used)
//		 *  step-n = step-12 (12th step)
//		 */
//		string::size_type pos = filename.find("step-",0);
//		char buffer[256];
//		string::size_type length = filename.copy(buffer,filename.length() - pos,pos+5);
//		buffer[length]='\0';
//		pSimPar->setStepOutputFile( atoi(buffer) );
//
////		printf("VTK step: %d\n",pSimPar->getStepOutputFile());
////		throw 1;
//
//		//cout << caseName << endl;
//
//		// "Elapsed simulation time (EST) based on last VTK file:
//		double step = pSimPar->getStepOutputFile();
//		double PVI_acc = step*pSimPar->getPVIincrement();
//		pSimPar->setPVIaccumulated( PVI_acc );
//		double EST = PVI_acc*pSimPar->getSimTime();
//		pSimPar->setAccumulatedSimulationTime( EST );
//
//		/*
//		 * Open simulation-monitor dat file to extract the sum of time steps when
//		 * the last vtk file was generated, i.e, at EST. Extract also CPU time
//		 * accumulated.
//		 */
//		char monitorCase[256]; sprintf(monitorCase,"%s_simulation-monitor-%d.dat",pSimPar->getOutputPathName().c_str(),P_size());
//		ifstream fid2; fid2.open(monitorCase);
//		if ( !fid2.is_open() ) failOpeningFile(monitorCase,__LINE__,__FILE__);
//		string line, CPU_time;
//		int count = 0;
//		char buffer2[256]; sprintf(buffer2,"Sum. tSteps: %.5f",EST);
//		cout << "buffer2: " << buffer2 << endl; //throw 1;
//		if (P_pid()==-1){
//			do{
//				getline(fid2,line);
//
//				// number of time-steps
//				pos = line.find("Sum. tSteps:");
//				if (pos!=string::npos) count++;
//
//				// accumulated CPU-time
//				pos = line.find("accumulated:");
//				if (pos!=string::npos) CPU_time = line;
//
//				cout << "line: " << line << endl;
//
//			}while ( line.compare(buffer2) );
//		}
//		pos = CPU_time.find(":",0);
//		char buffer3[256];
//		length = CPU_time.copy(buffer3,CPU_time.length() - pos,pos+1);
//		buffer3[length]='\0';
//		double cpu = (double)atof(buffer3);
//
//		cout << "cpu: " << cpu << endl;
//
//		pSimPar->setCPU_time(0);
//		pSimPar->setTStepNumber(0);
//
//		/*
//		 * Read saturation field from file and use it as input data to define
//		 * initial condition.
//		 */
//		pEntity node;
//		double Sw;
//
//		// set pointer file position and then start reading pressure data
//		pSimPar->setPositionToRead(fid,"SCALARS Saturation float 1");
//		getline(fid,line); // read one more line to skip "LOOKUP_TABLE default"
//		VIter vit = M_vertexIter(theMesh);
//		while ( (node = VIter_next(vit)) ){
//			fid >> Sw;
//			pPPData->setSaturation(node,Sw);
//		}
//		//throw 1;
//	}
//}
