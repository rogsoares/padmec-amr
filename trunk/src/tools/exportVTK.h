#ifndef EXPORT_VTK_H_
#define EXPORT_VTK_H_

#include "PhysicPropData.h"
#include "ErrorAnalysis.h"


void exportSolutionToVTK(pMesh,void*,void*,void*,string);
void printVTKHeader(ofstream &, pMesh);
void printVerticesCoordenates(ofstream &, pMesh);
void printElementConnectivities(ofstream&, pMesh, int, int);
void printElementConnectivities(ofstream&, pEntity, int);
void printCellTypeList(ofstream &, int, int);
void printPressure(ofstream &, pMesh, PRS::PhysicPropData*);
void printSaturation(ofstream &, pMesh, PRS::PhysicPropData*);
void printDegreeOfRefinement(ofstream &, pMesh, ErrorAnalysis*);
void printElementError(ofstream &, pMesh, ErrorAnalysis*);
void printCharac_Lenth(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis);
void print_Swgrad(ofstream &fid, pMesh theMesh, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_pw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_Sw_GradientNorm(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);
void print_Sw_GradientNorm2(ofstream &fid, pMesh theMesh, ErrorAnalysis *pErrorAnalysis, PRS::SimulatorParameters *pSimPar, PRS::PhysicPropData *pPPData);


#endif /*EXPORT_VTK_H_*/
