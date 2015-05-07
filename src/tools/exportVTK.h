#ifndef EXPORT_VTK_H_
#define EXPORT_VTK_H_

#include "PhysicPropData.h"
#include "ErrorAnalysis.h"
#include "GeomData.h"


void exportSolutionToVTK(pMesh,void*,void*,void*,void*,string);
void printVTKHeader(ofstream &, pMesh);
void printVerticesCoordenates(ofstream &, pMesh);
void printElementConnectivities(ofstream&, pMesh, int, int);
void printElementConnectivities(ofstream&, pEntity, int);
void printCellTypeList(ofstream &, int, int);
void printPressure(ofstream &, pMesh, PRS::PhysicPropData*);
void printSaturation(ofstream &, pMesh, PRS::PhysicPropData*);
void printSaturation(ofstream &, pMesh, PRS::PhysicPropData*, PRS::GeomData*);
void printPressureGradient(ofstream& fid, PRS::GeomData* pGCData, PRS::PhysicPropData *pPPData);

//void printElementError(ofstream &, pMesh, ErrorAnalysis*);
void printElementError(ofstream&, GeomData*, ErrorAnalysis*);
void print_h_ratio(ofstream &, GeomData*, ErrorAnalysis*);
void print_singular_regions(ofstream &, GeomData*, ErrorAnalysis*);

// New VTK source file
// --------------------------------------------------------------
void print_headers(ofstream&,int);
void print_coodinates(ofstream&,const GeomData*);
void print_connectivities(ofstream&, const GeomData*);
void print_celltype(ofstream&,int,int);
void print_pressure(ofstream&, const PhysicPropData*);
void print_saturation(ofstream&, const PhysicPropData*);


#endif /*EXPORT_VTK_H_*/
