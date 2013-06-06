#include "ErrorAnalysis.h"
#include <string>


class ErrorEstimator
{
public:
	// Contructor and Destructor
	ErrorEstimator(pMesh theMesh);
	~ErrorEstimator();


	map<int,float> VCL;
	set<pVertex> AllNodes;
	set<pFace> ListedFaces;

	void VertexList (pMesh theMesh);
	void LeitorDeLista(ErrorAnalysis *pErrorAnalysis, pMesh theMesh);

private:
	char meshFilename[];

};

