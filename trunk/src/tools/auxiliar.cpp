#include "auxiliar.h"

void failOpeningFile(string filename, int line, const char* cppfile){
	char msg[256]; sprintf(msg,"File '%s' could not be opened or do not exist.\n ",filename.c_str());
	throw Exception(line,cppfile,msg);
}

double F_area(const double *ptn1, const double *ptn2, const double *ptn3){

	double xyz[3][3], n[3];
	for (int i=0;i<3;i++){
		xyz[0][i] = ptn1[i];
		xyz[1][i] = ptn2[i];
		xyz[2][i] = ptn3[i];
	}
	double a[3] = { xyz[1][0]-xyz[0][0], xyz[1][1]-xyz[0][1], xyz[1][2]-xyz[0][2] } ;
	double b[3] = { xyz[2][0]-xyz[0][0], xyz[2][1]-xyz[0][1], xyz[2][2]-xyz[0][2] } ;
	computeCrossProduct(a,b,n);

	return .5*norm(n);
}

double norm(const double *n){
	return sqrt( pow(n[0],2) + pow(n[1],2) + pow(n[2],2) );
}


void getEdgeVector(pEdge edge, dblarray &edgVec){
	double coord1[3];
	double coord2[3];

	V_coord(edge->get(0,0),coord1);			// extract coordenates from ondes
	V_coord(edge->get(0,1),coord2);			// extract coordenates from ondes

	double signal = (EN_id(edge->get(0,0)) < EN_id(edge->get(0,1)))?1.0:-1.0;
	for (int i=0; i<3; i++) edgVec[i] = signal*(coord2[i] - coord1[i]);
}

double E_length(pEdge edge){
	dblarray v(3,.0);
	getEdgeVector(edge,v);
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// multiply a tensor K by a vector vec
void multTensorVector2D(const double *tensor, const double *vec, double *mult){
	mult[0] = tensor[0]*vec[0]+tensor[1]*vec[1];
	mult[1] = tensor[2]*vec[0]+tensor[3]*vec[1];
}

// converts a string read from file to a double
double strToDouble(string &str){
	return strtod(getSubString(str), 0);
}

// converts a string read from file to an interger
int strToInteger(string &str){
	return atoi(getSubString(str));
}

// filters string the numeric part - internal use
const char* getSubString(string &str){
	string::size_type loc = str.find( "=", 0 );
	string numberstr = str.substr(loc+1, str.size()-loc);
	return numberstr.c_str();
}

// Make a vector given to nodes
void makeVector(const double *A, const double *B, double *v){
	for (int i=0; i<3; i++) v[i] = B[i] - A[i];
}

bool isEdgeExternal(int flag){
	return ( flag>=1000 && flag<=1499 );
}

bool isEdgeInternal(int flag){
	return ( flag>=1500 && flag<=1999 );
}

bool isFaceExternal(int flag){
	return ( flag>=2000 && flag<=2499 );
}

bool isFaceInternal(int flag){
	return ( flag>=2500 && flag<=2999 );
}


double* getIdentityMatrix(const int &dim){
	int i,j,k = 0;
	double *Identity = new double[dim*dim];
	for (i=0; i<dim; i++)
		for (j=0; j<dim; j++)
			Identity[k++] = (i==j)?1.0:.0;
	return Identity;
}

int getVertexFlag(pVertex node){
	return getEntityFlag(0,node);
}

int getEdgeFlag(pEdge edge){
	return getEntityFlag(1,edge);
}

int getFaceFlag(pFace face){
	return getEntityFlag(2,face);
}

int getTetraFlag(pRegion tetra){
	return getEntityFlag(3,tetra);
}

int getEntityFlag(int i, pEntity ent){
	pGEntity gEntity;
	switch (i)
	{
	case 0:
		gEntity = (pGEntity)V_whatIn(ent);
		break;
	case 1:
		gEntity = (pGEntity)E_whatIn(ent);
		break;
	case 2:
		gEntity = (pGEntity)F_whatIn(ent);
		break;
	case 3:
		gEntity = (pGEntity)R_whatIn(ent);
		break;
	}
	return (gEntity)?GEN_tag(gEntity):0;
}

double getSmallestEdgeLength(pMesh theMesh){
	pEntity edge;
	double delta_x = 1e30; // infinity delta_x

	// loop over all edges to get the smallest one
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		double length = E_length(edge);
		if (delta_x>length)	delta_x = length;
	}
	EIter_delete(eit);
	//printf("smallest length = %f\n",delta_x);
	return P_getMinDbl(delta_x);
}

void replaceAllOccurencesOnString(string &theString, string::size_type size,
		string seekFor, string replaceFor){
	string::size_type pos = 0;
	while (true){
		pos = theString.find(seekFor,pos);
		if (pos != string::npos)
			theString.replace(pos,size,replaceFor);
		else
			break;
	}
}

void getIJnodes(pEdge edge, std::vector<pEntity> &vertex){
	M_GetVertices(edge,vertex);
	if (EN_id(vertex[0])>EN_id(vertex[1])) std::swap(vertex[0],vertex[1]);
}

void F_getEdges(pMesh theMesh, pEntity face, std::vector<pEntity> &edges){
	std::vector<pVertex> vertex;
	M_GetVertices(face,vertex);
	int IDs[3]={EN_id(vertex[0]),EN_id(vertex[1]),EN_id(vertex[2])};
	edges[0] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[0]),
			theMesh->getVertex(IDs[1]));
	edges[1] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[0]),
			theMesh->getVertex(IDs[2]));
	edges[2] = (pEntity)theMesh->getEdge(theMesh->getVertex(IDs[1]),
			theMesh->getVertex(IDs[2]));

	vertex.clear();
}

VPoint E_getVertexPoint(pEntity edge, int i){
	return ((mVertex*)edge->get(0,i))->point();
}


void E_getVerticesCoord(pEntity edge, double *I, double *J){
	VPoint point_I = ((mVertex*)edge->get(0,0))->point();
	VPoint point_J = ((mVertex*)edge->get(0,1))->point();
	for(int i=0; i<3; i++){
		I[i] = point_I(i);
		J[i] = point_J(i);
	}
}

void E_IJvector(double *I, double *J, double *IJ){
	for (int i=0; i<3; i++) IJ[i] = J[i]-I[i];
}

void E_vertices(pEntity edge, std::vector<pEntity> & vertex){
	vertex[0] = (pEntity)edge->get(0,0);
	vertex[1] = (pEntity)edge->get(0,1);
}

void printSimulationHeader(){
	if (!P_pid())
		std::cout<< "\n\n\t\t\t###################################################################\n"
		"\t\t\t\t\tUNIVERSIDADE FEDERAL DE PERNAMBUCO\n"
		"\t\t\t\t\tPADMEC - DECIV - DEMEC - PRH-26 - CNPq\n"
		"\t\t\t\t\t\t\t2008-2010\n"
		"\n\t\t\t\t\tParallel Reservoir Simulator - v1.1.5\n"
		"\t\t\t\t\tAuthor: Silva, R. S.\n"
		"\t\t\t\t\tNO WARRANTY!!! \n"
		"\t\t\t\t\tUSE IT AT YOUR OWN RISK!!!\n"
		"\t\t\t###################################################################\n\n";
}


double ratio(double a, double b);

//double dot(const double* a, const double *b){
//	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
//}

double inner_product(const dblarray &a1, const dblarray &a2){
	return std::inner_product(a1.begin(),a1.end(),a2.begin(),.0);
}

double ratio(double a, double b){
	return (double)(a/b);
}

void unitary_vector(const dblarray &a1, dblarray &a2){
	std::copy(a1.begin(),a1.end(),a2.begin());
	dblarray norm(a1.size(),sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) ));
	std::transform(a2.begin(),a2.end(),norm.begin(),a2.begin(),ratio);
}

void unitary_vector(dblarray &a1){
	dblarray norm(a1.size(),sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) ));
	std::transform(a1.begin(),a1.end(),norm.begin(),a1.begin(),ratio);
}

double norm2(const dblarray &a1){
	return sqrt( std::inner_product(a1.begin(),a1.end(),a1.begin(),.0) );
}

double dot(const double *u, const double *v, int dim){
	double _dot = .0;
	for (int i=0;i<dim;i++){
		_dot += u[i]*v[i];
	}
	return _dot;
}

double norm_L2(const double *p, int dim){
	return sqrt( dot(p,p,dim) );
}

void setBarrier(){
	MPI_Barrier(MPI_COMM_WORLD);
	throw 1;
}

void checklinepassing(int line, const char* filepath, std::string somemessage){
	static ofstream fid1;
	static bool openonce = true;
	static int count = 0;

	if (openonce){
		openonce=false;
		char filename[256]; sprintf(filename,"Debbugging-for-rank__%d.txt",P_pid());
		fid1.open(filename);
	}

	fid1 << "File	: " << filepath << std::endl;
	fid1 << "Line	: " << line << std::endl;
	fid1 << "Counter	: " << ++count << std::endl;
	fid1 << "Msg		: " << somemessage << std::endl;
	fid1 << std::endl << std::endl;
}

void convertSecToTime(double t, double *h, double *m, double *s){
	double frac;
	frac = modf(t/3600.,h);
	frac = modf(frac*60.,m);
	frac = modf(frac*60.,s);
}


void LogFiles(LOG_FILES LG, double t1, double t2, double timeStep, double accSimTime,
		      string path, bool hasRestart, int last_step, double CPU_time){

	static int count = (hasRestart)?last_step:0;
	static double acc_CPUtime = (hasRestart)?CPU_time:0;
	static ofstream fid1;
	static ofstream fid2;

	switch (LG){
	case OPENLG:{
		// print results on file
		if (!P_pid()){
			char fname1[256]; sprintf(fname1,"%s_simulation-monitor-%d.dat",path.c_str(),P_size());
			if (hasRestart){
				fid1.open(fname1,ios_base::app);
				fid1 << "* * * * * * * * * * * * * * * * * * * * * * * * *\n"
						"*                Restart required               *\n"
						"* * * * * * * * * * * * * * * * * * * * * * * * *\n\n\n";
			}
			else{
				fid1.open(fname1);
				fid1 << "SIMULATION PROGRESS MONITOR\nNumber of processors: "<<P_size()<<"\n\n\n";
			}
		//	cout << __LINE__ << endl;
			char fname2[256]; sprintf(fname2,"%s_CPU-process-time-%d.xls",path.c_str(),P_size());
			fid2.open(fname2);
			fid2 << "Elliptic Hyperbolic total-step Total(accumulated)\n";
		//	cout << __LINE__ << endl;
			std::cout << "\n--------------------------------------------------\n"
			"Start Simulation\n"
			"--------------------------------------------------\n\n";
		}
	}
	break;
	case UPDATELG:{
		// take average time from all processors
		double total = P_getSumDbl(t1+t2)/((double)P_size());

		// rank 0 is in charge to print the output CPU-time
		if (!P_pid()){
			double h, m, s;

			acc_CPUtime += total;

			fid1.precision(0);
			fid1 << "          CPU-time[sec]  percentual\n";
			fid1 << "------------------------------------------------\n";
			convertSecToTime(t1,&h,&m,&s);
			fid1 << "elliptic   : " << fixed << h <<"(h)  "<< m <<"(m)  "<< s <<"(s) " << 100.0*(t1/total) <<"%\n";
			convertSecToTime(t2,&h,&m,&s);
			fid1 << "hyperbolic : " << h <<"(h)  "<< m <<"(m)  "<< s <<"(s) " << 100.0*(t2/total) <<"%\n";
			convertSecToTime(t1+t2,&h,&m,&s);
			fid1 << "ellip+hyper: " << h <<"(h)  "<< m <<"(m)  "<< s <<"(s)\n";
			convertSecToTime(acc_CPUtime,&h,&m,&s);
			fid1 << "accumulated: " << h <<"(h)  "<< m <<"(m)  "<< s <<"(s)\n\n";
			cout << "Accumulated Simulation time: " << h <<"(h)  "<< m <<"(m)  "<< s <<"(s)\n\n";
			fid1.precision(5);
			fid1 << "Step       : " << ++count  <<"\n";
			fid1 << "timeStep   : " << timeStep <<"\n";
			fid1 << fixed << "Sum. tSteps: " << accSimTime  <<"\n\n\n";

			char cString[256]; sprintf(cString,"%f %f %f %f\n",t1,t2,t1+t2,acc_CPUtime);
			string theString(cString);
			replaceAllOccurencesOnString(theString,1,".",",");
			fid2 << theString;
		}
	}
	break;
	case CLOSELG:{
		if (!P_pid()){
			fid1.close();
			fid2.close();
		}
	}
	}
}

void STOP(){
	MPI_Barrier(MPI_COMM_WORLD); throw 1;
}

int printMatrixToFile(Mat& m,const char* filename){
	 PetscViewer viewer;
	 PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	 MatView(m,viewer);
	 PetscViewerDestroy(viewer);
	 return 0;
}

int printVectorToFile(Vec& v,const char* filename){
	 PetscViewer viewer;
	 PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer);
	 VecView(v,viewer);
	 PetscViewerDestroy(viewer);
	 return 0;
}

void getEdgesChildren(pEdge edge, mEdge** edgesChildren){
	std::list<mEntity*> edgeLeaves;
	std::list<mEntity*>::iterator iter;
	edge->getLeaves(edgeLeaves);
	iter = edgeLeaves.begin();
	edgesChildren[0] = (mEdge*)(*iter);
	iter++;
	edgesChildren[1] = (mEdge*)(*iter);
	edgeLeaves.clear();
}

void getFacesAroundFace(pFace face, pFace *faces, int &numFaces){
	numFaces = 0;
	// If a face has one or two boundary edge(s), this face will have two or one neighbor(s).
	for (int i=0; i<3; i++){
		pEdge edge = face->get(1,i);
		// if true, take face's neighbor
		if (E_numFaces(edge) == 2){
			if (edge->get(2,0)==face)
				faces[numFaces] = edge->get(2,1);
			else
				faces[numFaces] = edge->get(2,0);
			numFaces++;
		}
	}
}

int EN_getFlag(pEntity ent){
	if (!ent->getClassification()){
		pEntity root = ent->root();
		if (!root){
			throw Exception(__LINE__,__FILE__,"Null root?!\n");
		}
		return (!root->getClassification())?0:GEN_tag( root->getClassification() );
	}
	else{
		return GEN_tag( ent->getClassification() );
	}
}

void getTriVerticesIDs(pFace face, int* IDs){
	IDs[0] = EN_id( face->get(0,0) );
	IDs[1] = EN_id( face->get(0,1) );
	IDs[2] = EN_id( face->get(0,2) );
}

void printFaceIDs(pEntity face){
	int IDs[3];
	getTriVerticesIDs(face,IDs);
	cout << "Face IDs: " << IDs[0] << " " << IDs[1] << " "<< IDs[2] << "\n";
}

void getNeighboursFace(pEntity face, std::vector<pEntity> &neighbousFace, int &size){
	if (!face){
		throw Exception(__LINE__,__FILE__,"Null pointer!.");
	}
	int k = 0;
	for (int i=0; i<3; i++){
		pEdge edge = (pEdge)face->get(1,i);
		if (E_numFaces(edge)==2){
			if (edge->get(2,0) != face){
				neighbousFace[k++] = edge->get(2,0);
			}
			else{
				neighbousFace[k++] = edge->get(2,1);
			}
		}
	}
//	for (int i=0; i<k; i++){
//		printFaceIDs(neighbousFace[i]);
//	}
//	if (k==2){
//		cout << "Face is on boundary\n";
//	}
	size = k;
}
