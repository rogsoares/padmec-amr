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

void makeMeshCopy2(pMesh m, PADMEC_mesh *pm, double(*pGetPressure)(pEntity), double(*pGetSaturation)(pEntity)){
	pm->numVertices = M_numVertices(m);
	pm->numElements = M_numFaces(m);
	pm->ID = new int[M_numVertices(m)];
	pm->coords = new double[3*M_numVertices(m)];
	pm->elements = new int[3*M_numFaces(m)];
	pm->field1 = new double[M_numVertices(m)];
	pm->field2= new double[M_numVertices(m)];

	int i;
	int k = 0;
	int iter = 0;
	double coord[3];
	pEntity ent;
	VIter vit = M_vertexIter(m);
	while ( (ent = VIter_next(vit)) ){
		pm->ID[iter] = EN_id(ent);
		V_coord(ent,coord);
		for(i=0;i<3;i++){
			pm->coords[k] = coord[i];
		//	cout << pm->coords[k] << "\t";
			k++;
		}
		//cout << endl;
		pm->field1[iter] = pGetPressure(ent);
		pm->field2[iter] = pGetSaturation(ent);
		iter++;
	}
	VIter_delete(vit);

	k = 0;
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		for(i=0;i<3;i++){
			pm->elements[k] = EN_id((mVertex*)ent->get(0,i));
			k++;
		}
	}
	FIter_delete(fit);
}

void makeMeshCopy2(PADMEC_mesh* pm,pMesh m, void(*pSetPressure)(pEntity,double), void(*pSetSaturation)(pEntity,double)){
	int i,j, k = 0;
	double x,y,z;
	for (i=0; i<pm->numVertices; i++){
		x = pm->coords[k++];
		y = pm->coords[k++];
		z = pm->coords[k++];
	//	cout << "x = " << x << "\ty = " << x << "\tz = " << z << endl;
		mVertex* v = m->createVertex(pm->ID[i],x,y,z,0);
		pSetPressure((pEntity)v,pm->field1[i]);
		pSetSaturation((pEntity)v,pm->field2[i]);
	}

	k = 0;
	mVertex* vertices[3];
	for (i=0; i<pm->numElements; i++){
		vertices[0] = m->getVertex( pm->elements[k++] );
		vertices[1] = m->getVertex( pm->elements[k++] );
		vertices[2] = m->getVertex( pm->elements[k++] );
	//	cout << EN_id(vertices[0]) << "\t" << EN_id(vertices[1]) << "\t" << EN_id(vertices[2]) << endl;
		m->createFaceWithVertices(vertices[0],vertices[1],vertices[2],0);
	}
	m->modifyState(2,0);
	//m->modifyState(0,2);
	//STOP();
}

void makeMeshCopy(pMesh m1, pMesh m2, 
				  void(*pSetPressure)(pEntity,double), double(*pGetPressure)(pEntity),
				  void(*pSetSaturation)(pEntity,double), double(*pGetSaturation)(pEntity)){
	mEntity* ent;
	mVertex* v;
	mVertex* vertices[4];
	pGEntity pGEnt;
	double coord[3];
	int i,j,k,ID[4];
	double val;

	VIter vit = M_vertexIter(m1);
	while ( (ent = VIter_next(vit)) ){
		pGEnt = 0;
		if (ent->getClassification()){
			pGEnt = ent->getClassification();
		}
		V_coord(ent,coord);
		mVertex* v = m2->createVertex(EN_id(ent),coord[0],coord[1],coord[2],pGEnt);
		val = pGetPressure(ent);
		pSetPressure((pEntity)v,val);
		val = pGetSaturation(ent);
		pSetSaturation((pEntity)v,val);
		
		//cout << pGetPressure(v) << "\t" << pGetSaturation(v) << endl;
	}
	VIter_delete(vit);
	//STOP();

	EIter eit = M_edgeIter(m1);
	while ( (ent = EIter_next(eit)) ){
		pGEnt = 0;
		if (ent->getClassification()){
			pGEnt = ent->getClassification();
		}
		m2->createEdge((mVertex*)ent->get(0,0), (mVertex*)ent->get(0,1),pGEnt);
	}
	EIter_delete(eit);

	if (m1->getDim()==2){
		FIter fit = M_faceIter(m1);
		while ( (ent = FIter_next(fit)) ){
			pGEnt = 0;
			if (ent->getClassification()){
				pGEnt = ent->getClassification();
			}
			m2->createFaceWithVertices((mVertex*)ent->get(0,0), (mVertex*)ent->get(0,1),(mVertex*)ent->get(0,2),pGEnt);
		}
		FIter_delete(fit);
	}
	else{
	}

#ifdef _SEEKFORBUGS2_
/*	pEntity ent1, ent2;
	int IDs1[4],IDs2[4];

	VIter vit1 = M_vertexIter(m1);
	VIter vit2 = M_vertexIter(m2);
	while ( (ent1 = VIter_next(vit1)) ){
		ent2 = VIter_next(vit2);
		printf("%d [%d]\t%d [%d]\n",EN_id(ent1),GEN_tag(ent1->getClassification()),EN_id(ent2),GEN_tag(ent2->getClassification()));
	}
	VIter_delete(vit1);
	VIter_delete(vit2);
	//exit(1);
	EIter eit1 = M_edgeIter(m1);
	EIter eit2 = M_edgeIter(m2);
	while ( (ent1 = EIter_next(eit1)) ){
		ent2 = EIter_next(eit2);
		printf("%d %d [%d]\t %d %d [%d]\n",
				EN_id((mVertex*)ent1->get(0,0)),EN_id((mVertex*)ent1->get(0,1)),GEN_tag(ent1->getClassification()),
				EN_id((mVertex*)ent2->get(0,0)),EN_id((mVertex*)ent2->get(0,1)),GEN_tag(ent2->getClassification()));
	}
	EIter_delete(eit1);
	EIter_delete(eit2);
	//exit(1);

	// let's check if copy was made successfully
	if (m1->getDim()==2){
		FIter fit1 = M_faceIter(m1);
		FIter fit2 = M_faceIter(m2);
		while ( (ent1 = FIter_next(fit1)) ){
			ent2 = FIter_next(fit2);
			getTriVerticesIDs(ent1,IDs1);
			getTriVerticesIDs(ent2,IDs2);
			printf("%d %d %d [%d]\t%d %d %d [%d]\n",
					IDs1[0],IDs1[1],IDs1[2],GEN_tag(ent1->getClassification()),
					IDs2[0],IDs2[1],IDs2[2],GEN_tag(ent2->getClassification()));
			//			printf("%d %d %d [%d]\t%d %d %d [%d]\n",
			//								IDs1[0],IDs1[1],IDs1[2],0,IDs2[0],IDs2[1],IDs2[2],0);
		}
		FIter_delete(fit1);
		FIter_delete(fit2);


		// let's delete m1 and see if m2 still alive
		FIter fit = M_faceIter(m1);
		while ( (ent = FIter_next(fit)) ){
			m1->DEL(ent);
		}
		FIter_delete(fit);

		EIter eit = M_edgeIter(m1);
		while ( (ent = EIter_next(eit)) ){
			m1->DEL(ent);
		}
		EIter_delete(eit);

		VIter vit = M_vertexIter(m1);
		while ( (ent = VIter_next(vit)) ){
			m1->DEL(ent);
		}
		VIter_delete(vit);

		cout << "m1     m2\n";
		cout << "V: " << M_numVertices(m1);
		cout << "\tV: " << M_numVertices(m2) << endl;
		cout << "E: " << M_numEdges(m1);
		cout << "\tE: " << M_numEdges(m2) << endl;
		cout << "F: " << M_numFaces(m1);
		cout << "\tF: " << M_numFaces(m2) << endl;


		throw Exception(__LINE__,__FILE__,"Finish him!\n");


	}
	else{
		throw Exception(__LINE__,__FILE__,"To be implemented soon...\n");
		//		RIter rit = M_regionIter(m1);
		//		while ( (ent = RIter_next(rit)) ){
		//			v[0] = (mVertex*)ent->get(0,0);
		//			v[1] = (mVertex*)ent->get(0,1);
		//			v[2] = (mVertex*)ent->get(0,2);
		//			v[3] = (mVertex*)ent->get(0,3);
		//			pGEnt = ent->getClassification();
		//			m2->createTetWithVertices(v[0],v[1],v[2],v[3],pGEnt);
		//		}
		//		RIter_delete(rit);
	}*/
#endif
}
void deleteMesh(PADMEC_mesh* pm){
	delete[] pm->ID; pm->ID = 0;
	delete[] pm->coords; pm->coords = 0;
	delete[] pm->elements; pm->elements = 0;
	delete[] pm->field1; pm->field1 = 0;
	delete[] pm->field2; pm->field2 = 0;
}

void deleteMesh(pMesh m){
	pEntity ent;
	cout << "Deleting mesh: \n";
	if (m->getDim()==3){
		cout << "Deleting tetrahedras\n";
		RIter rit = M_regionIter(m);
		while ( (ent = RIter_next(rit)) ){
			m->DEL(ent);
		}
		RIter_delete(rit);
	}

	cout << "Deleting triangles\n";
	FIter fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		m->DEL(ent);
	}
	FIter_delete(fit);

	cout << "Deleting edges\n";
	EIter eit = M_edgeIter(m);
	while ( (ent = EIter_next(eit)) ){
		m->DEL(ent);
	}
	EIter_delete(eit);
}


/*
 * 
 * seta as entidade das malhas com os flags necessarios para uma nova simulação (condições de contorno). 
 * isso deve vir ja do adaptador. mas , enquanto isso nao é resolvido por saulo, vamos "adpatando" as coisas.
 * Rogério: 23/05/2013
 * 
 * ---------------------------------------------------------------------------------------------------------------
 * ATENÇÃO: ESTAMOS CONSIDERANDO UM FIVE-SPOT HOMOGENEO, OU SEJA, O MESMO FLAG 3300 PARA TODOS OS TRIANGULOS,
 *          FLAG 2000 PARA AS ARESTAS DO CONTORNO E 10/51 PARA OS POÇOS !!!!!!!
 * * -------------------------------------------------------------------------------------------------------------
 */
void PADMEC_GAMBIARRA(pMesh m){
	pEntity ent;
	
	m->modifyState(2,1);
	m->modifyState(1,2);
	m->modifyState(2,0);
	m->modifyState(0,2);

	// insere flag nas aresta do  contorno da geometria
 	EIter eit = M_edgeIter(m);
// 	while ( (ent = EIter_next(eit)) ){
// 	  if ( E_numFaces(ent)==1 )
// 		cout << GEN_tag( ent->getClassification() ) << endl;
// 	}
// 	EIter_delete(eit);
	
	std::list<pEntity> belist;
	int c1 = 0, c2 = 0;
	eit = M_edgeIter(m);
	while ( (ent = EIter_next(eit)) ){
	  if ( E_numFaces(ent)==1 ){
		ent->classify( m->getGEntity(2000,1) );
		c1++;
		belist.push_back(ent);
	  }
	  c2++;
	}
	EIter_delete(eit);
	
	
// 	eit = M_edgeIter(m);
// 	while ( (ent = EIter_next(eit)) ){
// 	  if ( E_numFaces(ent)==1 ){
// 	    cout << GEN_tag( ent->getClassification() ) << endl;
// 	  }
// 	}
// 	EIter_delete(eit);
	
	cout << " All edges: " << c2 << "\t\tBoundary Edges = " << c1 << endl;
	
	if (!belist.size()){
		throw Exception(__LINE__,__FILE__,"Num of boundary edges NULL!\n");
	}

	// insere flag nos nós dos poços
	// primeiro temos que encontar estes nós que estão nas quinas inferior esquerda e superior direita
	std::list<pEntity>::iterator iter = belist.begin();
	mVertex *v1, *v2;
	double coord1[3], coord2[3];
	double xmax = 0.0, xmin = 1e3;
	double ymax = 0.0, ymin = 1e3;
	
	mVertex* injectionWell;
	mVertex* productionWell;
	for (; iter!=belist.end(); iter++){
	  // pegue coordenadas dos nos e guarde a que tiver menor/maior valor para x e y.
	  ent = *iter;
	  v1 = (mVertex*)ent->get(0,0);
	  v2 = (mVertex*)ent->get(0,1);
	  V_coord(v1,coord1);
	  V_coord(v2,coord2);
	  
	  
	  
	  if (coord1[0] < xmin){
	    xmin = coord1[0];
	    injectionWell = v1;
	  }
	  if (coord1[1] < ymin){
	    ymin = coord1[1];
	  //  injectionWell = v1;
	  }

	  if (coord1[0] > xmax){
	    xmax = coord1[0];
	    //productionWell = v1;
	  }
	  if (coord1[1] > ymax){
	    ymax = coord1[1];
	    //productionWell = v1;
	  }

	  if (coord2[0] < xmin){
	    xmin = coord2[0];
	    //injectionWell = v2;
	  }
	  if (coord2[1] < ymin){
	    ymin = coord2[1];
	   // injectionWell = v2;
	  }

	  if (coord2[0] > xmax){
	    xmax = coord2[0];
	    //productionWell = v2;
	  }
	  if (coord2[1] > ymax){
	    ymax = coord2[1];
	    //productionWell = v2;
	  }
	  
	  if (coord1[0] >= xmax && coord1[1] >= ymax){
	    productionWell = v1;
	  }
	  if (coord2[0] >= xmax && coord2[1] >= ymax){
	    productionWell = v2;
	  }

	  
	  if (coord1[0] <= xmin && coord1[1] <= ymin){
	    injectionWell = v1;
	  }
	  if (coord2[0] <= xmin && coord2[1] <= ymin){
	    injectionWell = v2;
	  }
	  
	  
	}
	
	injectionWell->classify( m->getGEntity(10,0) );
	productionWell->classify( m->getGEntity(51,0) );
	cout << "min: " << xmin <<"\t"<< ymin << "injectionWell ID: " << EN_id(injectionWell) << "  flag " << GEN_tag( injectionWell->getClassification() ) << endl;
	cout << "max: " << xmax <<"\t"<< ymax << "productionWell ID: " << EN_id(productionWell) << "  flag " << GEN_tag( productionWell->getClassification() ) << endl;
	belist.clear();
	//STOP();
	
	// insere flag nos elementos
 	FIter fit = M_faceIter(m);
// 	while ( (ent = FIter_next(fit)) ){
// 		cout << GEN_tag( ent->getClassification() ) << endl;
// 	}
// 	FIter_delete(fit);
	
	fit = M_faceIter(m);
	while ( (ent = FIter_next(fit)) ){
		ent->classify( m->getGEntity(3300,2) );
	}
	FIter_delete(fit);
	
// 	fit = M_faceIter(m);
// 	while ( (ent = FIter_next(fit)) ){
// 		cout << GEN_tag( ent->getClassification() ) << endl;
// 	}
// 	FIter_delete(fit);
}
