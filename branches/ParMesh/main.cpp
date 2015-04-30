#include "mesh.h"

int main(int argc, char** argv){
	MPI_Init(&argc, &argv);

	MeshDB::Mesh* mesh = new MeshDB::Mesh;

	if (argc!=4){
		cout << "\nERROR\n\nYou MUST type: executable mesh_file(IN) #ref_level mesh_file(OUT)\nExiting...\n";
		exit(1);
	}

	mesh->read(argv[1]);
	mesh->createEdgeDataStructure();

	// codigo vem aqui.
	//mesh->quad_preprocessor();

	//mesh->refine_mesh( atoi(argv[2]) );
	//mesh->write(argv[3]);

	mesh->printMeshStatistic(argv[3]);
	MPI_Finalize();
	return 0;
}
