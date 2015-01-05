/*
 Original element: I-J-K
 New elements: I-R-T, S-R-T, R-S-J, S-T-K
                     I
 * *
 *   *
 *     *
 *       *
                R - - - - T
 * \       / *
 *   \     /   *
 *     \   /     *
 *       \ /       *
           J * * * * S * * * * K
 */

#include "mesh.h"

namespace MeshDB{
	void  Mesh::refine_mesh(int refLevel){
		cout << "Refinement level: " << endl;
		for(int i=0; i<refLevel; i++){
			cout << "Before refinement: \n";
			this->calculate_volume();
			createEdgeDataStructure();
			cout << "                  " << i+1 << "/" << refLevel << endl;
			switch (this->getElemType()){
			case TRI:
				refine_TRI();
				break;
			case QUAD:
				refine_QUAD();
				break;
			case TETRA:
				refine_TETRA();
				break;
			}
			cout << "After refinement: \n";
			this->calculate_volume();
			deleteEdgeDataStructure();
		}
		printMeshStatistic();
	}
}
