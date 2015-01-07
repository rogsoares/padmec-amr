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
		_refLevel = refLevel;
		mesh_partition();
		mesh_distribution();
		cout << "Refinement level: " << endl;
		for(int i=0; i<refLevel; i++){
//			cout << "Before refinement: \n";
//			calculate_volume();
			createEdgeDataStructure();
			cout << "                  " << i+1 << "/" << refLevel << endl;
			switch (getElemType()){
			case TRI:
				refine_TRI();
				break;
			case QUAD:
				refine_QUAD();
				break;
			case TETRA:
				refine_TRI();	// boundary element faces
				refine_TETRA();
				break;
			}
//			cout << "After refinement: \n";
//			calculate_volume();
			deleteEdgeDataStructure();
		}
		bdry_linkSetup();
	}
}
