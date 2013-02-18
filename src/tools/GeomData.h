#ifndef GEOMETRICCOEFFICIENTSDATA_H_
#define GEOMETRICCOEFFICIENTSDATA_H_

#include "AttachData.h"

// NASTY THING
// I need a variable to store pressure gradients to be used into a static member function
// If this variable is a class member, it will not compile!
// Global variable <argh!>. I swear I would not do that!
//extern Matrix<double> Cij_matrix;
//extern Matrix<double> Cij_norm_matrix;
//extern Matrix<double> Dij_matrix;
//extern Matrix<double> volume_matrix;


namespace PRS
{
	typedef std::vector<double> DataArray;
	typedef map<int,DataArray> MapDataArray;
	typedef map<int,double> MapData;
	typedef map<int,double>::iterator MData_Iter;

	/*! \struct Coefficients GeomData.h
	 *  \brief Set of data that is associated to some specific mesh entity.
	 */
	struct Coefficients{
		MapDataArray Cij;
		DataArray Dij;
		MapData volume;
		double weightedvolume;

		DataArray sumIJ;

		// domains flag for Dij
		int dom1;
		int dom2;

		DataArray edgeVector;
		DataArray edgeVectorUnitary;
		double length;
		double Cij_norm;

		int numRemoteCopies;
		std::map<int,int> numRC; // number of remote copies associated to a domain
		int flag;	// may have several meanings
		bool isBoundary;
	};


	/*! \class GeomData GeomData.h
	 *  \brief GeomData is designed to set/get set of data associated to mesh entities.
	 *  The set of data is defined by the structure Coefficients above.
	 */
	class GeomData : public AttachData
	{
	public:

		GeomData();
		~GeomData();

//		/*
//		 * New set/get functions for Cij, Cij_norm, Dij, volume coefficients. More efficient.
//		 */
//		void setCij(int dom, int row, const double* grad){
//			int col = 3*dom - 3;
//			Cij_matrix(row,col) = grad[0];
//			Cij_matrix(row,col+1) = grad[1];
//			Cij_matrix(row,col+2) = grad[2];
//		}
//
//		void getCij(int dom, int row, double* grad){
//			int col = 3*dom - 3;
//			grad[0] = Cij_matrix(row,col);
//			grad[1] = Cij_matrix(row,col+1);
//			grad[2] = Cij_matrix(row,col+2);
//		}
//
//		void setCij_norm(int dom, int row, double norm){
//			Cij_matrix(row,0) = norm;
//		}
//
//		double getCij_norm(int dom, int row){
//			Cij_matrix(row,0) = norm;
//		}


        /*
         * Geometric coefficients
         * ----------------------------------------------------------------------------------
         */

		// get Dij for 3-D meshes (the old one generate memory linking)
		bool getDij(pEntity face, int dom, double *Dij);
		void setDij(pEntity face, int dom1, int dom2, double *Dij);


		void getCij(pEntity , const int &, DataArray &);
		void getCij(pEntity edge, const int &dom, double *Cij);

		double getCij_norm(pEntity , const int &);
		bool getDij(pEntity , const int &, DataArray &);
		double getVolume(pEntity , const int &);
		double getWeightedVolume(pEntity);
		void getVolume(pEntity node, int dom, double &vol);
		int getFlag(pEntity);
		void setCij(pEntity , const int &, DataArray );
		void setCij_norm(pEntity , const int &, double );
		void setDij(pEntity , const int &, const int &, const DataArray &);
		void setVolume(pEntity , const int &, const  double &);
		void setWeightedVolume(pEntity, const  double &);
		void setFlag(pEntity, int);

		void setEdgeLength(pEntity , double );
		double getEdgeLength(pEntity );

		void setEdgeVector(pEntity, dblarray);
		void getEdgeVector(pEntity, dblarray&);

		void setEdgeVec_Unitary(pEntity, dblarray);
		void getEdgeVec_Unitary(pEntity, dblarray&);


		// check if an edge belongs to domain
		// -------------------------------------------------------------------------------------
		bool nodeBelongToDomain(pEntity , const int &);
		bool edgeBelongToDomain(pEntity , const int &);
		int getDomainFlag(pEntity node);

		// check if an edge belongs to boundary domain (any domain)
		// -------------------------------------------------------------------------------------
		void set_belongsToBoundary(pEdge,bool);
		bool belongsToBoundary(pEdge);


		void setTotalReservoirVolume(double );
		double getReservoirVolume() const;

		static int getNumRemoteCopies(pEntity node){
			GeomData pGCData;
			Coefficients* pCoeffnt = pGCData.getAttachedData_pointer<Coefficients>(node);
			return pCoeffnt->numRemoteCopies;
		}
		void setNumRemoteCopies(pEntity node, int nrc);

		// set the number of copies of an entity associated to a domain
		// An edge can have several remote copies but associated to a unique domain
		// this number can be small.
		void setNumRC(pEntity, int, int);

		static int getNumRC(pMesh theMesh, pEntity ent){
			return M_numRemoteCopies(theMesh,ent);
		}

		void setMeshDim(int mdim) { dim = mdim; }
		int getMeshDim() const { return dim; }

		void setNumGEdges(int nge) { numGEdges = nge; };
		int getNumGEdges() const { return numGEdges; }

		void setSmallestEdgeLength(double sel) { smallestEdge = sel; }
		double getSmallestEdgeLength() const { return smallestEdge; }

		double getReservoirHeight() const { return reservoirHeight; }
		void setReservoirHeight(double h) { reservoirHeight = h; }

		// used to verify coefficients summation
		void setSumIJ(pEntity,DataArray);

		// used to verify coefficients summation
		void getSumIJ(pEntity,DataArray&);

	private:
		double reservoirHeight;
		double reservoirVolume;
		int dim;
		double smallestEdge;
		int numGEdges;	// number or global edges
	};
}

#endif /*GEOMETRICCOEFFICIENTSDATA_H_*/
