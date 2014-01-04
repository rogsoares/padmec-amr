#ifndef GEOMETRICCOEFFICIENTSDATA_H_
#define GEOMETRICCOEFFICIENTSDATA_H_

#include "AttachData.h"
#include "Matrix.h"


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

        /*
         * Geometric coefficients
         * ----------------------------------------------------------------------------------
         */

		// get Dij for 3-D meshes (the old one generate memory linking)
		bool getDij(pEntity face, int dom, double *Dij);
		void setDij(pEntity face, int dom1, int dom2, double *Dij);


		void getCij(pEntity , const int &, DataArray &);
		void getCij(pEntity edge, const int &dom, double *Cij);


		bool getDij(pEntity , const int &, DataArray &);
		double getVolume(pEntity , const int &);
		double getWeightedVolume(pEntity);
		void getVolume(pEntity node, int dom, double &vol);
		int getFlag(pEntity);
		void setCij(pEntity , const int &, DataArray );
		void setCij_norm(pEntity , const int &, double );
		double getCij_norm(pEntity , const int &);
//		void setCij_norm(int dom, int row, double norm){
//			Cij_norm[dom].setValue(row,norm);
//		}
//
//		double getCij_norm(int dom, int row) const{
//			return Cij_norm[dom].getValue(row);
//		}

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

		void getCij(int dom, int row, double* cij);
		void setCij(int dom, int row, double* cij);

		void calculateNumEdges(pMesh, int, int*);
//		void calculateNumBDRYEdges(pMesh, int, int*);
//		void calculateNumNodes(pMesh, int, int*);
		void allocatePointers(int);
		void deallocatePointers(int);

		int getNumEdgesPerDomain(int i) const{
			return numDomEdges[i];
		}

	private:
		double reservoirHeight;
		double reservoirVolume;
		int dim;
		double smallestEdge;
		int numGEdges;	// number or global edges

		int* numDomEdges;	// number of edges per domain
		Matrix<double>* Cij;		// Cij vector
		//Matrix<double>* Cij_norm;
	};
}

#endif /*GEOMETRICCOEFFICIENTSDATA_H_*/
