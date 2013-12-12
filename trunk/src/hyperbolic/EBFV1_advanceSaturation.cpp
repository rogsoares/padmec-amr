#include "EBFV1_hyperbolic.h"

namespace PRS{
	
	// Euler Foward - time advance:  (compute saturation explicitly). For each new time step compute a new saturation.
	double EBFV1_hyperbolic::calculateExplicitAdvanceInTime(pMesh theMesh, double delta_T){
		double startt = MPI_Wtime();
		saveSwField(theMesh);			// save Sw field (Sw_t) before calculate Sw_t+1
		nodeWithOut_Wells(theMesh,delta_T);
		nodeWith_Wells(theMesh,delta_T);
		double endt = MPI_Wtime();
		return endt-startt;
	}
	
	// Advance time for all free node not located in wells
	int EBFV1_hyperbolic::nodeWithOut_Wells(pMesh theMesh, double delta_T){
		cout << "nodeWithOut_Wells\n";
		double Sw,Sw_old,nonvisc,wp,volume,phi;
		pEntity node;
		
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int id = EN_id(node);
			int well_flag = GEN_tag(node->getClassification());
			
			if ( Sw != 1.0 ){
				if ( !pStruct->pSimPar->isInjectionWell(well_flag) && !pStruct->pSimPar->isProductionWell(well_flag) ){
					Sw_old = pStruct->pPPData->getSaturation(node);
					nonvisc = pStruct->pPPData->getNonViscTerm(node);
					wp = .0;
					for (SIter dom=pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end();dom++){
						volume = pGCData->getVolume(node,*dom);
						phi = pStruct->pSimPar->getPorosity(*dom);
						wp += volume*phi;
					}
					// calculate new saturation value
					Sw = Sw_old - (delta_T/wp)*nonvisc;
					if ( Sw<1.0e-6 ) Sw = .0;
					
// 					if (id == 22){
// 						printf("delta_T: %f\twv: %f\tnonvisc: %f\n",delta_T,wp,nonvisc);
// 					}
					
					pStruct->pPPData->setSaturation(node,Sw);
					
					if (Sw > 1.01 || Sw < .0){
						char msg[256]; sprintf(msg,"Water saturation is out-of-bound [0 1]. Sw = %.10E on node %d.  Sw_old = %f\n",Sw,id,Sw_old);
						throw Exception(__LINE__,__FILE__,msg);
					}
				}
			}
		}
		VIter_delete(vit);
		return 0;
	}
	
	
	// Advance time for nodes in production well
	// =========================================================================
	void EBFV1_hyperbolic::nodeWith_Wells(pMesh theMesh, double delta_T){
		cout << "nodeWith_Wells\n";
		const int N = pStruct->pSimPar->getWellTimeDiscretion(); ///  1000;								// magic number :)
		int node_ID,well_flag,i;
		double Sw,Sw0,Sw_old,Vt,Vi,Qi,Qwi,Qt,fw,nonvisc,wp;
		double dt_well = (double)delta_T/N;			// time step for well nodes saturation
		double oilFlow = .0;
		pVertex node;
		
		if ( pStruct->pSimPar->rankHasProductionWell() ){
			
			// FOR EACH PRODUCTION WELL
			map<int,set<int> >::iterator miter;
			for ( miter=pStruct->pSimPar->mapNodesOnWell.begin(); miter!=pStruct->pSimPar->mapNodesOnWell.end(); miter++){
				well_flag = miter->first;
				if ( pStruct->pSimPar->isProductionWell(well_flag) ){
					// source/sink term
					Qt = pStruct->pSimPar->getFlowrateValue(well_flag);
					
					// for node i, Qi is a fraction of total well flow rate (Qt)
					Vt = pStruct->pSimPar->getWellVolume(well_flag);
					
					// FOR EACH NODE ON PRODUCTION WELL
					//int II = 0;
					for (SIter siter = miter->second.begin(); siter!=miter->second.end();siter++){
						node_ID = *siter;
						node = (mEntity*)theMesh->getVertex( node_ID );
						Sw_old = pStruct->pPPData->getSaturation(node);
						nonvisc = pStruct->pPPData->getNonViscTerm(node);
						int geomFlag = pGCData->getDomainFlag(node);
						wp = pGCData->getVolume(node,geomFlag)*pStruct->pSimPar->getPorosity(geomFlag);
						
						Vi = .0;
						double Vi = .0, volume, nrc;
						for (SIter_const dom=pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end(); dom++){
							pVertex node = (mEntity*)theMesh->getVertex( node_ID );
							volume = pGCData->getVolume(node,*dom);
							nrc = 1.0;//pGCData->getNumRemoteCopies(node) + 1.0;
							Vi += volume/nrc;
						}
						
						#ifdef _SEEKFORBUGS_
						if ( fabs(wp)==0.0 || fabs(Vi)==0.0 || fabs(Vt)==0.0 || fabs(Qt)==0.0 )
							//throw Exception(__LINE__,__FILE__,"Volume cannot be null.\n");
							#endif
							
							// Fluid (water+oil) flow rate through node i
							Qi = Qt*(Vi/Vt);
						
						// FOR 'N' WELL TIME-STEPS
						// =========================================================================
						for (i=0; i<N; i++){
							Sw0 = Sw_old - (dt_well/wp)*(nonvisc);
							fw = pStruct->pPPData->getFractionalFlux(Sw0);
							Qwi = fabs(fw*Qi);
							Sw = Sw0 - dt_well*(Qwi/wp);
							if ( Sw < 1.0e-6 ) Sw = .0;
							Sw_old = Sw;
						}
						
						if (Sw > 1.01 || Sw < .0){
							char msg[256]; sprintf(msg,"Water saturation is out-of-bound [0 1]. Sw = %.10E on node %d.  Sw_old = %f\n",Sw,node_ID,Sw_old);
						}
						
						pStruct->pPPData->setSaturation(node,Sw);
						
						// Oil production
						double Qoi = fabs(Qi) - fabs(fw*Qi);
						oilFlow += fabs(Qoi);
					}
					setRecoveredOilValue(oilFlow);
				}
			}
		}
	}
	
	void EBFV1_hyperbolic::saveSwField(pMesh theMesh){
		pEntity node;
		VIter vit = M_vertexIter(theMesh);
		while ( (node = VIter_next(vit)) ){
			int id = EN_id(node);
			double Sw_old = pStruct->pPPData->getSaturation(node);
			pStruct->pPPData->setSaturation_Old(node,Sw_old);
//			cout << "Sw_t = " << pStruct->pPPData->getSaturation_Old(node) << endl;
		}
		VIter_delete(vit);
	}
}
