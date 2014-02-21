#include "EBFV1_hyperbolic.h"

namespace PRS{

	// Euler Foward - time advance:  (compute saturation explicitly). For each new time step compute a new saturation.
	void EBFV1_hyperbolic::calculateExplicitAdvanceInTime(double delta_T){
		saveSwField();					// save Sw field (Sw_t) before calculate Sw_t+1
		nodeWithOut_Wells(delta_T);
		nodeWith_Wells(delta_T);
	}

	// Advance time for all free node not located in wells
	void EBFV1_hyperbolic::nodeWithOut_Wells(double delta_T){
		//cout << "nodeWithOut_Wells\n";
		double Sw,Sw_old,nonvisc,volume;
		int idx;
		int nnode = pPPData->getNumNodeFreeWells();
		for(int i=0; i<nnode; i++){
			pPPData->getFreeIndex(i,idx);
			pPPData->getSaturation(idx,Sw_old);
			pPPData->getNonvisc(idx,nonvisc);
			pGCData->getVolume(idx,volume);
			Sw = Sw_old - (delta_T/(0.2*volume))*nonvisc;
			if ( Sw<1.0e-6 ) Sw = .0;
			pPPData->setSaturation(idx,Sw);
	#ifdef _SEEKFORBUGS_
			if (Sw > 1.01 || Sw < .0){
				char msg[256]; sprintf(msg,"Sw[%d/%d] = %.4f. 0.0 <= Sw <-1.0 \n",i,nnode,Sw);
				throw Exception(__LINE__,__FILE__,msg);
			}
	#endif
		}
	}

	// Advance time for nodes in production well
	void EBFV1_hyperbolic::nodeWith_Wells(double delta_T){
		//cout << "nodeWith_Wells\n";
		const int N = pSimPar->getWellTimeDiscretion(); ///  1000;								// magic number :)
		int i,j, well_idx;
		double Sw,Sw0,Sw_old,Vt,Qi,Vi,Qwi,Qt,fw,nonvisc,wp,volume, nrc, porosity;
		double dt_well = (double)delta_T/N;			// time step for well nodes saturation
		double cml_oil,Qo,Qw;

		// todo: REMOVER ESSE BACALHO!!!!
		Qt = pSimPar->getFlowrateValue(51);	// source/sink term
		Vt = pSimPar->getWellVolume(51);		// for node i, Qi is a fraction of total well flow rate (Qt)
		cml_oil = .0;
		Qo = .0;
		Qw = .0;
		int nnodes = pPPData->getNumNodesWells();
		for (i=0; i<nnodes; i++){
			pPPData->getNeumannIndex(i,well_idx);
			pPPData->getNonvisc(well_idx,nonvisc);
			pPPData->getSaturation(well_idx,Sw_old);
			pGCData->getVolume(well_idx,Vi);
			Qi = Qt*(Vi/Vt);							// Fluid (water+oil) flow rate through node i
			wp = 0.2*Vi;
			for (j=0; j<N; j++){
				Sw0 = Sw_old - (dt_well/wp)*(nonvisc);
				fw = pPPData->getFractionalFlux(Sw0);
				Qwi = fabs(fw*Qi);
				Sw = Sw0 - dt_well*(Qwi/wp);
				Sw_old = Sw;
			}
	#ifdef _SEEKFORBUGS_
			if (Sw > 1.01 || Sw < .0){
				char msg[256]; sprintf(msg,"Water saturation is out-of-bound [0 1]. Sw = %.4f\n",Sw);
			}
	#endif

			pPPData->setSaturation(well_idx,Sw);
			double fw = pPPData->getFractionalFlux(Sw);	// oil fractional flux
			double fo = pPPData->getOilFractionalFlux(Sw);	// oil fractional flux
			Qo += fabs(Qi*fo);
			Qw += fabs(Qi*fw);
			cml_oil += Qo;
		}
		setRecoveredOil(Qo/(Qo+Qw));
		cml_oil = cml_oil*delta_T + getCumulativeOil();
		setCumulativeOil(cml_oil);
	}

	void EBFV1_hyperbolic::saveSwField(){
		double Sw;
		int i, nnodes;
		pGCData->getMeshNodes(nnodes);
		for(i=0; i<nnodes; i++){
			pPPData->getSaturation(i,Sw);
			pPPData->setSw_old(i,Sw);
		}
	}
}
