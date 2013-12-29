#include "EBFV1_hyperbolic.h"

namespace PRS{

#ifdef _SEEKFORBUGS_
static bool check_IAT = false;
#endif //_SEEKFORBUGS_

// computes total velocity on edges for domain dom
double EBFV1_hyperbolic::calculateIntegralAdvectiveTerm(pMesh theMesh, int dom){
	cout << "calculateIntegralAdvectiveTerm\n";
	int dim = pGCData->getMeshDim();
	double startt = MPI_Wtime();

	pEntity edge;
	dblarray Cij(dim,.0), vel(dim,0);
	double dt = 1e+10, Sw_I, Sw_J, fwII, fwJJ, fwIJ, df_dsIJ, n, alpha;
	double  non_visc_fv, non_visc_ad, nonvisc_I, nonvisc_J;
	double courant, phi, length;
	
	// loop over edges
	EIter eit = M_edgeIter(theMesh);
	while ( (edge = EIter_next(eit)) ){
		if ( pGCData->edgeBelongToDomain(edge,dom) ){
			pGCData->getCij(edge,dom,Cij);

			// get nodes I and J
			pEntity I = (pVertex)edge->get(0,0);
			pEntity J = (pVertex)edge->get(0,1);

			// todo: colocar este codigo em algum luar de forma que so seja feito uma unica vez.
			if (EN_id(edge->get(0,0)) > EN_id(edge->get(0,1))){
				std::swap(I,J);
			}

			// considering up-wind approximation for saturation field
			Sw_I = pStruct->pPPData->getSaturation(I);
			Sw_J = pStruct->pPPData->getSaturation(J);

			// get high order approximation for saturation
			if ( pStruct->pSimPar->useHOApproximation() ){
				pHOApproximation->getSw_HighOrderApproximation(edge,dom,Sw_I,Sw_J,dim);
			}

			fwII = pStruct->pPPData->getFractionalFlux(Sw_I);
			fwJJ = pStruct->pPPData->getFractionalFlux(Sw_J);
			fwIJ = 0.5*(fwII + fwJJ);

			// mid-edge total velocity
			pStruct->pPPData->getVelocity_new(edge,dom,vel);

			// Numerical Flux Function
			const double FluxIJ[3] = { fwIJ*vel[0], fwIJ*vel[1], fwIJ*vel[2] };

			// Fractional Flux Flow Function Derivative (with respect to saturation)
			df_dsIJ = ( fabs(Sw_I-Sw_J) > 1.e-12 )?fabs((fwJJ-fwII)/(Sw_J-Sw_I)) : .0;

			// Approximate Eigenvalue (Note that we are using the linearized form of df_dsIJ)
			n = .0;
			for (int i=0; i<dim; i++){
				n += pow(vel[i],2);
			}
			alpha = sqrt(n)*df_dsIJ;

			// get the maximum alpha to compute the time step
			alpha_max = std::max(alpha,alpha_max);

			// Central difference Contribution
			non_visc_fv = .0;
			for (int i=0; i<dim; i++){
				non_visc_fv +=  FluxIJ[i]*Cij[i];
			}

			// Numerical Diffusion
			non_visc_ad = 0.5*pGCData->getCij_norm(edge,dom)*alpha*(Sw_J - Sw_I);

			#ifdef _SEEKFORBUGS_
			if ( fabs(non_visc_ad) > 0.0 || fabs(non_visc_fv) > 0.0 ) check_IAT = true;
			#endif

			// Computing "Non-Viscous" Terms
			nonvisc_I = pStruct->pPPData->getNonViscTerm(I) + (non_visc_fv - non_visc_ad);///nrc;
			nonvisc_J = pStruct->pPPData->getNonViscTerm(J) - (non_visc_fv - non_visc_ad);///nrc;

			// update nonvisc term. it will be set to 0 at the next time iteration
			pStruct->pPPData->setNonViscTerm(I,nonvisc_I);
			pStruct->pPPData->setNonViscTerm(J,nonvisc_J);
		}
	}
	EIter_delete(eit);
	#ifdef _SEEKFORBUGS_
	   if (!check_IAT){
		   //char msg[256]; sprintf(msg,"For domain %d, non_visc_fv and non_visc_ad are always null!\n",dom);
		   //throw Exception(__LINE__,__FILE__,msg);
	   }
	#endif //_SEEKFORBUGS_

	//alpha_max = P_getMaxDbl(alpha_max);
	courant = pStruct->pSimPar->CFL();
	phi = pStruct->pSimPar->getPorosity(dom);
	length = pGCData->getSmallestEdgeLength();
	
	#ifdef _SEEKFORBUGS_
	if ( fabs(alpha_max) < 1e-8){
		//char msg[256]; sprintf(msg,"alpha_max NULL in domain %d\n",dom);
		//throw Exception(__LINE__,__FILE__,msg);
	}
	#endif //_SEEKFORBUGS_
	
	dt = (courant*length*phi)/alpha_max;
	cout << setprecision(8) << fixed << " ##### TIME STEP = " << dt << endl;

	// use dt as time step to advance saturation
	timeStepByDomain.insert(dt);

	double endt = MPI_Wtime();
	return endt-startt;
}

void EBFV1_hyperbolic::resetNodalNonviscTerms(pMesh theMesh){
	pEntity node;
	alpha_max = .0;
	VIter vit = M_vertexIter(theMesh);
	while ( (node = VIter_next(vit)) ){
		pStruct->pPPData->setNonViscTerm(node,.0);
	}
	VIter_delete(vit);
}
}
