/*
 * EBFV1-hyperbolic.cpp
 *
 *  Created on: 09/01/2009
 *      Author: rogerio
 */

#include "EBFV1_hyperbolic_adaptative.h"
#include "exportVTK.h"

namespace PRS{

	EBFV1_hyperbolic_adaptative::EBFV1_hyperbolic_adaptative(){
	}

	EBFV1_hyperbolic_adaptative::EBFV1_hyperbolic_adaptative(pMesh mesh,
			PhysicPropData *ppd, SimulatorParameters *sp,
			GeomData *gcd,MeshData *md, OilProductionManagement *popm, ErrorAnalysis *pEA){

		pMData = md;
		pOPManager = popm;

		/*
		 * Experimental:
		 */
		pStruct = new PointerStruct;
		pStruct->pPPData = ppd;
		pGCData = gcd;
		pStruct->pSimPar = sp;
		pStruct->theMesh = mesh;
		pStruct->pErrorAnalysis = pEA;
		pHOApproximation = new HighOrderApproximation(pStruct);

		initial_vel_ts_old = true;
		print_ATE_Data = true;
		vel_ts_counter = 0;
		sat_ts_counter = 0;
		DVTOL = pStruct->pSimPar->getDVTOL();
		createDimensionLessFactors();
		//rowToImport = 0;
	}

	EBFV1_hyperbolic_adaptative::~EBFV1_hyperbolic_adaptative(){
	}

	double EBFV1_hyperbolic_adaptative::solver(double &timeStep){
		if (!P_pid()) std::cout << "Adaptative Hyperbolic solver...\n";

		// reset logical variables which control programming flux.
		resetVariables();

		double hyp_time = .0;		// CPU time usage counter for hyperbolic equation
		while ( allowAdaptativeTimeAdvance ){

			// Set flux to zero before evaluate saturation equation.
			resetNodalNonviscTerms();

			// calculate saturation gradient
			if (pStruct->pSimPar->useHOApproximation()) hyp_time += calculateSaturationGradient();

			/*
			 * Loop over domains:
			 * Avoid physical inconsistent between domains boundary.
			 */
			int dom_counter = 0;
			SIter dom=pStruct->pSimPar->setDomain_begin();
			for (; dom!=pStruct->pSimPar->setDomain_end();dom++){

				if (allowVelocityCalculation) hyp_time += calculateVelocityField(*dom,dom_counter);

				// If a high order approximation for saturation equation is required...
				if (pStruct->pSimPar->useHOApproximation()){
					//NodeSlopeLimiter* pNodeSL = pHOApproximation->getNodeSL_Ptr();
					hyp_time += pHOApproximation->getNodeSL_Ptr()->defineSlopeLimiters();
				}
				hyp_time += calculateIntegralAdvectiveTerm(*dom);
			}
			pMData->unifyScalarsOnMeshNodes(PhysicPropData::getNonViscTerm,PhysicPropData::setNonViscTerm,pStruct->pGCData,0);

			// calculate explicit time-step (CFL constrained)
			sat_ts = getTimeStep();
			if (!P_pid()) printf("Time-step: %f\n",sat_ts);

			// correct time-step value to print out the desired simulation moment
			pStruct->pSimPar->correctTimeStep(sat_ts);

			// calculate implicit time-step
			vel_ts = calculateNewImplicitTS();

			/*
			 * Compute saturation field for n+1 (explicit saturation).
			 * If saturation is bigger that 1.0 an exception will be thrown.
			 */
			try{
				hyp_time += calculateExplicitAdvanceInTime(sat_ts);
			}
			catch (Exception excp){
				excp.showExceptionMessage();
				pStruct->pSimPar->stopSimulation();
				allowAdaptativeTimeAdvance = false;
			}
			sat_ts_counter++;

			// Output data (VTK)
			pStruct->pSimPar->printOutVTK(pStruct->theMesh,pStruct->pPPData,pStruct->pErrorAnalysis,pStruct->pSimPar,exportSolutionToVTK);

			/*
			 * The following lines below will be condensed to a function member call
			 * and it will belong to EBFV1_hyperbolic
			 */
			if (pStruct->pSimPar->rankHasProductionWell()){
				pOPManager->printOilProduction(sat_ts,
						pStruct->pSimPar->getAccumulatedSimulationTime(),
						pStruct->pSimPar->getSimTime(),
						getRecoveredOilValue());
			}
		}
		vel_ts_counter++;
		timeStep = vel_ts;

		if (!P_pid()){
			printAdatptativeTimeEvaluationData();
			std::cout << "done.\n------------------------------------------------------------------\n\n";
		}
		return hyp_time;
	}

	/*
	 * Compute a new implicit time-step (vel_ts) if necessary. This happens
	 * either at the beginning of simulation or when sum of sat_ts is greater
	 * than DT.
	 */
	double EBFV1_hyperbolic_adaptative::calculateNewImplicitTS(){
		if (allowVelocityCalculation){
			/*
			 *  stop velocity field be evaluated during saturation advance.
			 */
			allowVelocityCalculation = false;

			/*
			 *  set vel_ts_old as sat_ts just once at the beginning of simulation.
			 */
			if (initial_vel_ts_old){
				vel_ts_old         = sat_ts;
				initial_vel_ts_old = false;
			}

			// make implicit vel_ts_old variable dimensionless
			vel_ts_old *= time_factor;

			// implicitTS returns a dimensionless delta T
			dv_norm = calculateVelocityVariationNorm();

			// Update new implicit time step: vel_ts
			vel_ts = (DVTOL/dv_norm)*vel_ts_old;
			if (!P_pid()) printf("DVTOL: %f dv_norm: %f  DVTOL/dv_norm: %f time_factor: %f\tvel_ts_old: %f\tvel_ts: %f\n",
					DVTOL,dv_norm,DVTOL/dv_norm,time_factor,vel_ts_old,vel_ts);

			// dimensionless DT ratio
			double RT = (double)(vel_ts/vel_ts_old);
			if (RT>1.25)
				vel_ts = 1.25*vel_ts_old;
			else if ( RT<0.75 )
				vel_ts = 0.75*vel_ts_old;

			// gives implicit vel_ts variable dimension of time
			vel_ts /= time_factor;
//			printf("DVTOL: %f dv_norm: %f  DVTOL/dv_norm: %f  time_factor: %f\tvel_ts_old: %f\tvel_ts: %f\n",
//								DVTOL,dv_norm,DVTOL/dv_norm,time_factor,vel_ts_old,vel_ts);
			vel_ts_old = vel_ts;

			// It doesn't make sense vel_ts be less than sat_ts
			if (vel_ts < sat_ts){
				vel_ts = sat_ts;
				allowAdaptativeTimeAdvance = false;
				//if (!P_pid()) std::cerr << "\tWARNNING: vel_ts("<< vel_ts << ") < " << "sat_ts("<<sat_ts<<")." << std::endl;
			}
		}
		/*
		 * Verify explicit advance compared to implicit advance:
		 * sum(sat_ts)<=vel_ts
		 */
		cumulativeExplicitTS += sat_ts;				// summation of all sat_ts within the implicit time-step (vel_ts)
		verifyExplicitAdvanceForImplicitTS();
		//if (!P_pid()) printf("TS=>  implicit: %f  explicit: %f(sat_ts) %f(accumulated)\n",vel_ts,sat_ts,cumulativeExplicitTS);

		/*
		 * Verify explicit advance compared to simulation time.
		 * AccSimTime = AccSimTime + timeStep
		 * Summation of all sat_ts for the whole simulation
		 */
		pStruct->pSimPar->setAccumulatedSimulationTime(sat_ts);
		verifyExplicitAdvanceForSimulationTime();

		return vel_ts;
	}

	/*
	 * Compute a new implicit time-step as function of the last implicit time-step,
	 * velocity norm and DVTOL. DVTOL comes from the input data file 'numeric.dat'.
	 */
	double EBFV1_hyperbolic_adaptative::calculateVelocityVariationNorm(){
		const int dim = pGCData->getMeshDim();
		const int ndom = pStruct->pSimPar->getNumDomains();
		const int numGEdges = pGCData->getNumGEdges();
		SIter domIterator = pStruct->pSimPar->setDomain_begin();

		pEntity edge;
		dblarray v_new(dim), v_old(dim), Cij(dim), Dv_dom(ndom,.0);
		dblarray v_var(dim); // velocity variation

		int i, k = 0;
		//double sumDV;

		// loop over domains
		for (SIter_const dom=pStruct->pSimPar->setDomain_begin(); dom!=pStruct->pSimPar->setDomain_end();dom++){
			Dv_dom[k] = .0;
			// loop over domains edges

			EIter eit = M_edgeIter(pStruct->theMesh);
			while ( (edge = EIter_next(eit)) ){
			// we are supposing only one domain for while
				if ( pGCData->edgeBelongToDomain(edge,*dom) ){
					pGCData->getCij(edge,*dom,Cij);
					pStruct->pPPData->getVelocity_new(edge,*dom,v_new);
					pStruct->pPPData->getVelocity_old(edge,*dom,v_old);

					// make dimensionless velocity
					for (i=0; i<dim; i++){
						v_new[i] *= vel_factor;
						v_old[i] *= vel_factor;
						v_var[i] = v_new[i]-v_old[i];
					}
					// get Cij norm
					double norm = pGCData->getCij_norm(edge,*dom);

					// number of remote copies
					//double nrc = (double)pGCData->getNumRC(edge,*dom) + 1.0;
					double nrc = (double)pGCData->getNumRC(pStruct->theMesh,edge) + 1.0;
					cout << nrc << endl;


					double inner_prod = .0;
					for (i=0; i<dim; i++) inner_prod += (v_new[i]-v_old[i])*Cij[i];
					Dv_dom[k] += pow(fabs(inner_prod)/norm,2)/nrc;

					// pre- Flux variation norm
					//Dv_dom[k] += pow(fabs(inner_product(v_var,Cij))/norm,2)/nrc;
				}
			}
			EIter_delete(eit);


			Dv_dom[k] = P_getSumDbl(Dv_dom[k]);
			Dv_dom[k] /= (double)numGEdges;
			k++;
		}

		// take an average of velocity norm for all domains
		double vnorm = sqrt( (double)(std::accumulate(Dv_dom.begin(),Dv_dom.end(),.0)/ndom) );
		printf("VNORM: %f\n",vnorm);
		return vnorm;
	}

	void EBFV1_hyperbolic_adaptative::resetVariables(){
		cumulativeExplicitTS = .0;
		allowVelocityCalculation = true;
		allowAdaptativeTimeAdvance = true;
	}

	void EBFV1_hyperbolic_adaptative::verifyExplicitAdvanceForImplicitTS(){
		if ( cumulativeExplicitTS >= vel_ts ){
			sat_ts = sat_ts - (cumulativeExplicitTS - vel_ts);
			allowAdaptativeTimeAdvance = false;
		}
	}

	void EBFV1_hyperbolic_adaptative::verifyExplicitAdvanceForSimulationTime(){
		double accSimTime = pStruct->pSimPar->getAccumulatedSimulationTime();
		double simTime = pStruct->pSimPar->getSimTime();
		if ( accSimTime >= simTime ){
			sat_ts = sat_ts - (accSimTime - simTime);	// correct sat_ts (last step)
			allowAdaptativeTimeAdvance = false;
			pStruct->pSimPar->stopSimulation();
		}
	}

	void EBFV1_hyperbolic_adaptative::printAdatptativeTimeEvaluationData(){

		if (print_ATE_Data){

			int size = P_size();									// number fo processes
			double dvtol = pStruct->pSimPar->getDVTOL();			// dvtol value
			string path = pStruct->pSimPar->getOutputPathName();	// where file will be placed

			char fname[256]; sprintf(fname,"%s_AdaptTimeData-DVTOL__%.2E__%dp.xls",path.c_str(),dvtol,size);
			fid.open(fname);

			if ( !fid.is_open() ){
				char msg[256]; sprintf(msg,"%s could not be opened or it does not exist.\n",fname);
				throw Exception(__LINE__,__FILE__,msg);
			}

			fid << "Time-step ImpExp_TimeStep-ratio velocity_norm\n";
			print_ATE_Data = false;
		}

		if (!P_pid()) {
			double accSimTime = pStruct->pSimPar->getAccumulatedSimulationTime();
			double simTime = pStruct->pSimPar->getSimTime();
			double time_ratio = (double)sat_ts_counter/vel_ts_counter;
			// convert numeric data into a string
			char lineStr[256];
			sprintf(lineStr,"%f %f %f",accSimTime/simTime,time_ratio,dv_norm);

			cout << "\n\n\n" << lineStr << endl;



			string theString(lineStr);
			replaceAllOccurencesOnString(theString,1,".",",");

			cout << "\n\n\n" << theString << endl;


			fid <<  theString << std::endl;
		}
	}

	void EBFV1_hyperbolic_adaptative::createDimensionLessFactors(){
		double L = .0;
		double H = .0;
		pStruct->pSimPar->getReservoirGeometricDimensions(L,H);
		double Q = pStruct->pSimPar->getTotalInjectionFlowRate();
		time_factor = Q/(pStruct->pSimPar->getPorosity(3300)*pow(L,2)*H);
		vel_factor = L*H/Q;
		//printf("vel_factor = %f\n",vel_factor); throw 1;
	}
}
