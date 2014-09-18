# ======================================================
# 			UNIVERSIDADE FEDERAL DE PERNANMBUCO
#           DEPARTAMENTO DE ENGENHARIA MECANICA (DEMEC-CTG-UFPE)
#			NUCLEO DE TECNOLOGIA (NT-CAA-UFPE)
#
# Authors: 	Rogerio Soares(rogerio.soaress@ufpe.br)
#			Saulo
#			Guilherme Caminha
#			Paulo Lyra
# Created:	2008-2013
# ======================================================

# FLAGS LEGENDS:
# 	FVPO                   : For Visualization Purposes Only
# 	-g:                    : Debugging compilation (Use it to debug program by DDD)
# 	__ADAPTATION_DEBUG__   : Debug Adaptation functions only
#	__ERROR_ANALYSIS_DEBUG__: Debug error analysis functions only
#	_SEEKFORBUGS_           : Debuggin in general (whole code)
#    NOADAPTATION			: Compile code without mesh adaptation feature
#
CXXFLAGS=-DPARALLEL -g -DNOADAPTATION -Wunused-local-typedefs -D_SEEKFORBUGS_ -DTRACKING_PROGRAM_STEPS
#CXXFLAGS=-DPARALLEL -g -Wall -Wunused -D__ADAPTATION_DEBUG__ -D__ERROR_ANALYSIS_DEBUG__ -D_SEEKFORBUGS_ -DTRACKING_PROGRAM_STEPS -DFVPO
# scientific_application folder is where libraries are installed 
APP_DIR=$(HOME)/scientific_applications
PROJ_DIR=$(HOME)/projetos/GOOGLE_CODES/padmec-amr/trunk
GMSH_DIR=$(APP_DIR)/gmshGMSH

# compilers
CXX=mpicxx

# paths for headers
INCLUDES=-I$(PROJ_DIR)/include \
	-I$(PROJ_DIR)/src/adaptation -I$(PROJ_DIR)/src/adaptation/adaptive-remeshing -I$(PROJ_DIR)/src/adaptation/h-refinement \
	-I$(PROJ_DIR)/src/adaptation/rh-refinement -I$(APP_DIR)/MAdLib_no_gmsh/include \
	-I$(PROJ_DIR)/src/elliptic -I$(PROJ_DIR)/src/hyperbolic \
	-I$(PROJ_DIR)/src/error-analysis -I$(PROJ_DIR)/src/interpolation -I$(PROJ_DIR)/src/pre-processor \
	-I$(PROJ_DIR)/src/SIMULATION_core -I$(PROJ_DIR)/src/tools \
	-I$(GMSH_DIR)/build/Common -I$(GMSH_DIR)/Common -I$(GMSH_DIR)/Geo -I$(GMSH_DIR)/Mesh -I$(GMSH_DIR)/Numeric -I$(GMSH_DIR)/Parser -I$(GMSH_DIR)/Plugin -I$(GMSH_DIR)/Post \
 	-I$(APP_DIR)/autopack/include -I$(APP_DIR)/FMDB-2011/include

LIBS1=-L$(APP_DIR)/FMDB-2011/lib -lFMDB-O \
     -L$(APP_DIR)/ParMetis-3.1/lib -lparmetis -lmetis \
     -L$(APP_DIR)/autopack/lib -lautopack-O

ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS)))
LIBS2=-L$(APP_DIR)/gmshGMSH/lib -lGmsh
endif

LIBS=$(LIBS1) $(LIBS2)

# path for where the objects must be put
# ======================================
OBJ_DIR=$(PROJ_DIR)/objs

# paths for each program source part
# ==================================
SRC_DIR1=$(PROJ_DIR)/src
SRC_DIR2=$(PROJ_DIR)/src/elliptic
SRC_DIR4=$(PROJ_DIR)/src/hyperbolic
SRC_DIR5=$(PROJ_DIR)/src/pre-processor
SRC_DIR6=$(PROJ_DIR)/src/SIMULATION_core
SRC_DIR7=$(PROJ_DIR)/src/tools
SRC_DIR8=$(PROJ_DIR)/src/error-analysis
SRC_DIR9=$(PROJ_DIR)/src/adaptation
SRC_DIR10=$(PROJ_DIR)/src/interpolation
SRC_DIR11=$(PROJ_DIR)/src/adaptation/adaptive-remeshing
SRC_DIR12=$(PROJ_DIR)/src/adaptation/h-refinement
#SRC_DIR13=$(PROJ_DIR)/src/adaptation/rh-refinement

OBJS1=$(OBJ_DIR)/main.o $(OBJ_DIR)/SIMULATION_core.o $(OBJ_DIR)/SIMULATION_core__solvers.o $(OBJ_DIR)/SIMULATION_adaptation.o \
	$(OBJ_DIR)/EBFV1_elliptic.o $(OBJ_DIR)/set_SOE.o $(OBJ_DIR)/EBFV1_E.o $(OBJ_DIR)/EBFV1_G.o $(OBJ_DIR)/EBFV1_F_omega.o $(OBJ_DIR)/EBFV1_F_gamma.o \
	$(OBJ_DIR)/EBFV1__DefectCorrectionSolver.o $(OBJ_DIR)/EBFV1__MatrixFreeSolver.o $(OBJ_DIR)/EBFV1_pressure_gradient.o $(OBJ_DIR)/EBFV1_velocityField.o \
	$(OBJ_DIR)/EBFV1_hyperbolic.o $(OBJ_DIR)/EBFV1_hyperbolic-adaptative.o $(OBJ_DIR)/EBFV1_advectiveTerm.o $(OBJ_DIR)/EBFV1_advanceSaturation.o \
	$(OBJ_DIR)/SaturationGradient.o \
	$(OBJ_DIR)/MeshData.o $(OBJ_DIR)/OilProductionManagement.o \
	$(OBJ_DIR)/EBFV1_preprocessor.o $(OBJ_DIR)/Calculate-Cij-parallel.o $(OBJ_DIR)/Calculate-Vi-parallel.o $(OBJ_DIR)/EBFV1__pre-processors.o \
	$(OBJ_DIR)/EBFV1-2D-pp.o $(OBJ_DIR)/EBFV1-3D-pp.o $(OBJ_DIR)/setCorrectNumberOfRemoteCopies.o $(OBJ_DIR)/validate-EBFV1.o \
    $(OBJ_DIR)/LoadMeshAdaptationParameters.o  $(OBJ_DIR)/load_EBFV1_PreProcessorData.o $(OBJ_DIR)/auxiliar.o $(OBJ_DIR)/SimulatorParameters.o \
    $(OBJ_DIR)/Exception.o $(OBJ_DIR)/exportVTK.o $(OBJ_DIR)/PhysicPropData.o $(OBJ_DIR)/LoadSimulatorParameters.o \
    $(OBJ_DIR)/GeomData.o $(OBJ_DIR)/GeomData_initialize.o $(OBJ_DIR)/GeomData_calculate.o $(OBJ_DIR)/GeomData_Mapping.o $(OBJ_DIR)/GeomData_AllocateMemory.o \
    $(OBJ_DIR)/GeomData_TransferData.o
	
ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS)))
OBJS2=$(OBJ_DIR)/CalculateDegreeOfRefinement_2D.o $(OBJ_DIR)/AdaptiveRemeshing.o \
    $(OBJ_DIR)/H_Refinement.o $(OBJ_DIR)/H_Refinement_2D.o \
    $(OBJ_DIR)/CalculateElementsError_2D.o $(OBJ_DIR)/CalculateGlobalError.o \
    $(OBJ_DIR)/ErrorAnalysis.o $(OBJ_DIR)/ErrorAnalysisAuxiliar.o $(OBJ_DIR)/CalculateSmoothGradientNorm_2D.o \
    $(OBJ_DIR)/MeshRegularization.o $(OBJ_DIR)/SpecialElements_2D.o \
    $(OBJ_DIR)/UnevenElements.o $(OBJ_DIR)/IhR_main.o $(OBJ_DIR)/interpolation.o \
	$(OBJ_DIR)/IAR_calculate.o $(OBJ_DIR)/IAR_gradients.o $(OBJ_DIR)/FileManager.o $(OBJ_DIR)/FileManagerFunctions.o \
	$(OBJ_DIR)/MRE.o $(OBJ_DIR)/Rebuilder2.o $(OBJ_DIR)/ErrorEstimator.o $(OBJ_DIR)/Rebuilder3.o
endif	
OBJECTS=$(OBJS1) $(OBJS2) 

# how the executable will be named
EXEC=PADMEC_AMR.exe
all:	$(EXEC)
include ${PETSC_DIR}/bmake/common/base

$(EXEC):	$(OBJECTS) chkopts
	@echo "Linking objects..."
	-$(CXX) -o $(EXEC) $(OBJECTS) $(LIBS) $(PETSC_LIB)	

$(OBJ_DIR)/%.o:	$(SRC_DIR1)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR2)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR4)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR5)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR6)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR7)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)

ifeq (,$(findstring NOADAPTATION,$(CXXFLAGS)))
$(OBJ_DIR)/%.o:	$(SRC_DIR8)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
		
$(OBJ_DIR)/%.o:	$(SRC_DIR9)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
	
$(OBJ_DIR)/%.o:	$(SRC_DIR10)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR11)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)

$(OBJ_DIR)/%.o:	$(SRC_DIR12)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
	@mv *.o $(OBJ_DIR)
endif

#$(OBJ_DIR)/%.o:	$(SRC_DIR13)/%.cpp
#	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $< ${PETSC_INCLUDE}
#	@mv *.o $(OBJ_DIR)
	
rebuild:
	@rm -f $(OBJ_DIR)/*.o $(EXEC)
	@echo "Limpeza concluida"
	@make
