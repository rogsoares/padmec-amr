
# ======================================================
# 		UNIVERSIDADE FEDERAL DE PERNANMBUCO
#		NUCLEO DE TECNOLOGIA (NT-CAA-UFPE)
#Author: 	Julio Cezar
#Contact:	julio-cezar08@hotmail.com
# ======================================================


CXXFLAGS=-g -Wall

APP_DIR=$(HOME)/applications

#MG 1D directory
1D_DIR=/home/julio/Desktop/Programming/libraries/MG_lib/MG_1D

CXX=mpicxx

INCLUDES=-I$(1D_DIR)/include -I$(PETSC_DIR)/include 

OBJ_DIR=$(1D_DIR)/objs
SRC_DIR=$(1D_DIR)/src

OBJECTS_VCYCLE=$(OBJ_DIR)/Controller_V.o 
OBJECTS_GRID1D=$(OBJ_DIR)/grid1D.o $(OBJ_DIR)/SetOperator.o $(OBJ_DIR)/gridGenerator.o 
OBJECTS_SOLVER=$(OBJ_DIR)/smoother.o
OBJECTS_TOOLS=$(OBJ_DIR)/residualCalculator.o
OBJECTS= $(OBJECTS_VCYCLE) $(OBJECTS_GRID1D) $(OBJECTS_SOLVER) $(OBJECTS_TOOLS)

STATICLIB=$(1D_DIR)/lib/libmultigrid.a

all:$(STATICLIB)

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules 

#chkopts

$(STATICLIB):	$(OBJECTS)
	@echo 'Linking....'
	ar crsu $(STATICLIB) $(OBJECTS) 
	@echo "ok. Done!"

#V-cycle 
$(OBJ_DIR)/%.o:	$(SRC_DIR)/V-cycle/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(OBJ_DIR)

#grid's components and operators
$(OBJ_DIR)/%.o: $(SRC_DIR)/grid/1D/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(OBJ_DIR)

#solver
$(OBJ_DIR)/%.o: $(SRC_DIR)/smoother/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(OBJ_DIR)

#tools
$(OBJ_DIR)/%.o: $(SRC_DIR)/tools/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(OBJ_DIR)



cleanup:
	rm  $(OBJ_DIR)/*.o
	rm  $(1D_DIR)/lib/*.a
