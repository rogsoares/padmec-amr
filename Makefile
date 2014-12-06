CXXFLAGS=-g -Wall

APP_DIR=$(HOME)/applications

#MG 1D directory
1D_DIR=/home/julio/Desktop/Programming/libraries/MG_lib/MG_1D

CXX=mpicxx

INCLUDES=-I$(1D_DIR)/include -I$(PETSC_DIR)/include 

OBJ_DIR=$(1D_DIR)/objs
SRC_DIR=$(1D_DIR)/src
OBJECTS_1D=$(OBJ_DIR)/Controller_V.o $(OBJ_DIR)/grid1D.o $(OBJ_DIR)/SetOperator.o
STATICLIB=$(1D_DIR)/lib/libmultigrid.a


all:$(STATICLIB)


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules 

#chkopts

$(STATICLIB):	$(OBJECTS_1D) 
	@echo 'Linking....'
	ar crsu $(STATICLIB) $(OBJECTS_1D)
	@echo "ok. Done!"

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(OBJ_DIR)



cleanup:
	rm  $(OBJ_DIR)/*.o
	rm  $(1D_DIR)/lib/*.a
