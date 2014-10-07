CXXFLAGS=-g -Wall

APP_DIR=$(HOME)/applications

#MG 1D directory
1D_DIR=/home/julio/Desktop/Programming/libraries/MG_lib/MG_1D

CXX=mpicxx

INCLUDES=-I$(1D_DIR)/include -I$(PETSC_DIR)/include 

1D_OBJ_DIR=$(1D_DIR)/objs
1D_SRC_DIR=$(1D_DIR)/src
1D_OBJECTS=$(1D_OBJ_DIR)/MG_1D.o
1D_SRC=$(1D_SRC_DIR)/MG_1D.cpp

STATICLIB=libmultigrid.a


all:$(STATICLIB)


include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules 

$(STATICLIB):	$(1D_OBJECTS) chkopts
	@echo 'Linking....'
	ar crsu $(STATICLIB) $(1D_OBJECTS)
	@echo "ok. Done!"

$(1D_OBJ_DIR)/%.o:	$(1D_SRC_DIR)/%.cpp 
	$(CXX) $(CXXFLAGS) -c $(INCLUDES)  $<
	@mv *.o $(1D_OBJ_DIR)

cleanup:
	rm *.o $(1D_OBJ_DIR)/*.o *.exe
