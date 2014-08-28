CXXFLAGS=-g -Wall

APP_DIR=$(HOME)/applications
PROJ_DIR=/home/julio/Desktop/Programming/prototypes/SIMULATORPROJ1

CXX=mpicxx

INCLUDES=-I$(PROJ_DIR)/include -I$(PETSC_DIR)/include

OBJ_DIR=$(PROJ_DIR)/objs
SRC_DIR=$(PROJ_DIR)/src
OBJECTS=$(OBJ_DIR)/main.o
SRC=$(SRC_DIR)/main.cpp

EXEC=simulator.exe

all:	$(EXEC)

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules 

$(EXEC):	$(OBJECTS) chkopts
	@echo 'Linking....'
	$(CXX) -o $(EXEC) $(OBJECTS) $(PETSC_KSP_LIB) $(LIBS)
	@echo "ok. Done!"

$(OBJ_DIR)/%.o:	$(SRC_DIR)/%.cpp
	$(CXX) -c $(INCLUDES) $<
	@mv *.o $(OBJ_DIR)

cleanup:
	rm *.o $(OBJ_DIR)/*.o *.exe








