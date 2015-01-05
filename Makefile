# ======================================================
# 			UNIVERSIDADE FEDERAL DE PERNANMBUCO
#           DEPARTAMENTO DE ENGENHARIA MECANICA (DEMEC-CTG-UFPE)
#			NUCLEO DE TECNOLOGIA (NT-CAA-UFPE)
#
# Authors: 	Rogerio Soares(rogerio.soaress@ufpe.br)
# Created:	2015
# ======================================================


CXXFLAGS=-g -Wall
PROJ_DIR=$(HOME)/projetos/mesh-data
CXX=mpicxx
INCLUDES=-I$(PROJ_DIR)
SOURCES=$(PROJ_DIR)
OBJ=main.o read.o write.o mesh.o refine.o refine_TRI.o refine_QUAD.o refine_TETRA.o utilities.o
EXEC=mesh-data.exe

all:	$(EXEC)

$(EXEC):	$(OBJ)
	@echo "Linking objects..."
	-$(CXX) -o $(EXEC) $(OBJ)	

%.o:	$(SOURCES)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES)  -c $<	
	
rebuild:
	@echo "Cleaning all object file..."
	@rm -f *.o $(EXEC)
	@echo "OK. Now compiling all source files..."
	@make
	@echo "Done!"
