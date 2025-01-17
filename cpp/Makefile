#
# This program is licensed granted by the University of Brasília - UnB (the “University”)
# for use of SIMULATION OF MAGNETIC SUSPENSIONS - SIMMSUS (“the Software”) through this website
# https://github.com/rafaelgabler/simmsus (the ”Website”).
#
# By downloading the Software through the Website, you (the “Licensee”) are confirming that you agree
# that your use of the Software is subject to the academic license terms.
#
# For more information about SIMMSUS please contact: rafael.gabler@unb.br (Rafael Gabler Gontijo)
# or leandro.zanotto@gmail.com (Leandro Zanotto).
#
#
#  For Debug CCFLAGS= -g -O3 -std=c++0x -fopenmp -march=native -fma -mtune=native -ffast-math -fopt-info-vec-optimized
SRC_DIR=src
CCFLAGS= -O0
INC_DIR = .
CC=g++
LDFLAGS=
OBJ= config.o simmsus.o globals.o boxSize.o particleDistribution.o randomic.o greenTable.o periodicStructure.o

all: simmsus

simmsus: ${OBJ}
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -o $@ ${OBJ}

config.o: $(SRC_DIR)/config.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

globals.o: $(SRC_DIR)/globals.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

boxSize.o: $(SRC_DIR)/boxSize.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

particleDistribution.o: $(SRC_DIR)/particleDistribution.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

randomic.o: $(SRC_DIR)/randomic.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<		

greenTable.o: $(SRC_DIR)/greenTable.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

periodicStructure.o: $(SRC_DIR)/periodicStructure.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<				

simmsus.o: $(SRC_DIR)/simmsus.cpp
	$(CC) $(CCFLAGS) -I$(INC_DIR) $(LDFLAGS) -c $<

clean:
	rm *.o simmsus