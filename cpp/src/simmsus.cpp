/*
 * This program is licensed granted by the University of Brasília - UnB (the “University”)
 * for use of SIMULATION OF MAGNETIC SUSPENSIONS - SIMMSUS (“the Software”) through this website
 * https://github.com/rafaelgabler/simmsus (the ”Website”).
 *
 * By downloading the Software through the Website, you (the “Licensee”) are confirming that you agree
 * that your use of the Software is subject to the academic license terms.
 *
 * For more information about SIMMSUS please contact: rafael.gabler@unb.br (Rafael Gabler Gontijo)
 * or leandro.zanotto@gmail.com (Leandro Zanotto).
 *
 */
#include <iostream>
#include <math.h>
#include <header/config.hpp>
#include <header/constants.hpp>
#include <header/globals.hpp>
#include <header/particleDistribution.hpp>
#include <header/boxSize.hpp>
// #include <greenTable.hpp>
// #include <brownian.hpp>
// #include <repulsion.hpp>
// #include <periodicStructure.hpp>
using namespace std;

int main(int argc, char **argv){


if (argc < 2){
    cout << "File Not Found - Please specify an input configuration file as paramenter\n";
}else{
//Read the Input File and Create the Object

Configuration configuration(argv[1]);

std::cout << configuration.getSte() << "\n";

shearratei = configuration.getDynincrshrate(); //! shear-rate -> Arquivo de configuracao
per = (4.0 / 3.0) * configuration.getBrownianpecletnum(); //! rotational Peclet number -> Arquivo de configuracao

phi = configuration.getVolumefracpart();
numParticles = configuration.getNumpart();
numRealizations = configuration.getNumreal();
totalRealParticle = numParticles * numRealizations;

if(configuration.getStatanalysis()){
    npast = 2;
}else{
    npast = configuration.getSimulationtime()  / configuration.getNumtimestep();
}

// Number of random numbers used in each time-step
nnr = 3 * totalRealParticle;



if(configuration.getAccounthi() || configuration.getPmf() || configuration.getPmt()){
    periodicity = true;
}

diam = new double[totalRealParticle]{};
betaVec = new double[totalRealParticle]{};

// Defining the size of the particles
particleDistribution(configuration.getMonopolidisp(), numRealizations, numParticles, diam, betaVec);

// Calculating the size of the simulation box based on the 
// number of particles and on the volume fraction defined
// by the user in the simconfig.dat

boxSize(numParticles, configuration.getVolumefracpart(), configuration.getBoxaspectratio(), configuration.getInitialspheraggr());

}

return 0;
}