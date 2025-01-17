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
#include <vector>
#include <math.h>
#include <header/config.hpp>
#include <header/constants.hpp>
#include <header/globals.hpp>
#include <header/particleDistribution.hpp>
#include <header/boxSize.hpp>
#include <header/greenTable.hpp>
#include <header/periodicStructure.hpp>

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
double *cof01 = NULL;
double *cof11 = NULL;
double *cof02 = NULL;
double *cof12 = NULL;
double *cof03 = NULL;
double *cof13 = NULL;
double *cof04 = NULL;
double *cof14 = NULL;
double *cof05 = NULL;
double *cof15 = NULL;
double *cof06 = NULL;
double *cof16 = NULL;
double *cof07 = NULL;
double *cof17 = NULL;
double *cof08 = NULL;
double *cof18 = NULL;
const int nGreen = 10000;

//Globals Initializations
X0 = new double[numRealizations * numParticles];
X1 = new double[numRealizations * numParticles];
X2 = new double[numRealizations * numParticles];
U0 = new double[numRealizations * numParticles];
U1 = new double[numRealizations * numParticles];
U2 = new double[numRealizations * numParticles];
XI0 = new double[nb * numRealizations * numParticles];
XI1 = new double[nb * numRealizations * numParticles];
XI2 = new double[nb * numRealizations * numParticles];
ILF0 = new double[125];
ILF1 = new double[125];
ILF2 = new double[125];
ILR0 = new double[9];
ILR1 = new double[9];
ILR2 = new double[9];


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
particleDistribution(configuration.getMonopolidisp(), diam, betaVec);

// Calculating the size of the simulation box based on the 
// number of particles and on the volume fraction defined
// by the user in the simconfig.dat

boxSize(numParticles, configuration.getVolumefracpart(), configuration.getBoxaspectratio(), configuration.getInitialspheraggr());
// //initial condition

double qsi = pow(M_PI,0.5)  / (pow(l * l * h,(1.0/3.0)));

// Building a table with all the Green-functions
// And building the periodic structure to compute
// the Ewald summations

if(periodicity){
    double wave = pow(nbr,(1.0/3.0));
    if(wave == 5.0){
        if(configuration.getAccounthi()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            greenTableLigaihWave5(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13);
        }else if(configuration.getPmt()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof05 = new double[nGreen];
            cof15 = new double[nGreen];
            cof07 = new double[nGreen];
            cof17 = new double[nGreen];
            greenTableTmagper5(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13,
            cof04, cof14, cof05, cof15, cof07, cof17);
        }else if(configuration.getPmf()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof06 = new double[nGreen];
            cof16 = new double[nGreen];
            cof08 = new double[nGreen];
            cof18 = new double[nGreen];
            greenTableFmagper3(nGreen, qsi, l, nb, wave, h, cof01,  cof11, cof02, cof12, cof03, cof13,
            cof06, cof16, cof08, cof18);
        }
    }else if(wave == 3.0){
        if(configuration.getAccounthi()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            greenTableLigaihWave3(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13);
        }else if(configuration.getPmt()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof05 = new double[nGreen];
            cof15 = new double[nGreen];
            cof07 = new double[nGreen];
            cof17 = new double[nGreen];
            greenTableTmagper3(nGreen, qsi, l, nb, wave, h, cof01,  cof11, cof02, cof12, cof03, cof13,
            cof04, cof14, cof05, cof15, cof07, cof17);
        }else if(configuration.getPmf()){
            cof01 = new double[nGreen];
            cof11 = new double[nGreen];
            cof02 = new double[nGreen];
            cof12 = new double[nGreen];
            cof03 = new double[nGreen];
            cof13 = new double[nGreen];
            cof04 = new double[nGreen];
            cof14 = new double[nGreen];
            cof06 = new double[nGreen];
            cof16 = new double[nGreen];
            cof08 = new double[nGreen];
            cof18 = new double[nGreen];
            greenTableFmagper3(nGreen, qsi, l, nb, wave, h, cof01, cof11, cof02, cof12, cof03, cof13,
            cof06, cof16, cof08, cof18);            

        }
    }
}
periodicStructure();
// for(int j = 0; j < numRealizations; j++){
//     for(int i = 0; i < numParticles; i++){
//         cout << X0[j * numRealizations + i] << "\t" << X1[j * numRealizations + i] << "\t" << X2[j * numRealizations + i];
//         cout << U0[j * numRealizations + i] << "\t" << U1[j * numRealizations + i] << "\t" << U2[j * numRealizations + i];
//     }
// }

cout << "END\n";


//free memory

delete[] X0;
delete[] X1;
delete[] X2;
delete[] U0;
delete[] U1;
delete[] U2;
delete[] XI0;
delete[] XI1;
delete[] XI2;
delete[] ILF0;
delete[] ILF1;
delete[] ILF2;
delete[] ILR0;
delete[] ILR1;
delete[] ILR2;

if(cof01 != NULL)
    delete[] cof01;
if(cof11 != NULL)
    delete[] cof11;
if(cof02 != NULL)
    delete[] cof02;
if(cof12 != NULL)
    delete[] cof12;
if(cof03 != NULL)
    delete[] cof03;
if(cof13 != NULL)
    delete[] cof13;
if(cof04 != NULL)
    delete[] cof04;
if(cof14 != NULL)
    delete[] cof14;
if(cof05 != NULL)
    delete[] cof05;
if(cof15 != NULL)
    delete[] cof15;
if(cof06 != NULL)
    delete[] cof06;
if(cof16 != NULL)
    delete[] cof16;
if(cof07 != NULL)
    delete[] cof07;
if(cof17 != NULL)
    delete[] cof17;
if(cof08 != NULL)
    delete[] cof08;
if(cof18 != NULL)
    delete[] cof18;

}
return 0;
}