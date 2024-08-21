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

#include <math.h>
#include <header/globals.hpp>
#include <header/randomic.hpp>
void particleDistribution(bool poliDispersidade, int numRealizations, int numParticles, double *diam, double *beta){


if(poliDispersidade){
    // Allocating variables in the memory
    double *diarand = new double[totalRealParticle];
    randomic(1.5,2.5, (numParticles * numRealizations), diarand);

    for(int j = 0; j < numRealizations; j++){
        for(int i = 0; i < numParticles; i++){
            diam[j * numParticles + i] = diarand[j * numParticles + i];
            beta[j * numParticles + i] = diam[j * numParticles + i] / diam[j * numParticles + 1];
        }
    }
    delete[] diarand;
}
else {
    for(int j = 0; j < numRealizations; j++){
        for(int i = 0; i < numParticles; i++){
            diam[j * numParticles + i] = 2.0;
            beta[j * numParticles + i] = diam[j * numParticles + i] / diam[j * numParticles + 1];
        }
    }
}
}