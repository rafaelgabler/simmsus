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
#include <header/globals.hpp>

int nnr;
int npast;
int numRealizations;
int numParticles;
int totalRealParticle;
double l, h;
int numBoxes, numRecBoxes;
double shearratei;
double per, phi;
bool periodicity;
double *diam;
double *betaVec;
double *X0, *X1, *X2;
double *U0, *U1, *U2;
double *XI0, *XI1, *XI2;
double *ILF0, *ILF1, *ILF2;
double *ILR0, *ILR1, *ILR2;