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

#ifndef SRC_HEADERS_GLOBALS_HPP_
#define SRC_HEADERS_GLOBALS_HPP_
#include <math.h>

extern int nnr;
extern int npast;
extern int numRealizations;
extern int numParticles;
extern int num;
extern double l;
extern double h;
extern int totalRealParticle;
extern int numBoxes;
extern int numRecBoxes;
extern double shearratei; //= shearrate; //! shear-rate -> Arquivo de configuracao
extern double per; // = (4.0 / 3.0) * Pe; //! rotational Peclet number -> Arquivo de configuracao
extern double phi;
extern bool periodicity;
extern double *diam;
extern double *betaVec;



#endif /* SRC_HEADERS_GLOBALS_HPP_ */