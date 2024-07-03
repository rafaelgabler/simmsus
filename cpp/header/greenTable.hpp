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

/**Header File
* Subroutine resposible for pre-calculating the Green functions used in the 
* calculation of periodic interactions, due to hydrodynamic forces and 
* magnetic forces and torques.           
* This procedure is done to decrease the computational cost, 
* since theses functions are pre-calculated only once and are simply
* called any time a force or torque on a particle is calculated 
* through the code in future loops.   
*/
#ifndef SRC_HEADERS_GREENTABLE_HPP_
#define SRC_HEADERS_GREENTABLE_HPP_
#include <math.h>
#include <header/config.hpp>

void greenTable(Configuration config, double sig, double e, int f, int g, double y, double **cof1,
            double **cof2, double **cof3, double **cof4, double **cof5, double **cof6, double **cof7, double **cof8);

#endif /* SRC_HEADERS_CONFIG_HPP_ */