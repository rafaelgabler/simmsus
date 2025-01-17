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

#include <header/config.hpp>

void greenTableLigaihWave3(int nGreen, double sig, double e, int f, double wave, double y, double *cof01, double *cof11,
            double *cof02, double *cof12, double *cof03, double *cof13);

void greenTableLigaihWave5(int nGreen, double sig, double e, int f, double wave, double y, double *cof01, double *cof11,
            double *cof02, double *cof12, double *cof03, double *cof13);

void greenTableTmagper3(int nGreen, double sig, double e, int f, double wave, double y, double *cof01,
            double *cof11, double *cof02, double *cof12, double *cof03, double *cof13, double *cof04, 
			double *cof14, double *cof05, double *cof15, double *cof07, double *cof17);

void greenTableTmagper5(int nGreen, double sig, double e, int f, double wave, double y, double *cof01,
            double *cof11, double *cof02, double *cof12, double *cof03, double *cof13, double *cof04, 
			double *cof14, double *cof05, double *cof15, double *cof07, double *cof17);

void greenTableFmagper3(int nGreen, double sig, double e, int f, double wave, double y, double *cof01,
            double *cof11, double *cof02, double *cof12, double *cof03, double *cof13,
			double *cof06, double *cof16, double *cof08, double *cof18);

void greenTableFmagper5(int nGreen, double sig, double e, int f, double wave, double y, double *cof01,
            double *cof11, double *cof02, double *cof12, double *cof03, double *cof13,
			double *cof04, double *cof14, double *cof06, double *cof16, double *cof08, double *cof18);


#endif /* SRC_HEADERS_GREENTABLE_HPP_ */