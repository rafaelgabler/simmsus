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


/***********************
* 		 SIMMSUS
*SUBROUTINE: boxSize
*
************************/
#include <math.h>
#include <globals.hpp>

void boxSize(int N, double phi, int razao, bool initialAgregate){

    if(initialAgregate){

        double ragreg = pow((N / phi),(1.0 / 3.0));
        l = 100.0 * ragreg;
        h = l;

    } else{
        l = pow(( N / (razao * phi)) * ((4.0 * M_PI) / (3.0)),(1.0 / 3.0));       
        h = razao * l;
    }

}

