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
#include <iomanip>
#include <randomic.hpp>
#include <random>
#include <math.h>

void randomic(double start, double end, int total, double *rndNArray){

// Seed with a real random value, if available
    std::random_device r;
// Choose a random mean between start and end
    std::default_random_engine e1(r());
    std::uniform_real_distribution<double> uniform_dist(start, end);
   
    // Generate a normal distribution around that mean
    std::seed_seq seed2{r()};
    std::mt19937 e2(seed2);   
  
    for (int n = 0; n < total; ++n)
        rndNArray[n] = std::round(uniform_dist(e2));
}