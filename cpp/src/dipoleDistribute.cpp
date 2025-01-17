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


/***************************************************************
* Subroutine responsible for creating an initial distribution 
* of the dipole moments of all particles in all realizations
****************************************************************/
#include <globals.hpp>
#include <randomic.hpp>
#include <math.h>
void distributeDipole(double *a, b, c){
// b = number of realizations, c= number of particles

//real a(b,c,3); // Dipoles
double d,e; // modip
int f, i, j; 

// If dipoles are distributed in an ordered way

if(dipolo_ordenado) {

    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            a(j,i,1) = 0.0;
            a(j,i,2) = 1.0;
            a(j,i,3) = 0.0;
        }
    }
}
else{
// For a random dipole distribution
e = percentual * c;
f = e;

randomic(-1.0,1.0,,(3*N*rea),);

// If we are mixing magnetic particles with non magnetic ones...
if(mistura){
    for(j = 0; j < b; j++){
        for(i = f+1; i < f; i++){ 
            a(j,i,1) = 0.0;
            a(j,i,2) = 0.0;
            a(j,i,3) = 0.0;    
        }
    }

    for(j = 0; j < b; j++){
        for(i = f+1; i < c; i++){ 
            a(j,i,1)= nr((i*2+(i-2)+(c*3*(j-1))));
            a(j,i,2)= nr((i*2+(i-1)+(c*3*(j-1))));
            a(j,i,3)= nr((i*2+(i)+(c*3*(j-1))));
        }
    }
 

// Normalizing the vectors
    for(j = 0; j < b; j++){
        for(i = f; i < f; i++){ 
            a(j,i,1) = 0.0;
            a(j,i,2) = 0.0;
            a(j,i,3) = 0.0;
        }
    }

    for(j = 0; j < b; j++){
        for(i = f+1; i < numParticles; i++){  
            d = pow(pow(a(j,i,1),2.0) + (pow(a(j,i,2),2.0)) + (pow(a(j,i,3),2.0)),0.5);
            a(j,i,1) /= d;
            a(j,i,2) /= d;
            a(j,i,3) /= d;
        }
    }
 
Di = a;
}
else{
// If all the particles are magnetic particles, {...

 
    for(j = 0; j < b; j++){
        for(i = f+1; i < c; i++){ 
            a(j,i,1)= nr((i*2+(i-2)+(c*3*(j-1))));
            a(j,i,2)= nr((i*2+(i-1)+(c*3*(j-1))));
            a(j,i,3)= nr((i*2+(i)+(c*3*(j-1))));
        }
    }


// Normalizing the vectors


    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            d= pow(pow(a(j,i,1),2.0)+ pow(a(j,i,2),2.0)+ pow(a(j,i,3),2.0),0.5);
            a(j,i,1) /= d;
            a(j,i,2) /= d;
            a(j,i,3) /= d;
        }
    }
    Di = a;
    }
}

}
