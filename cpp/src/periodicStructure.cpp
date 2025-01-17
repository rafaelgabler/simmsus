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

/**************************************************
* 		     SIMMSUS			   
* SUBROUTINE: estrutura_periodica		  
* Last update: 16/07/2023			  
*************************************************
* Suroutine responsible for creating a periodic   
* strucutre (lattice) containing copies of a phy- 
* sical cell with the interacting particles repli-
* cated through the real and reciprocal spaces    
* using the Ewald summation technique (1921) des- 
* cribed in details in the works of Beenakker     
* (1986), Abade (2002), Gontijo (2013)            
*************************************************/
#include <math.h>
#include <header/globals.hpp>
#include <header/constants.hpp>
#include <header/periodicStructure.hpp>
#include <iostream>
using namespace std;
void periodicStructure(){
    
// Creating the lattices indeces in the real space

int auxper[5]{};
int auxper2[5]{};
int auxper3[5]{};
int auxperi[3]{};
int auxperi2[3]{};
int auxperi3[3]{};

int a,b,c,s;

auxper[0] = 1;
auxper[1] = 2;
auxper[2] = 3;
auxper[3] = 4;
auxper[4] = 5;

auxperi[0] = 1;
auxperi[1] = 2;
auxperi[2] = 3;

auxper2[0] = 0;
auxper2[1] = 5;
auxper2[2] = 10;
auxper2[3] = 15;
auxper2[4] = 20;


auxperi2[0] = 0;
auxperi2[1] = 3;
auxperi2[2] = 6;

auxper3[0] = 0;
auxper3[1] = 25;
auxper3[2] = 50;
auxper3[3] = 75;
auxper3[4] = 100;

auxperi3[0] = 0;
auxperi3[1] = 9;
auxperi3[2] = 18;    

for(a = 0; a < 5; a++){
    for(b = 0; b < 5; b++){
        for(c = 0; c < 5; c++){

            // Number of physical boxes
            s = auxper3[a] +  auxper2[b] + auxper[c];
            ILF0[s - 1] = a;
            ILF1[s - 1] = b;
            ILF2[s - 1] = c;
        }
    }
}

// // Creating the initial configuration of all the physical lattices

// for(a = 0; a < nb; a++){
//     for(b = 0; b < numRealizations; b++){
//         for(c = 0; c < numParticles; c++){
//             XI0[a + numRealizations * (b + numParticles * c)] = X0[b * numParticles + c] + ILF0[a] * l;
//             XI1[a + numRealizations * (b + numParticles * c)] = X1[b * numParticles + c] + ILF1[a] * l;
//             XI2[a + numRealizations * (b + numParticles * c)] = X2[b * numParticles + c] + ILF2[a] * h;            
//         }
//     }
// }

// Creating the lattice's indeces in the reciprocal space
for(a =  0; a < 3; a++){
    for(b = 0; b < 3; b++){
        for(c = 0; c < 3; c++){
            // Number of reciprocal lattices
            s = auxperi3[a] + auxperi2[b] + auxperi[c];  
            cout << s << "\t";
            ILR0[s - 1] = a;
            ILR1[s - 1] = b;
            ILR2[s - 1] = c;
        }
    }
}
}