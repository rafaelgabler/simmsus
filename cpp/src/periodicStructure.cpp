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
#include <globals.hpp>
void periodicStructure(){

double powNb = pow(numBoxes,(1.0/3.0));
double powNbr = pow(numRecBoxes,(1.0/3.0));

// Creating the lattices indeces in the real space

int auxper = 0;

if(powNb == 5.0){
    auxper[0] = 1;
    auxper[1] = 2;
    auxper[2] = 3;
    auxper[3] = 4;
    auxper[4] = 5;

    auxper2[0] = 0;
    auxper2[1] = 5;
    auxper2[2] = 10;
    auxper2[3] = 15;
    auxper2[4] = 20;

    auxper3[0] = 0;
    auxper3[1] = 25;
    auxper3[2] = 50;
    auxper3[3] = 75;
    auxper3[4] = 100;

    for(int a = 0; a < powNb; a++){
        for(int b = 0; b < powNb; b++){
            for(int c = 0; c < powNb; c++){

                // Number of physical boxes
                double s = auxper3[a] +  auxper2[b] + auxper[c];                
                ILF[s][0] = a - 3;
                ILF[s][1] = b - 3;
                ILF[s][2] = c - 3;
            }
        }
    }
} else if(powNb == 3.0){
        auxper[0] = 1;
        auxper[1] = 2;
        auxper[2] = 3;

        auxper2[0] = 0;
        auxper2[1] = 3;
        auxper2[2] = 6;

        auxper3[0] = 0;
        auxper3[1] = 9;
        auxper3[2] = 18;
    

for(int a = 0; a < powNb; a++){
    for(int b = 0; b < powNb; b++){
        for(int c = 0; c < powNb; c++){

            // Number of physical boxes

            double s = auxper3[a] +  auxper2[b] + auxper[c];
            ILF[s][0] = a - 2;
            ILF[s][1] = b - 2;
            ILF[s][2] = c - 2;
            }

        }
    }
}

// Creating the initial configuration of all the physical lattices

for(int a = 0; a < numBoxes; a++){
    for(int b = 0; b < numRealizations; b++){
        for(int i = 0; i < numParticles; i++){
            XI(a,b,i,0)= X(b,i,0) + ILF[a][0] * l;
            XI(a,b,i,1)= X(b,i,1) + ILF[a][1] * l;
            XI(a,b,i,2)= X(b,i,2) + ILF[a][2] * h;
        }
    }
}

// if(k.eq.1){
// open(872,file='condicao_inicial_periodica.plt')
// write(872,*)'Variables= "X1","X2","X3"'
// do a=1,nb
// do b=1,N
// write(872,*) XI(a,1,b,1),XI(a,1,b,2),XI(a,1,b,3)
// }


// Creating the lattice's indeces in the reciprocal space

if(powNbr == 5.0){
    auxper[0] = 1;
    auxper[1] = 2;
    auxper[2] = 3;
    auxper[3] = 4;
    auxper[4] = 5;

    auxper2[0] = 0;
    auxper2[1] = 5;
    auxper2[2] = 10;
    auxper2[3] = 15;
    auxper2[4] = 20;

    auxper3[0] = 0;
    auxper3[1] = 25;
    auxper3[2] = 50;
    auxper3[3] = 75;
    auxper3[4] = 100;
    
    int powRecBox = pow(numRecBoxes,(1.0/3.0));
    for(int a = 0; a < powRecBox; a++){
        for(int b = 0; b < powRecBox; b++){
            for(int c = 0; c < powRecBox; c++){

                // Number of reciprocal lattices
                double s = auxper3[a] + auxper2[b] + auxper[c];            
                ILR[s][0] = a - 3;
                ILR[s][1] = b - 3;
                ILR[s][2] = c - 3;
            }
        }
    }

} else if(powNbr == 3.0){
    auxper[0] = 1;
    auxper[1] = 2;
    auxper[2] = 3;

    auxper2[0] = 0;
    auxper2[1] = 3;
    auxper2[2] = 6;

    auxper3[0] = 0;
    auxper3[1] = 9;
    auxper3[2] = 18;
    
    for(int a = 0; a < (pow(nbr,(1.0/3.0))); a++){
        for(int b = 0; b < (pow(nbr,(1.0/3.0))); b++){
            for(int c = 0; c < (pow(nbr,(1.0/3.0))); c++){

                // Number of reciprocal lattices
                double s = auxper3[a] + auxper2[b] + auxper[c];
                ILR[s][1] = a - 2;
                ILR[s][2] = b - 2;
                ILR[s][3] = c - 2;

            }
        }
    }   
}
}