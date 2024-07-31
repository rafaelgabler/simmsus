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
#include <globals.hpp>
#include <math.h>

void repulsion(){

double r;

int i,j,q;
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            for(q = 0; q < numParticles; q++){

                if(i < q){
                    // Calculating the distance between all the particles in the suspension (pair by pair)
                        r = pow(pow((X(j,i,1)-X(j,q,1)),2.0) + pow((X(j,i,2)-X(j,q,2)),2.0) + pow((X(j,i,3)-X(j,q,3)),2.0),0.5);

                    /* If the distance is less than 2.2, i.e. less then 0.2*a, where a denotes the particle radius, 
                    then the repulsive force is turned on. It is important to notice that this force is turned off 
                    in case of particle overlap (this condition is really difficult to occur for particles 
                    wihtout inertia (mobility problem), because in this context we use a contact (Hertz) force. */

                    if(r <= 2.2 && r >= 2.0){
                        FORCAS(1,j,i,0)= 10.0 * exp(-r/0.01)*(X(j,i,1)-X(j,q,1)/r);
                        FORCAS(1,j,i,1)= 10.0 * exp(-r/0.01)*(X(j,i,2)-X(j,q,2)/r);
                        FORCAS(1,j,i,2)= 10.0 * exp(-r/0.01)*(X(j,i,3)-X(j,q,3)/r);
                    }
                    if(r >= 2.2 && r <= 2.0){
                        FORCAS(1,j,i,0) = 0.0;
                        FORCAS(1,j,i,1) = 0.0;
                        FORCAS(1,j,i,2) = 0.0;
                    }

                    // Here we calculate contact forces for overlapped particles

                    if(r <= (beta(j,i) + beta(j,q))){
                        Eij=abs(2.0-r)
                        FORCAS(2,j,i,0) = pow((100.0 * Eij),(3.0/2.0))*(X(j,i,1)-X(j,q,1))/Eij;
                        FORCAS(2,j,i,1) = pow((100.0 * Eij),(3.0/2.0))*(X(j,i,2)-X(j,q,2))/Eij;
                        FORCAS(2,j,i,2) = pow((100.0 * Eij),(3.0/2.0))*(X(j,i,3)-X(j,q,3))/Eij;
                    }else{
                        FORCAS(2,j,i,0) = 0.0;
                        FORCAS(2,j,i,1) = 0.0;
                        FORCAS(2,j,i,2) = 0.0;
                    }
                }
            }
        }
    }
}