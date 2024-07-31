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
#include <randomic.hpp>
#include <math.h>
#include <globals.hpp>

void brownian(bool gravidade, bool shear, bool torque){
int total = 3 * numParticles * numRealizations;

double *nr = new double[total];
randomic(-1.0,1.0,total,nr);
double nr1,nr2,nr3;
int i,j;

if(gravidade || shear){
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            nr1 = nr[(i * 2 + (i-2) + (numParticles * 3 * (j-1)))];
            nr2 = nr[(i * 2 + (i-1) + (numParticles * 3 * (j-1)))];
            nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j-1)))];

            double modrand = pow(pow(nr1,2.0) + pow(nr2,2.0) + pow(nr3,2.0),0.5);
            nr1 /= modrand;
            nr2 /= modrand;
            nr3 /= modrand;        
            double resPow = (pow(beta(j,i),(-2.0))) * (pow((6.0 / (Pe * dt)),0.5));
            FORCAS(6,j,i,1)=  resPow * nr1;
            FORCAS(6,j,i,2)=  resPow * nr2;
            FORCAS(6,j,i,3)=  resPow  * nr3;
        }
    }
}else{ 
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            nr1 = nr[(i * 2 + (i-2) + (numParticles * 3 * (j-1)))];
            nr2 = nr[(i * 2 + (i-1) + (numParticles * 3 * (j-1)))];
            nr3 = nr[(i * 2 + (i) + (numParticles * 3 * (j-1)))];

            double modrand = pow(pow(nr1,2.0) + pow(nr2,2.0) + pow(nr3,2.0),0.5);
            nr1 /= modrand;
            nr2 /= modrand;
            nr3 /= modrand;        
            double resPow = (pow(beta(j,i),(-2.0))) * (pow((6.0 / dt),0.5));
            FORCAS(6,j,i,1) = resPow * nr1;
            FORCAS(6,j,i,2) = resPow * nr2;
            FORCAS(6,j,i,3) = resPow * nr3;
            }
        }
    }   

if(torque){
    randomic(-1.0,1.0,(3 * numParticles * numRealizations),diarand);
    
    if(gravidade || shear){
        for(j = 0; j < numRealizations; j++){
            for(i = 0; i < numParticles; i++){
                nr1 = nr((i*2+(i-2)+(N*3*(j-1))));
                nr2 = nr((i*2+(i-1)+(N*3*(j-1))));
                nr3 = nr((i*2+(i)+(N*3*(j-1))));
                double modrand = pow(pow(nr1,2.0) + pow(nr2,2.0) + pow(nr3,2.0),0.5);            
                nr1 /= modrand;
                nr2 /= modrand;
                nr3 /= modrand;
                double resPow = pow(beta(j,i),(-2.0)) * pow(9.0 / (2.0 * Pe * dt),0.5);            
                TORQUES(3,j,i,1)= resPow * nr1;
                TORQUES(3,j,i,2)= resPow * nr2;
                TORQUES(3,j,i,3)= resPow * nr3;
            }
        }
    }else{
        for(j = 0; j < numRealizations; j++){
            for(i = 0; i < numParticles; i++){
                nr1 = nr((i*2+(i-2)+(N*3*(j-1))));
                nr2 = nr((i*2+(i-1)+(N*3*(j-1))));
                nr3 = nr((i*2+(i)+(N*3*(j-1))));
                double modrand = pow(pow(nr1,2.0) + pow(nr2,2.0) + pow(nr3,2.0),0.5);            
                nr1 /= modrand;
                nr2 /= modrand;
                nr3 /= modrand;
                double resPow = pow(beta(j,i),(-2.0)) * pow(9.0 / (2.0 * dt),0.5);
                TORQUES(3,j,i,1)= resPow * nr1;
                TORQUES(3,j,i,2)= resPow * nr2;
                TORQUES(3,j,i,3)= resPow * nr3;
            }
        }
    }
}