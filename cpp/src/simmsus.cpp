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
#include <math.h>
#include <config.hpp>
#include <constants.hpp>
#include <globals.hpp>
#include <particleDistribution.hpp>
#include <boxSize.hpp>
#include <greenTable.hpp>
#include <brownian.hpp>
#include <repulsion.hpp>

int main(int argc, char **argv){

//Read the Input File and Create the Object
Configuration configuration(argv[1]);

// Initialization some important variables with 0.0
int i,j;
double X = 0.0;
double U = 0.0;
double aux1 = 0.0;
double aux2 = 0.0;
double aux3 = 0.0;
double Di = 0.0;
double nr = 0.0;
double FORCAS = 0.0;
double FT = 0.0;
double hidrodinamica_aux1 = 0.0;
double hidrodinamica_aux2 = 0.0;
double hidro1 = 0.0;
double hidro2 = 0.0;
double contribuicao_self = 0.0;
double contribuicao_fisico = 0.0;
double contribuicao_reciproco = 0.0;
double torquereal = 0.0;
double torquereciproco = 0.0;
double auxt = 0.0;
bool periodicity = false;
double wave;

// double phi = configuration.getPh
int numParticles = configuration.getNumpart();
numRealizations = configuration.getNumreal();
totalRealParticle = numParticles * numRealizations;
//"estatica" = TRUE: Monte carlo simulation (number of time-steps = 2)
//"estatica" = FALSE: Dynamic simulation (number of time-steps = time/dt)

if(configuration.getStatanalysis()){
    npast = 2;
}else{
    npast = configuration.getSimulationtime()  / configuration.getNumtimestep();
}

// Number of random numbers used in each time-step
nnr = 3 * totalRealParticle;

// Allocating variables in the memory
diarand = new double[totalRealParticle];
//call allocatevariables
 
// Defining number formats for the output files

//1012 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)
//2024 FORMAT(E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2,3x,E11.4e2)

if(configuration.getAccounthi() || configuration.getPmf() || configuration.getPmt()){
    periodicity = true;
}

// Defining the size of the particles

particleDistribution(poliDispersity, numRealizations, numParticles, DIAM, diarand, beta);

// Calculating the size of the simulation box based on the 
// number of particles and on the volume fraction defined
// by the user in the simconfig.dat

boxSize(numParticles, configuration.getVolumefracpart(), configuration.getBoxaspectratio(), configuration.getInitialspheraggr());

double qsi = 1.0 * (pow(M_PI,0.5)  / (pow(l * l * h,(1.0/3.0))));

// Creating the files to write the results

// if(!continua){
//  gera_arquivos(posicao,velocidade,numRealizations)
// }

// Defining local simulation parameters
// convergence parameter for the periodic sums
 

// Building a table with all the Green-functions
// And building the periodic structure to compute
// the Ewald summations

if(periodicity){
wave = pow(nbr,(1.0/3.0));

if(wave == 5.0){
    if(configuration.getAccounthi()){
        greenTableLigaihWave5(qsi, l, nb, wave, h, **cof1, **cof2, **cof3);
    }else if(configuration.getPmt()){
        greenTableTmagper5(qsi, l, nb, wave, h, **cof1, **cof2, **cof3, **cof4, **cof5, **cof7);
    }else if(configuration.getPmf()){
        greenTableFmagper3(qsi, l, nb, wave, h, **cof1, **cof2, **cof3, **cof4, **cof6, **cof8);
    }
}else if(wave == 3.0){
    if(configuration.getAccounthi()){
        greenTableLigaihWave3(qsi, l, nb, wave, h, **cof1, **cof2, **cof3);
    }else if(configuration.getPmt()){
        greenTableTmagper3(qsi, l, nb, wave, h, **cof1, **cof2, **cof3, **cof4, **cof5, **cof7);
    }else if(configuration.getPmf()){
        greenTableFmagper3(qsi, l, nb, wave, h, **cof1, **cof2, **cof3, **cof4, **cof6, **cof8);
    }
}
 
 //call estrutura_periodica
}

// Creating the initial distribution of the particles dipole moments

// if(magpart){
// if(!continua){
//  call distribui_dipolo(Di,numRealizations,numParticles)
// }
// }

// Creating the initial particle distribution

//  call condicao_inicial

std::cout << "*,'******************************************************************************'";
std::cout << "*,'*                                                                            *'";
std::cout << "*,'*                   INITIAL CONDITIONS SUCCESSFULLY GENERATED                *'";
std::cout << "*,'*                                                                            *'";
std::cout << "*,'******************************************************************************'";


//  call field_excitations

// Checking if we are continuing an old simulation or starting a new one

// if(continua){
// iter=iter
// auxiliar_continua=npast
// else
// iter=1
// auxiliar_continua=npast-1
// }
// aux_real=auxiliar_continua

// Printing initial fields of the local volume fraction

// if(printphi){
// call campo_phi(numRealizations,k)   
// }

// Beggining of the main loop

std::cout << "'******************************************************************************'";
std::cout << "'*                                                                            *'";
std::cout << "'*                                SIMULATING                                  *'";
std::cout << "'*                                                                            *'";
std::cout << "'******************************************************************************'";


for(int k = iter; k < auxiliar_continua; k++){

        int k_real = k;

        // Updating the instant shear-rate 

        if(shear || oscillatory || bifshear){
            shearrate = gpvetor(k) * sin(freq * dt * k);
        }else{
            shearrate *= sin(freq * dt * k);
        }        
}

// Calculating Brownian forces and torques
if(configuration.getBm()){
    brownian(configuration.getSedimentation(), configuration.getTurnonshrate(), configuration.getSte());
}

// Calculating gravitational forces
for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){    
            FORCAS(2,j,i,0)=0.0;
            FORCAS(2,j,i,1)=0.0;
            FORCAS(2,j,i,2)=0.0;
        }
    }
if(configuration.getSedimentation()){
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            FORCAS(2,j,i,2)= pow(-beta(j,i),3.0);
        }
    }
}

// Calculating repulsion forces between all particles and
// contact forces between overlapped particles

 call repulsion

// Calculating non-periodic long-range dipolar forces between
// the particles

if(magpart && !fmagper){
 call forca_magnetica
}


// Calculating magnetic forces due to an external magnetic field
// At this point we have been considering that magnetic field
// gradients are null within the suspension space due to the
// scale of the simulation box, which is a continuum volume point.
// This way the magnetic force due to an external magnetic field
// should be null, but if the user wants to simulate the effects
// of a non-null magnetic field gradient within the suspension
// space,{ the commented command "call campo_externo" should
// be enabled.

//if(externo){
// call campo_externo
//else
FORCAS(5,:,:,:) = 0.0;
//}

// Calculating the total force acting on each particle 
for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        FT(j,i,0) = FORCAS(0,j,i,0) + FORCAS(1,j,i,0) + FORCAS(2,j,i,0) + FORCAS(3,j,i,0) + FORCAS(4,j,i,0) + FORCAS(5,j,i,0);
        FT(j,i,1) = FORCAS(0,j,i,1) + FORCAS(1,j,i,1) + FORCAS(2,j,i,1) + FORCAS(3,j,i,1) + FORCAS(4,j,i,1) + FORCAS(5,j,i,1);
        FT(j,i,2) = FORCAS(0,j,i,2) + FORCAS(1,j,i,2) + FORCAS(2,j,i,2) + FORCAS(3,j,i,2) + FORCAS(4,j,i,2) + FORCAS(5,j,i,2);
    }
} 

// Calculating periodic interactions

if(periodicidade){
 periodic_interactions();
}

// Calculating the particles velocities
if(!ligaih){
    if(inertia){ 
        for(j = 0; j < numRealizations; j++){
            for(i = 0; i < numParticles; i++){
                resvel(U(j,i,0),dt,St,FT(j,i,0));
                resvel(U(j,i,1),dt,St,FT(j,i,1));
                resvel(U(j,i,2),dt,St,FT(j,i,2));
            }
        } 
}
else {
    U=FT;    
}

// Adding the shear contribution to the velocity of each particle

if(shear){ 
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            U(j,i,1) += shearrate * X(j,i,2);
        }
    }  
}

//***************************** FLUCTUATION MODE OPTION ***********************************//

// Here the velocity of the particles is subtracted from the average velocity of the system
// and the user will observe the behavior of the particles moving with a coordinate system
// at the average speed of the suspension and will only see flucuations induced by long-range
// interactions.

if(leito){
    for(i = 0; i < numParticles; i++){
        usistema(i,0)=sum(U(:,i,0)) / numRealizations;
        usistema(i,1)=sum(U(:,i,1)) / numRealizations;
        usistema(i,2)=sum(U(:,i,2)) / numRealizations;
    }
}
 
for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        U(j,i,0) -= usistema(i,0);
        U(j,i,1) -= usistema(i,1);
        U(j,i,2) -= usistema(i,2);
    }
} 

//*****************************************************************************************//

// Calculating the current position of the particles using a straighfoward Euler integration
 
for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        respos(X(j,i,0),dt,U(j,i,0));
        respos(X(j,i,1),dt,U(j,i,1));
        respos(X(j,i,2),dt,U(j,i,2));
    }
}

//*****************************************************************************************//
// Imposing periodic boundary conditions to avoid any dependence of the system with respect
// to physical walls

 
for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        if(X(j,i,1) >  l){
            X(j,i,1) = X(j,i,1)-l;
        }
        if(X(j,i,1) < 0.0){
            X(j,i,1) = l - X(j,i,1);
        }
        if(X(j,i,2) > l){
            X(j,i,2) = X(j,i,2) - l; 
        }
        if(X(j,i,2) < 0.0){
            X(j,i,2)=l-X(j,i,2);
        }
        if(X(j,i,3) > h){
            X(j,i,3) = X(j,i,3) - h;
        }
        if(X(j,i,3) < 0.0){
            X(j,i,3) = h - X(j,i,3);
        }
    }
}
 


//********************************************************************************************************//

// Writting in a data file the positions, velocities and orientations of the particles in the current time
// step

//  call writting_files(k,k_real)

// SOLUTION OF THE ROTATIONAL MOTION OF THE PARTICLES

if(torque){

// Calculating long-range magnetic torques on each particle due to particle interaction

if(magpart && !tmagper){
    torque_magnetico();
}
 
// Calculating the magnetic torques induced by an external field 

 if(externo){
    if(rotating){
        rotating_field(alpha, freqcampo*k*dt);
    }else{
        if(oscilacampo){
            torque_externo(alpha*campo(k));
        } else {
            torque_externo(alpha);
        }
        }
    }
 }
 }

// Brownian torques have already been computed in subroutine Browian

// Calculating the total torques acting on the particles

for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        if(browniano){
            Tt(j,i,0) = TORQUES(0,j,i,0) + TORQUES(1,j,i,0) + TORQUES(2,j,i,0);
            Tt(j,i,1) = TORQUES(0,j,i,1) + TORQUES(1,j,i,1) + TORQUES(2,j,i,1);
            Tt(j,i,2) = TORQUES(0,j,i,2) + TORQUES(1,j,i,2) + TORQUES(2,j,i,2);
        }
        else{
            Tt(j,i,0) = TORQUES(0,j,i,0) + TORQUES(1,j,i,0); 
            Tt(j,i,1) = TORQUES(0,j,i,1) + TORQUES(1,j,i,1);
            Tt(j,i,2) = TORQUES(0,j,i,2) + TORQUES(1,j,i,2); 
        }
    }
}


// Solving the angular velocity of the particles

if(mistura){
    for(j = 0; j < numRealizations; j++){
        for(i = (numParticles * percentual) + 1; i < numParticles; i++){
            resomega(W(j,i,0),dt,Str,Tt(j,i,0));
            resomega(W(j,i,1),dt,Str,Tt(j,i,1));
            resomega(W(j,i,2),dt,Str,Tt(j,i,2));
        }
    }
else{
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            resomega(W(j,i,0),dt,Str,Tt(j,i,0));
            resomega(W(j,i,1),dt,Str,Tt(j,i,1));
            resomega(W(j,i,2),dt,Str,Tt(j,i,2));
        }
    }
}

if(shear){
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            W(j,i,0) -= shearrate * 0.5;
        }
    }
}

// Evolving the dipole moment of the particles with their angular velocities

if(mistura){
    for(j = 0; j < numRealizations; j++){
        for(i = (numParticles * percentual) + 1; i < numParticles; i++){
            evoldip(Di(j,i,0),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt);
            evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,0),W(j,i,2),W(j,i,0),dt);
            evoldip(Di(j,i,2),Di(j,i,0),Di(j,i,1),W(j,i,0),W(j,i,1),dt);
        }
    }
}
else{    
    for(j = 0; j < numRealizations; j++){
        for(i = 0; i < numParticles; i++){
            evoldip(Di(j,i,0),Di(j,i,1),Di(j,i,2),W(j,i,1),W(j,i,2),dt);
            evoldip(Di(j,i,1),Di(j,i,2),Di(j,i,0),W(j,i,2),W(j,i,0),dt);
            evoldip(Di(j,i,2),Di(j,i,0),Di(j,i,1),W(j,i,0),W(j,i,1),dt);
        }
    }
}

// Normalizing the dipole moments

for(j = 0; j < numRealizations; j++){
    if(mistura){
        for(i = 0; i < (percentual * numParticles); i++){
            Di(j,i,0) = 0.0;
            Di(j,i,1) = 0.0;
            Di(j,i,2) = 0.0;
        }    
        for(i = (numParticles * percentual) + 1; i < numParticles; i++){
            modulodipolo = pow(pow(Di(j,i,0),2.0) + pow(Di(j,i,1),2.0) + pow(Di(j,i,2),2.0),0.5);
            Di(j,i,0) /= modulodipolo;
            Di(j,i,1) /= modulodipolo;
            Di(j,i,2) /= modulodipolo;
        }
    }
else{
    for(i = 0; i < numParticles; i++){
        modulodipolo = pow(pow(Di(j,i,0),2.0) + pow(Di(j,i,1),2.0) + pow(Di(j,i,2),2.0),0.5);
        Di(j,i,0) /= modulodipolo;
        Di(j,i,1) /= modulodipolo;
        Di(j,i,2) /= modulodipolo;
    }
}
}

// Calculating the magnetization of the suspension

if(grafmag){
    media_dipolo(Di,numParticles ,numRealizations,magtempo(1,k),1);
    media_dipolo(Di,numParticles,numRealizations,magtempo(2,k),2);
    media_dipolo(Di,numParticles,numRealizations,magtempo(3,k),3);

// Calculating the magnetization errorbar

for(j = 0; j < numRealizations; j++){
    for(i = 0; i < numParticles; i++){
        flutmag(i,j) = pow((Di(j,i,2)-magtempo(2,k)),2.0);
    }
}
erromag = ((1.0 / (numParticles * numRealizations)) * pow(sum(flutmag)),0.5);

// Calculating the magnetization derivatives in each direction

derivada1 = (magtempo(0,k)-magtempo(0,k-1))/dt;
derivada2 = (magtempo(1,k)-magtempo(1,k-1))/dt;
derivada3 = (magtempo(2,k)-magtempo(2,k-1))/dt;

// Writting the current magnetization components and their derivatives in a data file

// write(5*numRealizations,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt

// Writting aditional magnetization datafiles related to dynamical increase of the field's frequency

contfreqinteiro1 = ((k-1) * dt) / intervalo;
contfreqinteiro2 = (k * dt) / intervalo;

if(contfreqinteiro1.ne.contfreqinteiro2){
    multiplofreq++;
}

// if(bifurcation){ 
// write(400*numRealizations+multiplofreq+1,2024) campo(k),y(k), magtempo(1,k),magtempo(2,k),magtempo(3,k), derivada1, derivada2, derivada3, k*dt
// }
// }

tempototal(k)=k*dt
}
} //checar o for

// Computing the integral of M.dH for MHT information (to be fixed)
//if(grafmag){
//do k=2,npast
//trap(k-1)= ((magtempo(3,k-1)*y(k-1)) + (magtempo(3,k)*y(k)))*dt*0.5
//end do
//trapezio=sum(trap)
//print *,'******************************************************************************'
//print *,'*                                                                            *'
//print *,'*                                 MHT REPORT	                              *'
//print *,'*                                                                            *'
//print *,'******************************************************************************'
//write(*,*)''
//write(*,*)'INTEGRAL OF M.dH:',abs(trapezio)
//}

// Printing local volume fraction maps inside the simulation box
// if(printphi){
// call campo_phi(numRealizations,k)
// }

// Calculating the final structure factor of the suspension
if(fator){
  fator_estrutura(X,numParticles,l,h,dt,numRealizations);
}

// Closing files
// do j=1,2*numRealizations
//  close(j)
// end do

//  close(100*numRealizations)
//  close(300*numRealizations)

// Deallocating all matrices and vectors

// deallocate(X, STAT = DeAllocateStatus)
// deallocate(U, STAT = DeAllocateStatus)
// deallocate(FORCAS, STAT = DeAllocateStatus)
// deallocate(FT, STAT = DeAllocateStatus)
// deallocate(nr, STAT = DeAllocateStatus)
// deallocate(hidrodinamica_aux1, STAT = DeAllocateStatus)
// deallocate(hidrodinamica_aux2, STAT = DeAllocateStatus)
// deallocate(hidro1, STAT = DeAllocateStatus)
// deallocate(hidro2, STAT = DeAllocateStatus)
// deallocate(ILF, STAT = DeAllocateStatus)
// deallocate(ILR, STAT = DeAllocateStatus)
// deallocate(XI, STAT = DeAllocateStatus)
// deallocate(Tt, STAT = DeAllocateStatus)
// deallocate(Di, STAT = DeAllocateStatus)
// deallocate(aux1, STAT = DeAllocateStatus)
// deallocate(aux2, STAT = DeAllocateStatus)
// deallocate(aux3, STAT = DeAllocateStatus)
// deallocate(aux4, STAT = DeAllocateStatus)
// deallocate(contribuicao_self, STAT = DeAllocateStatus)
// deallocate(contribuicao_fisico, STAT = DeAllocateStatus)
// deallocate(contribuicao_reciproco, STAT = DeAllocateStatus)
// if(tmagper){
// deallocate(auxt, STAT = DeAllocateStatus)
// deallocate(torquereal, STAT = DeAllocateStatus)
// deallocate(torquereciproco, STAT = DeAllocateStatus)
// deallocate(cof4, STAT = DeAllocateStatus)
// deallocate(cof5, STAT = DeAllocateStatus)
// deallocate(cof7, STAT = DeAllocateStatus)
// }
// if(fmagper){
// deallocate(cof6, STAT = DeAllocateStatus)
// deallocate(cof8, STAT = DeAllocateStatus)
// deallocate(auxf, STAT = DeAllocateStatus)
// deallocate(forcareal, STAT = DeAllocateStatus)
// deallocate(forcareciproca, STAT = DeAllocateStatus)
// }
// deallocate(ILF, STAT = DeAllocateStatus)
// deallocate(ILR, STAT = DeAllocateStatus)
// deallocate(XI, STAT = DeAllocateStatus)
// deallocate(cof1, STAT = DeAllocateStatus)
// deallocate(cof2, STAT = DeAllocateStatus)
// deallocate(cof3, STAT = DeAllocateStatus)
// if(leito){
// deallocate(usistema, STAT = DeAllocateStatus)
// }
// if(grafmag){
// deallocate(magtempo, STAT = DeAllocateStatus)
// deallocate(flutmag, STAT = DeAllocateStatus)
// }
// deallocate(tempototal, STAT = DeAllocateStatus)
// if(agregado_inicial){
// deallocate(centro_massa, STAT = DeAllocateStatus)
// }
// deallocate(DIAM, STAT = DeAllocateStatus)
// deallocate(beta, STAT = DeAllocateStatus)
// deallocate(diarand, STAT = DeAllocateStatus)

// write(*,*) ''
// write(*,*) 'End of the processing module...'
// write(*,*) ''
// end subroutine main
return 0;

}