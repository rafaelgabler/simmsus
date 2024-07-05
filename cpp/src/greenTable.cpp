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

/**CPP File
* Subroutine resposible for pre-calculating the Green functions used in the 
* calculation of periodic interactions, due to hydrodynamic forces and 
* magnetic forces and torques.           
* This procedure is done to decrease the computational cost, 
* since theses functions are pre-calculated only once and are simply
* called any time a force or torque on a particle is calculated 
* through the code in future loops.   
*/
#include <math.h>
#include <header/config.hpp>


void greenTableLigaihWave5(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;

double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);
//ligaih ACCOUNT HYDRODYNAMIC INTERACTIONS......: FALSE
//fmagper PERIODIC MAGNETIC FORCES
//tmagper PERIODIC MAGNETIC TORQUES


for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate hydrodynamic interactions (in real space)

	cof1[1][i] = (0.75  / RR1 + 0.5 / RR3) * erfc(sig * RR1) + (4.0 * SIG7 * RR4 + 3.0 * SIG3 * RR2 - 
		20.0 * SIG5 * RR2 - 4.50 * sig + 14.0 * SIG3 + sig / RR2) * exp(-SIG2 * RR2) / sqrtpi; 
	cof2[1][i] = (0.75 / RR1 - 1.5 / RR3) * erfc(sig * RR1) + (-4.0 * SIG7 * RR4 - 3.0 * SIG3 * RR2 +
		16.0 * SIG5 * RR2 + 1.50 * sig - 2.0 * SIG3 - 3.0 * sig / RR2) * exp(-SIG2 * RR2) / sqrtpi;

// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers

	cof3[0][i] = twpie + ((i) * 4.0 * cof3pow);	

// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate hydrodynamic interactions (in reciprocal space)	

	cof3[1][i] = (1.0 / (e * e * y)) * (1.0 - 1.0 * RR2 / 3.0) * 
	(1.0 + 0.25  * RR2 / SIG2 + 0.125 * RR4 / SIG4) * 6.0 * M_PI * exp(-0.25 * RR2 / SIG2) / RR2;

	}
}

void greenTableLigaihWave3(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;

double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);
//ligaih ACCOUNT HYDRODYNAMIC INTERACTIONS......: FALSE
//fmagper PERIODIC MAGNETIC FORCES
//tmagper PERIODIC MAGNETIC TORQUES


for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate hydrodynamic interactions (in real space)

	cof1[1][i] = (0.75  / RR1 + 0.5 / RR3) * erfc(sig * RR1) + (4.0 * SIG7 * RR4 + 3.0 * SIG3 * RR2 - 
		20.0 * SIG5 * RR2 - 4.50 * sig + 14.0 * SIG3 + sig / RR2) * exp(-SIG2 * RR2) / sqrtpi; 
	cof2[1][i] = (0.75 / RR1 - 1.5 / RR3) * erfc(sig * RR1) + (-4.0 * SIG7 * RR4 - 3.0 * SIG3 * RR2 +
		16.0 * SIG5 * RR2 + 1.50 * sig - 2.0 * SIG3 - 3.0 * sig / RR2) * exp(-SIG2 * RR2) / sqrtpi;

// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers	
	cof3[0][i] = twpie + ((i) * 2.0 * cof3pow);


// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate hydrodynamic interactions (in reciprocal space)	

	cof3[1][i] = (1.0 / (e * e * y)) * (1.0 - 1.0 * RR2 / 3.0) * 
	(1.0 + 0.25  * RR2 / SIG2 + 0.125 * RR4 / SIG4) * 6.0 * M_PI * exp(-0.25 * RR2 / SIG2) / RR2;

	}
}

void greenTableTmagper5(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3, double **cof4, double **cof5, double **cof7){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;
double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);

for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate periodic magnetic torques (in real space)	

	cof4[0][i] = cof1[0][i]; 	    
	cof4[1][i] = (erfc(sig * RR1)+ ((2.0 * sig * RR1)/(sqrtpi)) * exp(-SIG2 * RR2)) / RR3; 
	cof5[0][i] =  cof1[0][i];
	cof5[1][i] = (3.0 * erfc(sig * RR1) + (((2.0 * sig * RR1) / (sqrtpi)) * (3.0 + (2.0 * SIG2 * RR2)) * exp(-SIG2 * RR2))) / RR5;		
    
// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers

	cof3[0][i] = twpie + ((i) * 4.0 * cof3pow);

// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate periodic magnetic torques (in reciprocal space)
	cof7[0][i] = cof3[0][i];
	cof7[1][i] = -(1.0/(e * e * y)) * (4.0 * M_PI * (exp(-(pow((M_PI/e),2.0)) * RR2 / SIG2)))/1.0;
	}
}

void greenTableTmagper3(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3, double **cof4, double **cof5, double **cof7){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;
double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);

for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate periodic magnetic torques (in real space)	

	cof4[0][i] = cof1[0][i]; 	    
	cof4[1][i] = (erfc(sig * RR1)+ ((2.0 * sig * RR1)/(sqrtpi)) * exp(-SIG2 * RR2)) / RR3; 
	cof5[0][i] =  cof1[0][i];
	cof5[1][i] = (3.0 * erfc(sig * RR1) + (((2.0 * sig * RR1) / (sqrtpi)) * (3.0 + (2.0 * SIG2 * RR2)) * exp(-SIG2 * RR2))) / RR5;		
    
// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers
	cof3[0][i] = twpie + ((i) * 2.0 * cof3pow);

// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate periodic magnetic torques (in reciprocal space)	

	cof7[0][i] = cof3[0][i];
	cof7[1][i] = -(1.0/(e * e * y)) * (4.0 * M_PI * (exp(-(pow((M_PI/e),2.0)) * RR2 / SIG2)))/1.0;
	}
}

void greenTableFmagper5(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3, double **cof4, double **cof6, double **cof8){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;
double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);

for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate periodic magnetic forces (in real space)	

	cof6[0][i] = cof1[0][i];
	cof6[1][i] = (15.0 * erfc(sig * RR1) + ((2.0 * sig * RR1) / (sqrtpi)) * (15.0 + (10.0 * SIG2 * RR2 + 4.0 * SIG4 * RR4)) * exp(-SIG2 * RR2)) / RR7;		


// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers

	if(wave == 5.0){
	    cof3[0][i] = twpie + ((i) * 4.0 * cof3pow) ;
	}else if(wave == 3.0){
	    cof3[0][i] = twpie + ((i) * 2.0 * cof3pow);
	}

// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate periodic magnetic forces (in reciprocal space)

	cof8[0][i] = cof3[0][i];	
	cof8[1][i]=((8.0 * pow(M_PI,2.0)) / (pow(e,4.0) * (RR2))) * exp(-(pow(M_PI,2.0) * RR2) / 
	(SIG2 * pow(e,2.0))) / RR2; 
	}
}

void greenTableFmagper3(double sig, double e, int f, double wave, double y, double **cof1,
            double **cof2, double **cof3, double **cof6, double **cof8){

// auxiliary parameters used to compute these functions
double RR1,RR2,RR3,RR4,RR7,RR5, SIG2,SIG3,SIG4,SIG5,SIG7;
double twpie =  (2.0 * M_PI / e);
double cof3pow = ((((pow(3.0,0.5)) * M_PI / e) - (2.0 * M_PI / e)) / 9999.0);
double cof1pow = ((pow(3.0,0.5) * e * (pow(f,(1.0 / 3.0)))) / 9999.0);
double sqrtpi = sqrt(M_PI);

for (int i = 0; i < 10000; i++){

    cof1[0][i] = (2.0 + i) * cof1pow;
    cof2[0][i] = cof1[0][i];

// Calculating the powers of all the possible distances between the particles

	RR1 = cof1[0][i];
	RR2 = pow(cof1[0][i],2.0);
	RR3 = pow(cof1[0][i],3.0);
	RR4 = pow(cof1[0][i],4.0);
	RR5 = pow(cof1[0][i],5.0);
	RR7 = pow(cof1[0][i],7.0);

	SIG2 = pow(sig,2.0);
	SIG3 = pow(sig,3.0);
	SIG4 = pow(sig,4.0);
	SIG5 = pow(sig,5.0);
	SIG7 = pow(sig,7.0);

// Computing the functions used to calculate periodic magnetic forces (in real space)	

	cof6[0][i] = cof1[0][i];
	cof6[1][i] = (15.0 * erfc(sig * RR1) + ((2.0 * sig * RR1) / (sqrtpi)) * (15.0 + (10.0 * SIG2 * RR2 + 4.0 * SIG4 * RR4)) * exp(-SIG2 * RR2)) / RR7;		


// Same procedure done now for the reciprocal space
// Calculating all possible wave numbers

	cof3[0][i] = twpie + ((i) * 2.0 * cof3pow);

// Determining the 2th and 4th powers of these possible wave numbers
	RR2 = pow(cof3[0][i],2.0);
	RR4 = pow(cof3[0][i],4.0);

// Computing the functions used to calculate periodic magnetic forces (in reciprocal space)

	cof8[0][i] = cof3[0][i];	
	cof8[1][i]=((8.0 * pow(M_PI,2.0)) / (pow(e,4.0) * (RR2))) * exp(-(pow(M_PI,2.0) * RR2) / 
	(SIG2 * pow(e,2.0))) / RR2; 
	}
}