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

/**Header File - Create the Configuration Object With the values from config.yaml input file */
#ifndef SRC_HEADERS_CONFIG_HPP_
#define SRC_HEADERS_CONFIG_HPP_

#include <string>

class Configuration {

	private:
	    //SOLVE TORQUE EQUATION
        bool ste;
        //MAGNETIC PARTICLES
        bool mp;
        //STATIC (T) OR DYNAMIC (F) SIMULATION
        bool sds;
        //SEDIMENTATION
        bool sedimentation;
        //FLUCTUATION MODE
        bool flucmode;
        //BROWNIAN MOTION
        bool bm;
        //ACCOUNT HYDRODYNAMIC INTERACTIONS
        bool accounthi;
        //CONTINUE AN OLD SIMULATION
        bool continueosim;
        //PARTICLE INERTIA
        bool pi;

        //DIPOLAR INTERACTION CALCULATIONS 

        //PERIODIC MAGNETIC TORQUES
        bool pmt;
        //PERIODIC MAGNETIC FORCES
        bool pmf;

        //RECORD DATA OPTIONS

        //RECORD POSITION IN FILE
        bool recordposfile;
        //RECORD VELOCITY IN FILE
        bool recordvelfile;
        //RECORD DIPOLE IN FILE
        bool recorddipfile;
        //OVITO FORMAT
        bool of;
        //TECPLOT FORMAT
        bool tf;

        //SIMULATION QUANTITIES

        //NUMBER OF PARTICLES
        int numpart;
        //NUMBER OF REALIZATIONS
        int numreal;
        

        //INITIAL CONFIGURATION INFORMATION

        //ORDERED (T) OR RANDOM (F) ARRANGEMENT
        bool ordrandarr;
        //MIX MAGNETIC AND NON MAGNETIC PARTICLES
        bool mixmagnonmagpart;
        //INITIAL SPHERICAL AGGREGATE
        bool initialspheraggr;
        //INITIALLY ORDERED DIPOLES
        bool ordereddipoles;

        //PARTICLE DISTRIBUTION DATA

        //MONO (F) OR POLIDISPERSITY (T)
        bool monopolidisp;
        //VOLUME FRACTION OF PARTICLES
        double volumefracpart;
        //BOX ASPECT RATIO
        int boxaspectratio;
        //PERCENTAGE OF NON-MAGNETIC PARTICLES
        double percentnonmagpart;

        // MAGNETIC FIELD INFORMATION

        //APPLY AN EXTERNAL MAGNETIC FIELD
        bool applyextmagfield;
        //POSITION OF THE EXTERNAL FIELD
        double posextfield;
        //OSCILLATORY FIELD
        bool oscillatoryfield;
        //ROTATING FIELD
        bool rotatingfield;
        //NON-LINEAR DUFFING FIELD EXCITATION
        bool nonlineardufffieldexci;
        //DOUBLE FREQUENCY FIELD EXCITATION
        bool doublefreqfieldexci;
        //DYNAMICAL INCREASE OF FIELD FREQUENCY
        bool dynamicincrff;
        //FREQUENCY 1 OF THE MAGNETIC FIELD
        double frequency1magfield;
        //FREQUENCY 2 OF THE MAGNETIC FIELD
        double frequency2magfield;
        //C1 PARAMETER FOR DUFFING EXCITATION
        double c1parduffexc;
        //C2 PARAMETER FOR DUFFING EXCITATION
        double c2parduffexc;
        //C3 PARAMETER FOR DUFFING EXCITATION
        double c3parduffexc;
        //C4 PARAMETER FOR DUFFING EXCITATION
        double c4parduffexc;
        //MAX FREQUENCY FOR DYNAMICAL INCREASE
        double maxfreqdynincr;
        //NUMBER OF INTERVALS FOR DYN. INCREASE
        int numberintdynincr;

        //SHEAR RATE INFORMATION

        //TURN ON SHEAR RATE
        bool turnonshrate;
        //OSCILLATORY SHEAR
        bool oscillatorysh;
        //DYNAMICAL INCREASE OF SHEAR RATE
        bool dynincrshrate;
        //DIMENSIONLESS SHEAR RATE
        double dimensionlessshrate;
        //FREQUENCY FOR THE OSCILLATORY SHEAR
        double frequencyoscillsh;

        // PHYSICAL PARAMETERS

        //LAMBDA
        double lambda;
        //ALPHA
        double alpha;
        //BROWNIAN PECLET NUMBER
        double brownianpecletnum;
        //TRANSLATIONAL STOKES NUMBER
        double translationalstokesnum;

        //NUMERICAL DATA

        //SIMULATION TIME
        double simulationtime;
        //NUMERICAL TIME-STEP
        double numtimestep;
        //STEP FOR STORING THE RESULTS
        int stepstoringresults;
        //CONTINUE FROM ITERACTION NUMBER
        int continuefitnum;

        //OPTIONAL DATA TREATMENT

        //STATISTICAL ANALYSIS
        bool statanalysis;
        //CALCULATE THE STRUCTURE FACTOR
        bool calcstructfactor;
        //PRINT LOCAL MAPS OF PHI
        bool printlocalmapsphi;

	public:
		//Constructor
		Configuration();

		~Configuration();
        
        Configuration(char *configFile);

		//Getters and Setters
		void setSte(bool s_ste);
        bool getSte(){ return ste;}
        
        void setMp(bool s_mp);
        bool getMp(){ return mp;}
        
        void setSds(bool s_sds);
        bool getSds(){ return sds;}
        
        void setSedimentation(bool s_sedimentation);
        bool getSedimentation(){ return sedimentation;}
        
        void setFlucmode(bool s_flucmode);
        bool getFlucmode(){ return flucmode;}

        void setBm(bool s_bm);
        bool getBm(){ return bm;}

        void setAccounthi(bool s_accounthi);
        bool getAccounthi(){ return accounthi;}

        void setContinueosim(bool s_continueosim);
        bool getContinueosim(){ return continueosim;}

        void setPi(bool s_pi);
        bool getPi(){ return pi;}
        
        void setPmt(bool s_pmt);
        bool getPmt(){ return pmt;}

        void setPmf(bool s_pmf);
        bool getPmf(){ return pmf;}
        
        void setRecordposfile(bool s_recordposfile);
        bool getRecordposfile(){ return recordposfile;}
        
        void setRecordvelfile(bool s_recordvelfile);
        bool getRecordvelfile(){ return recordvelfile;}
        
        void setRecorddipfile(bool s_recorddipfile);
        bool getRecorddipfile(){ return recorddipfile;}
        
        void setOf(bool s_of);
        bool getOf(){ return of;}

        void setTf(bool s_tf);
        bool getTf(){ return tf;}

        void setNumpart(int s_numpart);
        double getNumpart(){ return numpart;}
        
        void setNumreal(int s_numreal);
        double getNumreal(){ return numreal;}
        
        void setOrdrandarr(bool s_ordrandarr);
        bool getOrdrandarr(){ return ordrandarr;}

        void setMixmagnonmagpart(bool s_mixmagnonmagpart);
        bool getMixmagnonmagpart(){ return mixmagnonmagpart;}       
        
        void setInitialspheraggr(bool s_initialspheraggr);
        bool getInitialspheraggr(){ return initialspheraggr;}       

        void setOrdereddipoles(bool s_ordereddipoles);
        bool getOrdereddipoles(){ return ordereddipoles;}       

        void setMonopolidisp(bool s_monopolidisp);
        bool getMonopolidisp(){ return monopolidisp;}       

        void setVolumefracpart(double s_volumefracpart);
        double getVolumefracpart(){ return volumefracpart;}
        
        void setBoxaspectratio(int s_boxaspectratio);
        int getBoxaspectratio(){ return boxaspectratio;}
        
        void setPercentnonmagpart(double s_percentnonmagpart);
        double getPercentnonmagpart(){ return percentnonmagpart;}        
        
        void setApplyextmagfield(bool s_applyextmagfield);
        bool getApplyextmagfield(){ return applyextmagfield;}

        void setPosextfield(double s_posextfield);
        double getPosextfield(){ return posextfield;}        
        
        void setOscillatoryfield(double s_oscillatoryfield);
        double getOscillatoryfield(){ return oscillatoryfield;}   

        void setRotatingfield(bool s_rotatingfield);
        bool getRotatingfield(){ return rotatingfield;}   
        
        void setNonlineardufffieldexci(bool s_nonlineardufffieldexci);
        bool getNonlineardufffieldexci(){ return nonlineardufffieldexci;}   
        
        void setDoublefreqfieldexci(bool s_doublefreqfieldexci);
        bool getDoublefreqfieldexci(){ return doublefreqfieldexci;}   

        void setDynamicincrff(bool s_dynamicincrff);
        bool getDynamicincrff(){ return dynamicincrff;}   

        void setFrequency1magfield(double s_frequency1magfield);
        double getFrequency1magfield(){ return frequency1magfield;}   

        void setFrequency2magfield(double s_frequency2magfield);
        double getFrequency2magfield(){ return frequency2magfield;}   
        
        void setC1parduffexc(double s_c1parduffexc);
        double getC1parduffexc(){ return c1parduffexc;}   

        void setC2parduffexc(double s_c2parduffexc);
        double getC2parduffexc(){ return c2parduffexc;}   

        void setC3parduffexc(double s_c3parduffexc);
        double getC3parduffexc(){ return c3parduffexc;}   

        void setC4parduffexc(double s_c4parduffexc);
        double getC4parduffexc(){ return c4parduffexc;}   

        void setMaxfreqdynincr(double s_maxfreqdynincr);
        double getMaxfreqdynincr(){ return maxfreqdynincr;}   

        void setNumberintdynincr(int s_numberintdynincr);
        int getNumberintdynincr(){ return numberintdynincr;}   

        void setTurnonshrate(bool s_turnonshrate);
        bool getTurnonshrate(){ return turnonshrate;}          

        void setOscillatorysh(bool s_oscillatorysh);
        bool getOscillatorysh(){ return oscillatorysh;}   
        
        void setDynincrshrate(bool s_dynincrshrate);
        bool getDynincrshrate(){ return dynincrshrate;}   
        
        void setDimensionlessshrate(double s_dimensionlessshrate);
        double getDimensionlessshrate(){ return dimensionlessshrate;}

        void setFrequencyoscillsh(double s_frequencyoscillsh);
        double getFrequencyoscillsh(){ return frequencyoscillsh;}

        void setLambda(double lambda);
        double getLambda(){ return lambda;}

        void setAlpha(double s_alpha);
        double getAlpha(){ return alpha;}

        void setBrownianpecletnum(double s_brownianpecletnum);
        double getBrownianpecletnum(){ return brownianpecletnum;}

        void setTranslationalstokesnum(double s_translationalstokesnum);
        double getTranslationalstokesnum(){ return translationalstokesnum;}
        
        void setSimulationtime(double s_simulationtime);
        double getSimulationtime(){ return simulationtime;}

        void setNumtimestep(double s_numtimestep);
        double getNumtimestep(){ return numtimestep;}

        void setStepstoringresults(int s_stepstoringresults);
        int getStepstoringresults(){ return stepstoringresults;}

        void setContinuefitnum(int s_continuefitnum);
        int getContinuefitnum(){ return continuefitnum;}

        void setStatanalysis(bool s_statanalysis);
        bool getStatanalysis(){ return statanalysis;}

        void setCalcstructfactor(bool s_calcstructfactor);
        bool getCalcstructfactor(){ return calcstructfactor;}

        void setPrintlocalmapsphi(bool s_printlocalmapsphi);
        bool getPrintlocalmapsphi(){ return printlocalmapsphi;}
               
};
#endif /* SRC_HEADERS_CONFIG_HPP_ */