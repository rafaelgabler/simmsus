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

/**CPP File - Create the Configuration Object With the values from config.yaml input file */
// The methods are in alphabetical order

#include <header/config.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

using namespace std;
Configuration::Configuration(){ }

Configuration::~Configuration(){ }

Configuration::Configuration(char *configFile){

ifstream infile;
std::string string, recordName, value;

infile.open (configFile);

if (infile.is_open()){
	while(getline(infile,string)) // To get you all the lines.
	{
		istringstream sstream(string);
		sstream >> recordName >> value;
		char charStart = recordName.at(0);
		if (charStart != ' ' || charStart != '#') {
			if (recordName == "STE"){
				if(value == "TRUE"){
					this->ste = true;
				}else{
					this->ste = false;
				}
				continue;
			}
			if (recordName == "MP"){
				if(value == "TRUE"){
					this->mp = true;
				}else{
					this->mp = false;
				}
				continue;
			}
			if (recordName == "SDS"){
				if(value == "TRUE"){
					this->sds = true;
				}else{
					this->sds = false;
				}
				continue;
			}
			if (recordName == "SEDIMENTATION"){
				if(value == "TRUE"){
					this->sedimentation = true;
				}else{
					this->sedimentation = false;
				}
				continue;
			}
			if (recordName == "FLUCMODE"){
				if(value == "TRUE"){
					this->flucmode = true;
				}else{
					this->flucmode = false;
				}
				continue;
			}
			if (recordName == "BM"){
				if(value == "TRUE"){
					this->bm = true;
				}else{
					this->bm = false;
				}
				continue;
			}
			if (recordName == "ACCOUNTHI"){
				if(value == "TRUE"){
					this->accounthi = true;
				}else{
					this->accounthi = false;
				}
				continue;
			}
			if (recordName == "CONTINUEOSIM"){
				if(value == "TRUE"){
					this->continueosim = true;
				}else{
					this->continueosim = false;
				}
				continue;
			}
			if (recordName == "PI"){
				if(value == "TRUE"){
					this->pi = true;
				}else{
					this->pi = false;
				}
				continue;
			}
			if (recordName == "PMT"){
				if(value == "TRUE"){
					this->pmt = true;
				}else{
					this->pmt = false;
				}
				continue;
			}
			if (recordName == "PMF"){
				if(value == "TRUE"){
					this->pmf = true;
				}else{
					this->pmf = false;
				}
				continue;
			}
			if (recordName == "RECORDPOSFILE"){
				if(value == "TRUE"){
					this->recordposfile = true;
				}else{
					this->recordposfile = false;
				}
				continue;
			}
			if (recordName == "RECORDVELFILE"){
				if(value == "TRUE"){
					this->recordvelfile = true;
				}else{
					this->recordvelfile = false;
				}
				continue;
			}
			if (recordName == "RECORDDIPFILE"){
				if(value == "TRUE"){
					this->recorddipfile = true;
				}else{
					this->recorddipfile = false;
				}
				continue;
			}
			if (recordName == "OF"){
				if(value == "TRUE"){
					this->of = true;
				}else{
					this->of = false;
				}
				continue;
			}
			if (recordName == "TF"){
				if(value == "TRUE"){
					this->tf = true;
				}else{
					this->tf = false;
				}
				continue;
			}
			if (recordName == "NUMPART"){				
				this->numpart = stoi(value);			
				continue;
			}
			if (recordName == "NUMREAL"){				
				this->numreal = stoi(value);			
				continue;
			}
			if (recordName == "ORDRANDARR"){
				if(value == "TRUE"){
					this->ordrandarr = true;
				}else{
					this->ordrandarr = false;
				}
				continue;
			}
			if (recordName == "MIXMAGNONMAGPART"){
				if(value == "TRUE"){
					this->mixmagnonmagpart = true;
				}else{
					this->mixmagnonmagpart = false;
				}
				continue;
			}
			if (recordName == "INITIALSPHERAGGR"){
				if(value == "TRUE"){
					this->initialspheraggr = true;
				}else{
					this->initialspheraggr = false;
				}
				continue;
			}
			if (recordName == "ORDEREDDIPOLES"){
				if(value == "TRUE"){
					this->ordereddipoles = true;
				}else{
					this->ordereddipoles = false;
				}
				continue;
			}
			if (recordName == "MONOPOLIDISP"){
				if(value == "TRUE"){
					this->monopolidisp = true;
				}else{
					this->monopolidisp = false;
				}
				continue;
			}
			if (recordName == "VOLUMEFRACPART"){				
				this->volumefracpart = stod(value);			
				continue;
			}
			if (recordName == "BOXASPECTRATIO"){				
				this->boxaspectratio = stod(value);			
				continue;
			}
			if (recordName == "PERCENTNONMAGPART"){				
				this->percentnonmagpart = stod(value);			
				continue;
			}
			if (recordName == "APPLYEXTMAGFIELD"){
				if(value == "TRUE"){
					this->applyextmagfield = true;
				}else{
					this->applyextmagfield = false;
				}
				continue;
			}
			if (recordName == "POSEXTFIELD"){				
				this->posextfield = stod(value);			
				continue;
			}
			if (recordName == "OSCILLATORYFIELD"){
				if(value == "TRUE"){
					this->oscillatoryfield = true;
				}else{
					this->oscillatoryfield = false;
				}
				continue;
			}
			if (recordName == "ROTATINGFIELD"){
				if(value == "TRUE"){
					this->rotatingfield = true;
				}else{
					this->rotatingfield = false;
				}
				continue;
			}
			if (recordName == "NONLINEARDUFFFIELDEXCI"){
				if(value == "TRUE"){
					this->nonlineardufffieldexci = true;
				}else{
					this->nonlineardufffieldexci = false;
				}
				continue;
			}
			if (recordName == "DOUBLEFREQFIELDEXCI"){
				if(value == "TRUE"){
					this->doublefreqfieldexci = true;
				}else{
					this->doublefreqfieldexci = false;
				}
				continue;
			}
			if (recordName == "DYNAMICINCRFF"){
				if(value == "TRUE"){
					this->dynamicincrff = true;
				}else{
					this->dynamicincrff = false;
				}
				continue;
			}
			if (recordName == "FREQUENCY1MAGFIELD"){				
				this->frequency1magfield = stod(value);			
				continue;
			}
			if (recordName == "FREQUENCY2AGFIELD"){				
				this->frequency2magfield = stod(value);			
				continue;
			}
			if (recordName == "C1PARDUFFEXC"){				
				this->c1parduffexc = stod(value);			
				continue;
			}
			if (recordName == "C2PARDUFFEXC"){				
				this->c2parduffexc = stod(value);			
				continue;
			}
			if (recordName == "C3PARDUFFEXC"){				
				this->c3parduffexc = stod(value);			
				continue;
			}
			if (recordName == "C4PARDUFFEXC"){				
				this->c4parduffexc = stod(value);			
				continue;
			}
			if (recordName == "MAXFREQDYNINCR"){				
				this->maxfreqdynincr = stod(value);			
				continue;
			}
			if (recordName == "NUMBERINTDYNINCR"){				
				this->numberintdynincr = stoi(value);			
				continue;
			}
			if (recordName == "TURNONSHRATE"){
				if(value == "TRUE"){
					this->turnonshrate = true;
				}else{
					this->turnonshrate = false;
				}
				continue;
			}
			if (recordName == "OSCILLATORYSH"){
				if(value == "TRUE"){
					this->oscillatorysh = true;
				}else{
					this->oscillatorysh = false;
				}
				continue;
			}
			if (recordName == "DYNINCRSHRATE"){
				if(value == "TRUE"){
					this->dynincrshrate = true;
				}else{
					this->dynincrshrate = false;
				}
				continue;
			}
			if (recordName == "DIMENSIONLESSSHRATE"){				
				this->dimensionlessshrate = stod(value);			
				continue;
			}
			if (recordName == "FREQUENCYOSCILLSH"){				
				this->frequencyoscillsh = stod(value);			
				continue;
			}
			if (recordName == "LAMBDA"){				
				this->lambda = stod(value);			
				continue;
			}
			if (recordName == "ALPHA"){				
				this->alpha = stod(value);			
				continue;
			}
			if (recordName == "BROWNIANPECLETNUM"){				
				this->brownianpecletnum = stod(value);			
				continue;
			}
			if (recordName == "TRANSLATIONALSTOKESNUM"){				
				this->translationalstokesnum = stod(value);			
				continue;
			}
			if (recordName == "SIMULATIONTIME"){				
				this->simulationtime = stod(value);			
				continue;
			}
			if (recordName == "STEPSTORINGRESULTS"){				
				this->stepstoringresults = stoi(value);			
				continue;
			}
			if (recordName == "CONTINUEFITNUM"){				
				this->continuefitnum = stoi(value);			
				continue;
			}
			if (recordName == "STATANALYSIS"){
				if(value == "TRUE"){
					this->statanalysis = true;
				}else{
					this->statanalysis = false;
				}
				continue;
			}
			if (recordName == "CALCSTRUCTFACTOR"){
				if(value == "TRUE"){
					this->calcstructfactor = true;
				}else{
					this->calcstructfactor = false;
				}
				continue;
			}
			if (recordName == "PRINTLOCALMAPSPHI"){
				if(value == "TRUE"){
					this->printlocalmapsphi = true;
				}else{
					this->printlocalmapsphi = false;
				}
				continue;
			}
		}		
	}	
	infile.close();
}
}


void Configuration::setAccounthi(bool s_accounthi)
{
	accounthi = s_accounthi;
}

void Configuration::setAlpha(double s_alpha)
{
	alpha = s_alpha;
}

void Configuration::setApplyextmagfield(bool s_applyextmagfield)
{
	applyextmagfield = s_applyextmagfield;
}

void Configuration::setBm(bool s_bm)
{
	bm = s_bm;
}

void Configuration::setBoxaspectratio(int s_boxaspectratio)
{
	boxaspectratio = s_boxaspectratio;
}

void Configuration::setBrownianpecletnum(double s_brownianpecletnum)
{
	brownianpecletnum = s_brownianpecletnum;
}

void Configuration::setDynincrshrate(bool s_dynincrshrate)
{
	dynincrshrate = s_dynincrshrate;
}

void Configuration::setC1parduffexc(double s_c1parduffexc)
{
	c1parduffexc = s_c1parduffexc;
}

void Configuration::setC2parduffexc(double s_c2parduffexc)
{
	c2parduffexc = s_c2parduffexc;
}

void Configuration::setC3parduffexc(double s_c3parduffexc)
{
	c3parduffexc = s_c3parduffexc;
}

void Configuration::setC4parduffexc(double s_c4parduffexc)
{
	c4parduffexc = s_c4parduffexc;
}

void Configuration::setCalcstructfactor(bool s_calcstructfactor)
{
	calcstructfactor = s_calcstructfactor;
}

void Configuration::setContinuefitnum(int s_continuefitnum)
{
	continuefitnum = s_continuefitnum;
}

void Configuration::setContinueosim(bool s_continueosim)
{
	continueosim = s_continueosim;
}

void Configuration::setDimensionlessshrate(double s_dimensionlessshrate)
{
	dimensionlessshrate = s_dimensionlessshrate;
}

void Configuration::setDoublefreqfieldexci(bool s_doublefreqfieldexci)
{
	doublefreqfieldexci = s_doublefreqfieldexci;
}

void Configuration::setDynamicincrff(bool s_dynamicincrff)
{
	dynamicincrff = s_dynamicincrff;
}

void Configuration::setFlucmode(bool s_flucmode)
{
	flucmode = s_flucmode;
}

void Configuration::setFrequency1magfield(double s_frequency1magfield)
{
	frequency1magfield = s_frequency1magfield;
}

void Configuration::setFrequency2magfield(double s_frequency2magfield)
{
	frequency2magfield = s_frequency2magfield;
}

void Configuration::setFrequencyoscillsh(double s_frequencyoscillsh)
{
	frequencyoscillsh = s_frequencyoscillsh;
}

void Configuration::setInitialspheraggr(bool s_initialspheraggr)
{
	initialspheraggr = s_initialspheraggr;
}

void Configuration::setLambda(double s_lambda)
{
	lambda = s_lambda;
}

void Configuration::setMaxfreqdynincr(double s_maxfreqdynincr)
{
	maxfreqdynincr = s_maxfreqdynincr;
}

void Configuration::setMixmagnonmagpart(bool s_mixmagnonmagpart)
{
	mixmagnonmagpart = s_mixmagnonmagpart;
}

void Configuration::setMonopolidisp(bool s_monopolidisp)
{
	monopolidisp = s_monopolidisp;
}

void Configuration::setMp(bool s_mp)
{
	mp = s_mp;
}

void Configuration::setNonlineardufffieldexci(bool s_nonlineardufffieldexci)
{
	nonlineardufffieldexci = s_nonlineardufffieldexci;
}

void Configuration::setNumberintdynincr(int s_numberintdynincr)
{
	numberintdynincr = s_numberintdynincr;
}

void Configuration::setNumpart(int s_numpart)
{
	numpart = s_numpart;
}

void Configuration::setNumreal(int s_numreal)
{
	numreal = s_numreal;
}

void Configuration::setNumtimestep(double s_numtimestep)
{
	numtimestep = s_numtimestep;
}

void Configuration::setOf(bool s_of)
{
	of = s_of;
}

void Configuration::setOrdereddipoles(bool s_odereddipoles)
{
	ordereddipoles = s_odereddipoles;
}

void Configuration::setOrdrandarr(bool s_ordrandarr)
{
	ordrandarr = s_ordrandarr;
}

void Configuration::setOscillatoryfield(double s_oscillatoryfield)
{
	oscillatoryfield = s_oscillatoryfield;
}

void Configuration::setOscillatorysh(bool s_oscillatorysh)
{
	oscillatorysh = s_oscillatorysh;
}

void Configuration::setPercentnonmagpart(double s_percentnonmagpart)
{
	percentnonmagpart = s_percentnonmagpart;
}

void Configuration::setPi(bool s_pi)
{
	pi = s_pi;
}

void Configuration::setPmf(bool s_pmf)
{
	pmf = s_pmf;
}

void Configuration::setPmt(bool s_pmt)
{
	pmt = s_pmt;
}

void Configuration::setPosextfield(double s_posextfield)
{
	posextfield = s_posextfield;
}

void Configuration::setPrintlocalmapsphi(bool s_printlocalmapsphi)
{
	printlocalmapsphi = s_printlocalmapsphi;
}

void Configuration::setRecorddipfile(bool s_recorddipfile)
{
	recorddipfile = s_recorddipfile;
}

void Configuration::setRecordposfile(bool s_recordposfile)
{
	recordposfile = s_recordposfile;
}

void Configuration::setRecordvelfile(bool s_recordvelfile)
{
	recordvelfile = s_recordvelfile;
}

void Configuration::setRotatingfield(bool s_rotatingfield)
{
	rotatingfield = s_rotatingfield;
}

void Configuration::setSds(bool s_sds)
{
	sds = s_sds;
}

void Configuration::setSedimentation(bool s_sedimentation)
{
	sedimentation = s_sedimentation;
}

void Configuration::setSimulationtime(double s_simulationtime)
{
	simulationtime = s_simulationtime;
}

void Configuration::setStatanalysis(bool s_statanalysis)
{
	statanalysis = s_statanalysis;
}

void Configuration::setSte(bool s_ste)
{
	ste = s_ste;
}

void Configuration::setStepstoringresults(int s_stepstoringresults)
{
	stepstoringresults = s_stepstoringresults;
}

void Configuration::setTf(bool s_tf)
{
	tf = s_tf;
}

void Configuration::setTranslationalstokesnum(double s_translationalstokesnum)
{
	translationalstokesnum = s_translationalstokesnum;
}

void Configuration::setTurnonshrate(bool s_turnonshrate)
{
	turnonshrate = s_turnonshrate;
}

void Configuration::setVolumefracpart(double s_volumefracpart)
{
	volumefracpart = s_volumefracpart;
}