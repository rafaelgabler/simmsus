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

#include<config.hpp>
#include <string>
#include <stdio.h>
#include <yaml.h>
using namespace std;
Configuration::Configuration(){ }

Configuration::~Configuration(){ }

Configuration::Configuration(const char *configFile){

FILE *fh = fopen(configFile, "r");

yaml_parser_t parser;

/* Initialize parser */
if(!yaml_parser_initialize(&parser))
fputs("Failed to initialize parser!\n", stderr);
if(fh == NULL)
fputs("Failed to open file!\n", stderr);

/* Set input file */
yaml_parser_set_input_file(&parser, fh);

/* CODE HERE */

/* Cleanup */
yaml_parser_delete(&parser);
fclose(fh);


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

void Configuration::setBm(bool s_bm)
{
	bm = s_bm;
}

void Configuration::setBoxaspectratio(double s_boxaspectratio)
{
	boxaspectratio = s_boxaspectratio;
}

void Configuration::setBrownianpecletnum(double s_brownianpecletnum)
{
	brownianpecletnum = s_brownianpecletnum;
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

void Configuration::setContinuefitnum(double s_continuefitnum)
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

void Configuration::setC4parduffexc(double s_c4parduffexc)
{
	c4parduffexc = s_c4parduffexc;
}

void Configuration::setCalcstructfactor(bool s_calcstructfactor)
{
	calcstructfactor = s_calcstructfactor;
}

void Configuration::setContinuefitnum(double s_continuefitnum)
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

void Configuration::setDynincrshrate(bool s_dynincrshrate)
{
	dynincrshrate = s_dynincrshrate;
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

void Configuration::setNumberintdynincr(double s_numberintdynincr)
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

void Configuration::setStepstoringresults(double s_stepstoringresults)
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


