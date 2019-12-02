#include "../headers/Source.h"
#include "TMath.h"
#include <iostream>

/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
*/

double instant_release(double x){
    if(x==0){ return 1;} //easy case
    else { return 0;}
}

double constant_release(double x){
    double tmin = 0;
    double tmax =1.4*TMath::Power(10,6); //reference [1], 700 kg of water pumped at the rate of 1.80 Kg/hours 
    if(x<tmax){return 1/(tmax - tmin);} //return the activity at time x
    else {return 0;}
    }

double gauss(double t) {return TMath::Gaus(t, 2000, 500, kTRUE);} //just for fun

Source::Source():pName(NULL), pHalfLife(1*TMath::Power(10,9)), function(constant_release){pNuclei = (int)pHalfLife/TMath::Log(2);}
Source::Source(string name, double HalfLife, double(*f)(double)): pName(name), pHalfLife(HalfLife), function(f){ pNuclei = (int)pHalfLife/TMath::Log(2);}
Source::Source(string name, double HalfLife): pName(name), pHalfLife(HalfLife), function(constant_release){pNuclei = (int)pHalfLife/TMath::Log(2);}
Source::~Source(){}

double Source::Eval(double t){return function(t);}

void Source::setHalfLife(double HalfLife ){
    pHalfLife = HalfLife;
    pNuclei = (int)pHalfLife/TMath::Log(2);
    } 
void Source::setName(string name){pName=name;}
void Source::setReleaseFunction(double(*f)(double)){function = f;}

double Source::getDecayConst() const {return 0.693/pHalfLife;}
double Source::getHalfLife() const {return pHalfLife;}
int Source::getNuclei() const {return pNuclei;}
string Source::getName() const {return pName;}