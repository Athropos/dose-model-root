#include "../headers/Transport.h"
#include "../headers/Source.h"
#include "../headers/Utility.h"
#include "../headers/Output.h"
#include "../headers/OuterSource.h"
#include "../headers/Dispersion.h"
#include "../headers/Dose.h"
#include "../headers/Model.h"

#include <time.h>
#include <cassert>

#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH2F.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TMath.h"


/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] J.D.Harrison, Uncertainties in dose coefficients for intakes of tritiated water and organically bound forms of tritium by members of the public, 2001
[4] USA EPA, Metabolically Derived Human Ventilation Rates: A revised approach based upon Oxygen consuption rates (Final Report, 2009)
[5] NKS-386, Addressing off-site consequence criteria using Level 3 PSA
[6] D.M. Hamby, The Gaussian Atmospheric Transport Model and its Sensitivity to the Joint Frequency Distributions and Parametric Variability, Health Phys., 2002
[7] ESS - Atmospheric modeling, ESS - 0051604
*/


Model::Model(vector<double> parameters, bool example): par(parameters), final_dose(0.){

    if(example == false){

    double water_spilled,release_duration, activity_conc_H3,activity_conc_C14,deposition_vel_H3,deposition_vel_C14, sigma_z,wind_speed,lumped_par, inh_coeff, ing_coeff, inh_rate;
    vector<double> activity_conc, deposition_vel;
    water_spilled = par.at(0);
    OuterSource outer_source;
    activity_conc.push_back(par.at(1));
    activity_conc.push_back(par.at(2));
    outer_source.setActivity_conc(activity_conc);

    outer_source.setWater_spilled(water_spilled);
    release_duration = water_spilled*3600/1.8;
    outer_source.setRelease_duration(release_duration);//reference [1], water pumped at the rate of 1.80 Kg/hours 
    outer_source.Go();
    Dispersion dispersion(outer_source);
    int stability_class =3;
    dispersion.setStability_class(stability_class);

    double emission_height = 20;
    dispersion.setStack_height(emission_height);

    //vertical dispersion coefficient at 300 m lognormal PDF, ref [6]
    sigma_z =  par.at(3);
    dispersion.setSigma_z(sigma_z);
    
    wind_speed = par.at(4);
    if(wind_speed == 0) wind_speed = 0.5; //min wind speed
    dispersion.setWind_speed(wind_speed);
    // deposition velocity for H-3 and C-14, ref [7]
    deposition_vel.push_back(par.at(5));
    //deposition_vel.push_back(0);
    deposition_vel.push_back(par.at(6));
    dispersion.setDeposition_vel(deposition_vel);

    dispersion.Go();
    Dose dose(dispersion);
    lumped_par = par.at(7);
    dose.setLumped_translocation(lumped_par);
    inh_coeff = par.at(8)*TMath::Power(10,-11);
    ing_coeff = par.at(9)*TMath::Power(10,-11);
    dose.setInhal_coeff(inh_coeff);
    dose.setIng_coeff(ing_coeff);
    
    inh_rate = par.at(10)/(3600*24);

    dose.setInhal_rate(inh_rate);
    dose.Go();
    final_dose = 1000*dose.getFinal_dose(); //dose in mSv
    }

    if(example == true){ 
        //final_dose = 3*par[0]*par[1] + par[2]*par[3] + par[4]*par[5]; //a)
        
        //final_dose = par[0]*par[1] + par[2]*par[3] + par[4]*par[5]; //b)
        
        final_dose = par[0]*par[1] + par[2]*par[3] + 3*par[4]*par[5]; //c)
    }

}


double Model::getFinal_dose(){return final_dose;}


