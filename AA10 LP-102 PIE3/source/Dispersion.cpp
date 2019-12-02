#include "../headers/Transport.h"
#include "../headers/Source.h"
#include "../headers/Utility.h"
#include "../headers/Output.h"
#include "../headers/OuterSource.h"
#include "../headers/Dispersion.h"

#include "TMath.h"
#include <fstream>

/* LIST OF REFERENCES
[1] NRPB-R91, A Model for Short and Medium Range Dispersion of Radionuclides Released to the Atmosphere, R.H.Clarke
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] U. Bäverstam, "LENA 2003 Manual", 2003
[4] Hosker, Estimates of dry deposition and plume depletion over forests and grasslands (1974)
[5] Briggs, Diffusion estimates for small emission (1973)

Stability Class code: 0 - A, 1 - B, 2 - C, 3 - D, 4 - E, 5 - F
*/

/* DA INSERIRE: Dry and wet dep
                fattore di correzione per la durata */

Dispersion::Dispersion(OuterSource outer_source): STouter(outer_source),
                                            emission_time(0.),
                                            distance(0.),
                                            wind_speed(0.),
                                            stability_class(0),
                                            sigma_y(0.),
                                            sigma_z(0.),
                                            stack_height(0.),
                                            mixing_layer(0.),
                                            Q(),
                                            RC(),
                                            decay_depletion_factor(),
                                            names(),
                                            half_lives(),
                                            Hosker(6,4),
                                            Briggs(6,2),
                                            tau(0.),
                                            init_ver_broadening(0.),
                                            init_hor_broadening(0.){

    cout << "\n######### DISPERSION #########\n";
    emission_time = STouter.getEmission_time()/3600;  //from s to h                                       
    distance = 300; //m      
    init_ver_broadening = 20; //m Default on LENA
    init_hor_broadening = 20; //m Default on LENA  

    //Assume one set of conditions presented in [2]
    wind_speed = 0.84; //m/s
    stability_class = 5; 
    stack_height = 20; //m

    //Hosker paramaters for sigma_z, ref[3]
    //a, b, c and d (cols) as a function on Pasquill stability (rows)
    //from ref [4]
    Hosker(0,0) = 0.122;    Hosker(0,1) = 1.06;     Hosker(0,2) = 0.000538; Hosker(0,3) = 0.815;
    Hosker(1,0) = 0.13;     Hosker(1,1) = 0.95;     Hosker(1,2) = 0.000652; Hosker(1,3) = 0.75;
    Hosker(2,0) = 0.112;    Hosker(2,1) = 0.92;     Hosker(2,2) = 0.000905; Hosker(2,3) = 0.718;
    Hosker(3,0) = 0.098;    Hosker(3,1) = 0.889;    Hosker(3,2) = 0.00135;  Hosker(3,3) = 0.688;
    Hosker(4,0) = 0.0609;   Hosker(4,1) = 0.895;    Hosker(4,2) = 0.00196;  Hosker(4,3) = 0.684;
    Hosker(5,0) = 0.0638;   Hosker(5,1) = 0.783;    Hosker(5,2) = 0.00136;  Hosker(5,3) = 0.672;

    //Briggs paramaters for sigma_y, ref[3]
    //epsilon and gamma (cols) as a function on Pasquill stability (rows)
    //from ref [5]
    Briggs(0,0) = 0.22; Briggs(0,1) = 0.0001;
    Briggs(1,0) = 0.16; Briggs(1,1) = 0.0001;
    Briggs(2,0) = 0.11; Briggs(2,1) = 0.0001;
    Briggs(3,0) = 0.08; Briggs(3,1) = 0.0001;
    Briggs(4,0) = 0.06; Briggs(4,1) = 0.0001;
    Briggs(5,0) = 0.04; Briggs(5,1) = 0.0001;
    }

void Dispersion::Go(){

    switch (stability_class){
        case 0:
        mixing_layer = 1300; //m
        break;

        case 1:
        mixing_layer = 900; //m
        break;

        case 2:
        mixing_layer = 850; //m
        break;

        case 3:
        mixing_layer = 800; //m
        break;

        case 4:
        mixing_layer = 400; //m
        break;

        case 5:
        mixing_layer = 100; //m
        break;
    }

    ifstream in;
    in.open("Outer_Inventory");
    double activity, half_life;
    string name;
    int nlines = 0;

    while(in >> name >> activity >> half_life){
    	
        if(activity != 0.){
    	names.push_back(name);
        Q.push_back(activity);
        half_lives.push_back(half_life);
    	
        nlines++;
        }
    }
    in.close();

    cout << "Nuclides dispersed:\t" << nlines << "\n"; 


    //Evaluate depletion factor for each nuclide, see [3]
    for(int i=0; i<Q.size(); i++){
        decay_depletion_factor.push_back(TMath::Exp(-TMath::Log(0.5)/half_lives.at(i)*distance/wind_speed));
    }

    double a = Hosker(stability_class, 0);
    double b = Hosker(stability_class, 1);
    double c = Hosker(stability_class, 2);
    double d = Hosker(stability_class, 3);

    double epsilon = Briggs(stability_class, 0);
    double gamma = Briggs(stability_class, 1);

    //parameterization of the flutctuations in wind direction, see [1]
    tau = 0.065*TMath::Sqrt(7*emission_time/wind_speed);

    sigma_z = a*TMath::Power(distance,b)/(1+c*TMath::Power(distance,d));
    sigma_z = TMath::Sqrt(TMath::Power(sigma_z, 2) + TMath::Power(init_ver_broadening, 2));

    cout << "\nsigma z = " << sigma_z << " m\n";
    
    sigma_y = TMath::Power(epsilon*distance,2)/(1+gamma*distance) + TMath::Power(tau*distance,2);
    sigma_y = TMath::Sqrt(sigma_y + TMath::Power(init_hor_broadening, 2));

    cout << "sigma y = " << sigma_y << " m\n";

    //Plume reflections at ground level, see [1] or [3]

    double F = 2*TMath::Gaus(0, stack_height, sigma_z);         // 0 order reflection 
    F += 2*TMath::Gaus(2*mixing_layer, stack_height,sigma_z);   // 1st order reflection
    F += 2*TMath::Gaus(-2*mixing_layer, stack_height,sigma_z);  // 1st order reflection

    cout << "\nDilution from " << stack_height << " m, Weather: " << wind_speed << " m/s, " << stability_conv(stability_class) << " " << mixing_layer << " m\n";
    
    cout << "\nnuclide\t rel conc \t activity\n";

    //Gaussian Model for diluition at ground level on central axis for each nuclide
    for(int i=0; i<Q.size(); i++){
        double rel_conc = F*decay_depletion_factor.at(i)/(2*TMath::Pi()*sigma_z*sigma_y*wind_speed);
        RC.push_back(rel_conc);
        cout << names.at(i) << "\t" << rel_conc << "\t" << Q.at(i)*RC.at(i) << "\n";
    }
}

Dispersion::~Dispersion(){}

OuterSource Dispersion::getSTouter(){return STouter;}

vector<double> Dispersion::getRelativeConcentration(){return RC;}
vector<double> Dispersion::getActivityRealeased(){return Q;}

void Dispersion::setEmission_time(double emi_time){emission_time = emi_time;}
void Dispersion::setDistance (double d){distance = d;}
void Dispersion::setWind_speed(double w_speed){wind_speed = w_speed;}
void Dispersion::setStability_class(int stab){stability_class = stab;}
void Dispersion::setStack_height(double stack_h){stack_height = stack_h;}


