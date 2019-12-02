#include "../headers/OuterSource.h"
#include "../headers/Output.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMath.h"
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>


using namespace std;

/* LIST OF REFERENCES
[1] Source term and dose to worker in accident analyses, ESS-0094187
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] AA10_LP-102-PIE1 Publ UnMit - Dose dominating nuclides_a001_181018.xlsx
[4] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
*/


OuterSource::OuterSource():             release_duration(1.4*TMath::Power(10,6)),
                                        water_spilled(0.),
                                        start_of_emission(0.),
                                        emission_duration(0.),
                                        emission_time(0.),
                                        names(),
                                        activities_conc(),
                                        activities_emitted(), 
                                        half_lives(){

    //cout<< "\n######### TRANSPORT #########\n";

    ifstream in; 
    in.open("Inner_Inventory"); // open inner inventory, ref [1]
    double activity, half_life;
    string name;
    int nlines = 0;

    while(in >> name >> activity >> half_life){
    	

        //only relevant nuclides for the dose, see [3] or Table 11 in [4]
        if(name=="H-3" || name =="C-14"){ 

    	names.push_back(name);
        activities_conc.push_back(activity);
        half_lives.push_back(half_life);
    	
        }

        if(activity != 0) nlines++;
    }
    ////cout << "\nNuclides in the inventory: " << nlines << "\n";
    in.close();
    }

void OuterSource::Go(){
    
    //remove("Outer_Inventory"); //remove existing outer inventory file
    //create text file
    //ofstream out;


    for (int i=0; i<names.size(); i++){

        //Define a source for each nuclide in the inventory
        Source s(names.at(i), half_lives.at(i));
        s.setRelease_duration(release_duration);
        //Calculate transport
        Transport t(s, 100, release_duration+release_duration/10.);
        t.Go();
        TMatrixD N = t.getMatrixN(); //retrieve 1 Bq normalized N matrix
        for (int step=0; step<N.GetNrows(); step++){
            TVectorD row = TMatrixDRow(N,step); 
            row *= activities_conc.at(i)*water_spilled;//actual activity in Bq at each step
            TMatrixDRow(N,step) = row;
        } 
        t.setMatrixN(N);
        vector<double> emission_fraction = t.getEmission_fraction();
        vector<double> emission_times = t.getEmission_times();
        //for definitions see [2], pag 29
        //H-3 is a long lived nuclide present in all PIE.
        if(names.at(i)=="H-3"){
            start_of_emission = emission_times.at(0);
            emission_duration = emission_times.at(3)-emission_times.at(0);
            emission_time = emission_times.at(2)-emission_times.at(0);
            for (size_t j = 0; j < emission_times.size(); j++)
            {
                //cout << emission_fraction.at(j)*100 << "\% \t" << emission_times.at(j) << " s\n"; 
            }

            //cout << "\nStart of emission: " << start_of_emission << " s\n";
            //cout << "Emission duration: " << emission_duration << " s\n";
            //cout << "Emission time: " << emission_time << " s\n"; 
        }

        //Output decay_out;
        //decay_out.OuterInventory(t, emission_times.at(4), out);
    TVectorD time_elapsed = t.getTimeElapsed();
    for(int j=time_elapsed.GetNrows()-1; j>=0; --j){


        if(time_elapsed(j) < emission_times.at(4)){
            
            activities_emitted.push_back(N(j+1, t.getRoom_number() -1));
            break;
        }
    }


    }

    

}
OuterSource::~OuterSource(){}

void OuterSource::setWater_spilled(double water){water_spilled = water;}
void OuterSource::setActivity_conc( vector<double> ac){activities_conc = ac;}
void OuterSource::setRelease_duration( double t){release_duration = t;}


vector<double> OuterSource::getActivityConc() const{return activities_conc;}
vector<double> OuterSource::getActivity_emitted() const{return activities_emitted;}
vector<string> OuterSource::getNames() const{return names;}
vector<double> OuterSource::getHalfLives() const{return half_lives;}

double OuterSource::getRelease_duration() const {return release_duration;}
double OuterSource::getStart_of_emission() const {return start_of_emission;}
double OuterSource::getEmission_duration() const {return emission_duration;}
double OuterSource::getEmission_time() const {return emission_time;}
