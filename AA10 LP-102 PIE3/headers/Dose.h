#ifndef DOSE_H
#define DOSE_H

#include "Dispersion.h"
#include <vector>
#include "TMatrixD.h"

class Dose{
     
    public:
    Dose( Dispersion );
    virtual ~Dose();

    void Go();

    private:
    Dispersion dispersion;
    double emission_time;
    double  dry_soil_density;
    double  lumped_translocation;
    double  inhal_rate;
    double  food_consuption;
    double  time_delay;
    double  time_exposure;
    vector<string> names;
    vector<double> half_lives;
    vector<double> inhal_coeff;
    vector<double> ing_coeff;
    vector<double> cloudshine_coeff;
    vector<double> groundshine_coeff;
    vector<double> deposition_vel;
    vector<double> dose_inhal;
    vector<double> dose_ing;
    vector<double> dose_cloudshine;
    vector<double> dose_groundshine;
    double final_dose;


          
};

#endif


