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

    double getFinal_dose() const;
    double getInhal_dose() const;
    double getIng_dose() const;
    double getCS_dose() const;
    double getGS_dose() const;

    void setInhal_coeff( vector<double> );
    void setInhal_coeff( double );
    void setIng_coeff( vector<double> );
    void setIng_coeff( double );
    void setCloudshine_coeff( vector<double> );
    void setGroundshine_coeff( vector<double> );
    void setFood_consumption( double );
    void setInhal_rate( double );
    void setDeposition_vel( vector<double> );
    void setLumped_translocation( double );


    private:
    Dispersion dispersion;
    double  emission_time;
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
    vector<double> deposition_vel;
    vector<double> groundshine_coeff;
    vector<double> dose_inhal;
    vector<double> dose_ing;
    vector<double> dose_cloudshine;
    vector<double> dose_groundshine;

    double final_dose;
    double total_inhal_dose;
    double total_ing_dose;
    double total_cs_dose;
    double total_gs_dose;


          
};

#endif


