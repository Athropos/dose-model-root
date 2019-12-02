#ifndef DISPERSION_H
#define DISPERSION_H

#include "OuterSource.h"
#include <vector>
#include <string>
#include "TMatrixD.h"

class Dispersion{
     
    public:
    Dispersion(OuterSource );
    virtual ~Dispersion();

    OuterSource getSTouter();
    vector<double> getRelativeConcentration();
    vector<double> getActivityRealeased();

    void setEmission_time(double);
    void setTime_step(double);
    void setDistance (double);
    void setWind_speed(double);
    void setStability_class(int);
    void setStack_height(double);

    void Go();

    private:

    OuterSource STouter;

    double  emission_time;
    double  time_step;
    double  distance;
    double  wind_speed;
    int     stability_class;
    double  sigma_y;
    double  sigma_z;
    double  stack_height;
    double  mixing_layer;
    vector<double>  Q; 
    vector<double>  RC;
    vector<double> decay_depletion_factor;
    vector<string> names;
    vector<double> half_lives;
    TMatrixD Hosker;
    TMatrixD Briggs;
    double tau;
    double init_ver_broadening;
    double init_hor_broadening;
          
};

#endif


