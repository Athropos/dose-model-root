#ifndef DISPERSION_H
#define DISPERSION_H

#include "OuterSource.h"
#include <vector>
#include <string>
#include "TMatrixD.h"

Double_t F_func(Double_t *, Double_t *);
Double_t sigma_z_func(Double_t *, Double_t *);


class Dispersion{
     
    public:
    Dispersion(OuterSource);
    virtual ~Dispersion();

    OuterSource getSTouter();

    vector<double> getDeposition_vel();
    vector<double> getWash_out();
    vector<double> getRelative_conc();
    vector<double> getRelative_ground_conc();
    vector<double> getActivity_realeased();


    void setEmission_time(double);
    void setTime_step(double);
    void setDistance (double);
    void setWind_speed(double);
    void setStability_class(int);
    void setStack_height(double);
    void setSigma_z(double);
    void setDeposition_vel( vector<double> );

    void Go();

    private:

    OuterSource STouter;

     double  emission_time;
    double  time_step;
    double  distance;
    double  wind_speed;
    double  rain_intensity;
    int     stability_class;
    double  sigma_y;
    double  sigma_z;
    double  stack_height;
    double  mixing_layer;
    vector<double>  Q; 
    vector<double>  RC;
    vector<double>  RGC;
    vector<double> decay_depletion_factor;
    vector<double> deposition_vel;
    vector<double> wash_out_A;
    vector<double> wash_out_B;
    vector<double> wash_out;
    vector<double> dry_dep_depletion_factor;
    vector<string> names;
    vector<double> half_lives;
    TMatrixD Hosker;
    TMatrixD Briggs;
    double tau;
    double init_ver_broadening;
    double init_hor_broadening;
          
};

#endif


