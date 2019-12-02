#ifndef OUTERSOURCE_H
#define OUTERSOURCE_H

#include "../headers/Transport.h"

class OuterSource
{
    public:
    OuterSource();
    virtual ~OuterSource();

    void Go();
    void setWater_spilled(double);
    void setRelease_duration( double );
    void setActivity_conc( vector<double> );


    vector<double> getActivityConc() const;
    vector<double> getActivity_emitted() const;
    vector<string> getNames() const;
    vector<double> getHalfLives() const;
    double getRelease_duration() const;
    double getStart_of_emission() const;
    double getEmission_duration() const;
    double getEmission_time() const;

    private:

    vector<string> names;
    vector<double> activities_conc, half_lives, activities_emitted;

    double release_duration;
    double water_spilled;
    double start_of_emission;
    double emission_duration;
    double emission_time;

};

#endif