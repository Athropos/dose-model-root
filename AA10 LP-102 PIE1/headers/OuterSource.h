#ifndef OUTERSOURCE_H
#define OUTERSOURCE_H

#include "../headers/Transport.h"

class OuterSource
{
    public:
    OuterSource( bool=true);
    virtual ~OuterSource();


    double getStart_of_emission() const;
    double getEmission_duration() const;
    double getEmission_time() const;

    private:

    double start_of_emission;
    double emission_duration;
    double emission_time;

};

#endif