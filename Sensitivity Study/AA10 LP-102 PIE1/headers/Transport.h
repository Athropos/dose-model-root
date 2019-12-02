#ifndef TRANSPORT_H
#define TRANSPORT_H

#include "Source.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <vector>

using namespace std;

class Transport{
     
    public:
    Transport(Source, double ,double);
    virtual ~Transport();
    void Init();
    void Go();

    void setSource( Source ); 
    void setTime_step( double );
    void setFinal_time( double );
    void setMatrixN( TMatrixD );

    Source getSource() const;
    TMatrixD getMatrixN() const;
    TVectorD getTimeElapsed() const;
    TMatrixD getMatrixT() const;
    int getRoom_number() const;
    double getTime_step() const;
    vector<double> getEmission_fraction() ;
    vector<double> getEmission_times() ;
    

    private:
    int room_number;
    double pTime_step;
    double pFinal_time;
    int steps_number;
    vector<double> emission_fraction;
    vector<double> emission_times;


    Source src;
    TVectorD time_elapsed;
    TVectorD V;
    TMatrixD T;    
    TVectorD S;                 
    TVectorD D;  
    TMatrixD N;               
};

#endif


