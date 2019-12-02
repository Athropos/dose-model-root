#ifndef SOURCE_H
#define SOURCE_H

#include <string>

using namespace std;

double instant_release(double);
double dsgmf(double); 
double gauss(double);

class Source
{

    public:

    Source();
    Source(string, double, double(* )(double));
    Source(string, double );

    virtual ~Source();

    double Eval( double );
    double ConstEval( double );

    void setHalfLife( double ); 
    void setName( string );
    void setRelease_function( double(* )(double));
    void setRelease_duration( double );

    double getDecayConst() const;
    double getHalfLife() const;
    int getNuclei() const;
    string getName() const;
    double getRelease_duration() const;

    private:

    double release_duration;
    double (*function)(double);
    double pHalfLife; 
    int pNuclei;
    string pName;

};

#endif


