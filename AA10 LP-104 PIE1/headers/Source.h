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

    void setHalfLife( double ); 
    void setName( string );
    void setReleaseFunction( double(* )(double));

    double getDecayConst() const;
    double getHalfLife() const;
    int getNuclei() const;
    string getName() const;

    private:

    double (*function)(double);
    double pHalfLife; 
    int pNuclei;
    string pName;

};

#endif


