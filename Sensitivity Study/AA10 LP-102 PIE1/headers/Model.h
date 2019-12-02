#ifndef MODEL_H
#define MODEL_H

#include "Dispersion.h"
#include <vector>

class Model{
    public:
    Model( vector<double> , bool = false);

    double getFinal_dose();

    private:
    vector<double> par;
    double final_dose;

};

#endif


