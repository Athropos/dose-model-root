#ifndef UTILITY_H
#define UTILITY_H

using namespace std;

#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>
#include <fstream>
#include <iterator>


void print_row(TVectorD, std::ofstream& );

TVectorD Mult(TVectorD , TMatrixD );

const char* stability_conv(int);

#endif


