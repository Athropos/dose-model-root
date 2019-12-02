#ifndef UTILITY_H
#define UTILITY_H

using namespace std;

#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <iterator>


void print_row(TVectorD, std::ofstream& );

TVectorD Mult(TVectorD , TMatrixD );

TF1 Sigma_zPDF(int);

const char* stability_conv(int);


Double_t Deposition_vel_C14_PDF(Double_t *, Double_t *);




#endif


