#ifndef UTILITY_H
#define UTILITY_H

using namespace std;

#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "Sobol.h"
// #include "TUnuran.h"
// #include "TUnuranEmpDist.h"
// #include "TUnuranDiscrDist.h"
#include <iostream>
#include <fstream>
#include <iterator>

void print_row(TVectorD, std::ofstream& );

TVectorD Mult(TVectorD , TMatrixD );

Double_t WaterPDF(Double_t *, Double_t *);

MyTF1 Sigma_zPDF(int);

const char* stability_conv(int);

// TUnuranEmpDist windEmpDist();
// TUnuranDiscrDist stabilityDiscrDist();

Double_t Deposition_vel_C14_PDF(Double_t *, Double_t *);

#endif


