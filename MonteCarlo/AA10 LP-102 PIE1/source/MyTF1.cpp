#include "../headers/MyTF1.h"

using namespace std;

#include "TMath.h"
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <iterator>

////////////////////////////////////////////////////////////////////////////////
/// Return a random number following this function shape
///
/// Author: Rene Brun   18/08/95
/// Modified by: Nicola Rizzi 05/11/2019


Double_t MyTF1::GetRandom()
{
   //check if integral array must be build
    if (fIntegral.size() == 0) {
      fIntegral.resize(fNpx + 1);
      fAlpha.resize(fNpx + 1);
      fBeta.resize(fNpx);
      fGamma.resize(fNpx);
      fIntegral[0] = 0;
      fAlpha[fNpx] = 0;
      Double_t integ;
      Int_t intNegative = 0;
      Int_t i;
      Bool_t logbin = kFALSE;
      Double_t dx;
      Double_t xmin = fXmin;
      Double_t xmax = fXmax;
      if (xmin > 0 && xmax / xmin > fNpx) {
         logbin =  kTRUE;
         fAlpha[fNpx] = 1;
         xmin = TMath::Log10(fXmin);
         xmax = TMath::Log10(fXmax);
      }
      dx = (xmax - xmin) / fNpx;

      Double_t *xx = new Double_t[fNpx + 1];
      for (i = 0; i < fNpx; i++) {
         xx[i] = xmin + i * dx;
      }
      xx[fNpx] = xmax;
      for (i = 0; i < fNpx; i++) {
         if (logbin) {
            integ = Integral(TMath::Power(10, xx[i]), TMath::Power(10, xx[i + 1]), 0.0);
         } else {
            integ = Integral(xx[i], xx[i + 1], 0.0);
         }
         if (integ < 0) {
            intNegative++;
            integ = -integ;
         }
         fIntegral[i + 1] = fIntegral[i] + integ;
      }
      if (intNegative > 0) {
         Warning("GetRandom", "function:%s has %d negative values: abs assumed", GetName(), intNegative);
      }
      if (fIntegral[fNpx] == 0) {
         delete [] xx;
         Error("GetRandom", "Integral of function is zero");
         return 0;
      }
      Double_t total = fIntegral[fNpx];
      for (i = 1; i <= fNpx; i++) { // normalize integral to 1
         fIntegral[i] /= total;
      }
      //the integral r for each bin is approximated by a parabola
      //  x = alpha + beta*r +gamma*r**2
      // compute the coefficients alpha, beta, gamma for each bin
      Double_t x0, r1, r2, r3;
      for (i = 0; i < fNpx; i++) {
         x0 = xx[i];
         r2 = fIntegral[i + 1] - fIntegral[i];
         if (logbin) r1 = Integral(TMath::Power(10, x0), TMath::Power(10, x0 + 0.5 * dx), 0.0) / total;
         else        r1 = Integral(x0, x0 + 0.5 * dx, 0.0) / total;
         r3 = 2 * r2 - 4 * r1;
         if (TMath::Abs(r3) > 1e-8) fGamma[i] = r3 / (dx * dx);
         else           fGamma[i] = 0;
         fBeta[i]  = r2 / dx - fGamma[i] * dx;
         fAlpha[i] = x0;
         fGamma[i] *= 2;
      }
      delete [] xx;
    }
   // return random number
   Double_t r  = 0.5;
   Int_t bin  = TMath::BinarySearch(fNpx, fIntegral.data(), r);
   Double_t rr = r - fIntegral[bin];

   Double_t yy;
   if (fGamma[bin] != 0)
      yy = (-fBeta[bin] + TMath::Sqrt(fBeta[bin] * fBeta[bin] + 2 * fGamma[bin] * rr)) / fGamma[bin];
   else
      yy = rr / fBeta[bin];
   Double_t x = fAlpha[bin] + yy;
   if (fAlpha[fNpx] > 0) return TMath::Power(10, x);
   return x;
}


