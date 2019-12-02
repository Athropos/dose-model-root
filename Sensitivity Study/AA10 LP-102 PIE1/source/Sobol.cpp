#include "../headers/Sobol.h"

using namespace std;

#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <iterator>

////////////////////////////////////////////////////////////////////////////////
/// Return a random number following this function shape
///
/// Author: Rene Brun   18/08/95
/// Modified by: Nicola Rizzi 05/11/2019


Double_t MyTF1::GetRandom(Double_t Sobol)
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
   Double_t r  = Sobol;
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

Double_t MyTH1F::GetRandom(Double_t Sobol) const
{
   if (fDimension > 1) {
      Error("GetRandom","Function only valid for 1-d histograms");
      return 0;
   }
   Int_t nbinsx = GetNbinsX();
   Double_t integral = 0;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbinsx+1] != fEntries) integral = ((TH1*)this)->ComputeIntegral(true);
      else  integral = fIntegral[nbinsx];
   } else {
      integral = ((TH1*)this)->ComputeIntegral(true);
   }
   if (integral == 0) return 0;
   // return a NaN in case some bins have negative content
   if (integral == TMath::QuietNaN() ) return TMath::QuietNaN();

   Double_t r1 = Sobol;
   Int_t ibin = TMath::BinarySearch(nbinsx,fIntegral,r1);
   Double_t x = GetBinLowEdge(ibin+1);
   if (r1 > fIntegral[ibin]) x +=
      GetBinWidth(ibin+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   return x;
}

void MyTH2F::GetRandom2(Double_t &x, Double_t &y, Double_t Sobol) 
{
   Int_t nbinsx = GetNbinsX();
   Int_t nbinsy = GetNbinsY();
   Int_t nbins  = nbinsx*nbinsy;
   Double_t integral;
   // compute integral checking that all bins have positive content (see ROOT-5894)
   if (fIntegral) {
      if (fIntegral[nbins+1] != fEntries) integral = ComputeIntegral(true);
      else integral = fIntegral[nbins];
   } else {
      integral = ComputeIntegral(true);
   }
   if (integral == 0 ) { x = 0; y = 0; return;}
   // case histogram has negative bins
   if (integral == TMath::QuietNaN() ) { x = TMath::QuietNaN(); y = TMath::QuietNaN(); return;}

   Double_t r1 = Sobol;
   Int_t ibin = TMath::BinarySearch(nbins,fIntegral,(Double_t) r1);
   Int_t biny = ibin/nbinsx;
   Int_t binx = ibin - nbinsx*biny;
   x = fXaxis.GetBinLowEdge(binx+1);
   if (r1 > fIntegral[ibin]) x +=
      fXaxis.GetBinWidth(binx+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*gRandom->Rndm();
}

void MyTF2::GetRandom2(Double_t &xrandom, Double_t &yrandom, Double_t Sobol)
{
   //  Check if integral array must be build
   Int_t i,j,cell;
   Double_t dx   = (fXmax-fXmin)/fNpx;
   Double_t dy   = (fYmax-fYmin)/fNpy;
   Int_t ncells = fNpx*fNpy;
   if (fIntegral.empty()) {
      fIntegral.resize(ncells+1);
      fIntegral[0] = 0;
      Double_t integ;
      Int_t intNegative = 0;
      cell = 0;
      for (j=0;j<fNpy;j++) {
         for (i=0;i<fNpx;i++) {
            integ = Integral(fXmin+i*dx,fXmin+i*dx+dx,fYmin+j*dy,fYmin+j*dy+dy);
            if (integ < 0) {intNegative++; integ = -integ;}
            fIntegral[cell+1] = fIntegral[cell] + integ;
            cell++;
         }
      }
      if (intNegative > 0) {
         Warning("GetRandom2","function:%s has %d negative values: abs assumed",GetName(),intNegative);
      }
      if (fIntegral[ncells] == 0) {
         Error("GetRandom2","Integral of function is zero");
         return;
      }
      for (i=1;i<=ncells;i++) {  // normalize integral to 1
         fIntegral[i] /= fIntegral[ncells];
      }
   }

// return random numbers
   Double_t r,ddx,ddy,dxint;
   r     = Sobol;
   cell  = TMath::BinarySearch(ncells,fIntegral.data(),r);
   dxint = fIntegral[cell+1] - fIntegral[cell];
   if (dxint > 0) ddx = dx*(r - fIntegral[cell])/dxint;
   else           ddx = 0;
   ddy = dy*gRandom->Rndm();
   j   = cell/fNpx;
   i   = cell%fNpx;
   xrandom = fXmin +dx*i +ddx;
   yrandom = fYmin +dy*j +ddy;
}
