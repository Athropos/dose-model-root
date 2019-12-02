#include "headers/Transport.h"
#include "headers/Source.h"
#include "headers/Utility.h"
#include "headers/Output.h"
#include "headers/OuterSource.h"
#include "headers/Dispersion.h"
#include "headers/Dose.h"
#include "headers/Model.h"

#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>

#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH2F.h"
#include "TString.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TMath.h"


/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] J.D.Harrison, Uncertainties in dose coefficients for intakes of tritiated water and organically bound forms of tritium by members of the public, 2001
[4] USA EPA, Metabolically Derived Human Ventilation Rates: A revised approach based upon Oxygen consuption rates (Final Report, 2009)
[5] NKS-386, Addressing off-site consequence criteria using Level 3 PSA
[6] D.M. Hamby, The Gaussian Atmospheric Transport Model and its Sensitivity to the Joint Frequency Distributions and Parametric Variability, Health Phys., 2002
[7] ESS - Atmospheric modeling, ESS - 0051604
[8] J.Jacques, Sensitivity analysis in presence of model uncertainty and correlated inputs
[9] Andrea Saltelli, Paola Annoni - Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index

*/

//create A_B starting from A and B
vector<vector<double>> matrixA_B(vector<vector<double>> A,vector<vector<double>> B,size_t begin, size_t nPar){
  int nAccidents = A.size();
  //Substitute the i-th column of B with the one of A and estimate S(total)
  vector< vector<double> > A_B = A;
      for (size_t n = 0; n < nAccidents; n++)
    {
      for (size_t k = begin; k < begin+ nPar; k++)
      {
        A_B.at(n).at(k) = B.at(n).at(k);
      }
        
    }

  return A_B;
  
}

// S_i and S_Ti from [9]
double sensitivity_index(vector<double> doseA, vector< double > doseB, vector<double> doseA_B){

  int nAccidents = doseA.size();
  
    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Dj += doseB.at(n)*(doseA_B.at(n)-doseA.at(n));
    }
    return Dj/nAccidents;

}

double sensitivity_index(vector<double> doseA, vector< double > doseB, vector<double> doseA_B, vector<int> bootstrap){

  int nAccidents = bootstrap.size();
  
    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Dj += doseB.at(bootstrap.at(n))*(doseA_B.at(bootstrap.at(n))-doseA.at(bootstrap.at(n)));
    }
    return Dj/nAccidents;


}

double sensitivity_index_t(vector<double> doseA, vector< double > doseB, vector<double> doseA_B){

  int nAccidents = doseA.size();

    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Dj += (doseA.at(n)-doseA_B.at(n))*(doseA.at(n)-doseA_B.at(n));
    }
    Dj = Dj/(2*nAccidents);
    return Dj;

}

double sensitivity_index_t(vector<double> doseA, vector< double > doseB, vector<double> doseA_B,vector<int> bootstrap){

  int nAccidents = bootstrap.size();

    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Dj += (doseA.at(bootstrap.at(n))-doseA_B.at(bootstrap.at(n)))*(doseA.at(bootstrap.at(n))-doseA_B.at(bootstrap.at(n)));
    }
    Dj = Dj/(2*nAccidents);
    return Dj;

}

void StandaloneApplication(int argc, char** argv) {
  // ==>> here the ROOT macro is called
  clock_t tStart = clock();
  gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
  int nAccidents = 1000;
  vector<vector< double > > A, B;

  vector<double> doseA, doseB;

  ifstream inA, inB;
  inA.open("Example A");
  inB.open("Example B");
  string line;

  for(int i=0; i<nAccidents; i++){
    getline(inA, line);
    std::istringstream iss(line);
    std::vector<double> tmp;

    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));

    A.push_back(tmp);
  }
  for(int i=0; i<nAccidents; i++){
    getline(inB, line);
    
    std::istringstream iss(line);
    std::vector<double> tmp;

    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));

    B.push_back(tmp);
  }

  //numerical estimation of first order sensitivity indices, [8]
  
  int nParameters = A.at(1).size();
  cout << "nParameters = " << nParameters << "\tnAccidents = " << nAccidents << endl;
  double mean=0;
  double D=0;
  
  for (size_t i = 0; i < nAccidents; i++){

    Model dose_modelA(A.at(i), true);
    doseA.push_back(dose_modelA.getFinal_dose());
    mean += dose_modelA.getFinal_dose();
    D += dose_modelA.getFinal_dose()*dose_modelA.getFinal_dose();

    Model dose_modelB(B.at(i), true);
    doseB.push_back(dose_modelB.getFinal_dose());
  }
  
  mean = mean/nAccidents;
  D = D/nAccidents-mean*mean;

  cout << "mean: " << mean << "\tVar: " << D << "\n";

  int nIndices = 5;
  // doseA_B for each index!
  vector<vector<vector<double>>> A_Bs(nIndices);
  A_Bs.at(0).resize(nAccidents, vector<double>(nParameters));
  A_Bs.at(0) = matrixA_B(A, B, 0, 1); //matrix A_B for the first index;
  A_Bs.at(1).resize(nAccidents, vector<double>(nParameters));
  A_Bs.at(1) = matrixA_B(A, B, 1, 1);
  A_Bs.at(2).resize(nAccidents, vector<double>(nParameters));
  A_Bs.at(2) = matrixA_B(A, B, 0, 2);
  A_Bs.at(3).resize(nAccidents, vector<double>(nParameters));
  A_Bs.at(3) = matrixA_B(A, B, 2, 2);
  A_Bs.at(4).resize(nAccidents, vector<double>(nParameters));
  A_Bs.at(4) = matrixA_B(A, B, 4, 2);

  vector<vector<double>> doseA_B(nIndices);
  for(int i=0; i< nIndices; i++){
    for(int j=0; j<nAccidents; j++){

      Model dose_model(A_Bs.at(i).at(j), true);
      doseA_B.at(i).push_back(dose_model.getFinal_dose());
    }
  }



  //BCIs
  vector<TH1F> histVec; 
  vector<TH1F> hist_tVec; 

  for(int j=0; j<nIndices; j++){
    TString title = std::to_string(j);
    TString title_t = std::to_string(j) + "_t";
    TH1F hD(title, title, 1500, -1.5, 1.5);
    TH1F hD_t(title_t, title_t, 1500, -1.5, 1.5);
    histVec.push_back(hD);
    hist_tVec.push_back(hD_t);
  }

  int nSample=1000;
  int nBootstrap = 10000;

  for(int n=0; n<nBootstrap; n++){
    vector<int> bootstrap;

    for(int i=0; i<nSample; i++){
      double r =gRandom->Rndm() ;
      int index = r*nAccidents;
      bootstrap.push_back(index);
    }

    double mean_bootstrap = 0;
    double D_bootstrap = 0;

    for (size_t i = 0; i < nSample; i++){
    mean_bootstrap += doseA.at(bootstrap[i]);
    D_bootstrap += doseA.at(bootstrap[i])*doseA.at(bootstrap[i]);
    } 

    mean_bootstrap = mean_bootstrap/nSample;
    D_bootstrap = D_bootstrap/nSample-mean_bootstrap*mean_bootstrap;

    for(int j=0; j<nIndices; j++){
      double Dj = sensitivity_index(doseA, doseB, doseA_B.at(j), bootstrap);
      double Dj_t = sensitivity_index_t(doseA, doseB, doseA_B.at(j), bootstrap);
      histVec[j].Fill(1- Dj/D_bootstrap);
      hist_tVec[j].Fill(Dj_t/D_bootstrap);
    }

  }

  for(int j=0; j<nIndices; j++){
    double Dj = sensitivity_index(doseA, doseB, doseA_B.at(j));
    double Dj_t = sensitivity_index_t(doseA, doseB, doseA_B.at(j));
    double quantiles_mark[3] = {0.05, 0.5, 0.95};
    double quantiles[3];
    double quantiles_t[3];
    histVec[j].GetQuantiles(3, quantiles,quantiles_mark );
    hist_tVec[j].GetQuantiles(3, quantiles_t,quantiles_mark );
    cout << "Par "<< j <<"\n";
    cout << "S : " << Dj/D << " + " << quantiles[2]-quantiles[1] << " - " << quantiles[1]-quantiles[0]<< "\n";
    cout << "S_t : " << Dj_t/D <<  " + " << quantiles_t[2]-quantiles_t[1] << " - " << quantiles_t[1]-quantiles_t[0]<< "\n";
    cout << "\n";
  }

  TFile file("Results/exampleBCI(1000).root", "RECREATE");
  file.cd();
  for(int j=0; j<nIndices; j++){
    histVec[j].Write();
    hist_tVec[j].Write();
  }
  file.Close();

  cout << "\nElapsed Time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " s\n";

}
  // This is the standard "main" of C++ starting
  // a ROOT application
int main(int argc, char** argv) {

  TApplication app("ROOT Application", &argc, argv);
  StandaloneApplication(app.Argc(), app.Argv());
  app.Run();

  return 0;
}

 