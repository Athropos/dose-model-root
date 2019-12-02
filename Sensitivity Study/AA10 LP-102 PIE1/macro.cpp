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

// S_i and S_Ti from [9]
double sensitivity_index(vector< vector<double>> A, vector< vector<double>> B, size_t begin, size_t nPar){

  int nAccidents = A.size();
  
  //Substitute the i-th column of A with the one of B and estimate S
  vector< vector<double> > A_B = A;
      for (size_t n = 0; n < nAccidents; n++)
    {
      for (size_t k = begin; k < begin+nPar; k++)
      {
        A_B.at(n).at(k) = B.at(n).at(k);
      }
        
    }

    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Model f_B(B.at(n));
      Model f_A_B(A_B.at(n));
      //Dj += f_B.getFinal_dose()*f_A_B.getFinal_dose();
      Dj += (f_B.getFinal_dose()-f_A_B.getFinal_dose())*(f_B.getFinal_dose()-f_A_B.getFinal_dose());
    }
    return Dj/(2*nAccidents);

}

double sensitivity_index_t(vector< vector<double>> A, vector< vector<double>> B, size_t begin, size_t nPar){

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

    double Dj =0;
    for(size_t n=0; n< nAccidents; n++){
      Model f_A(A.at(n));
      Model f_A_B(A_B.at(n));
      Dj += (f_A.getFinal_dose()-f_A_B.getFinal_dose())*(f_A.getFinal_dose()-f_A_B.getFinal_dose());
    }
    Dj = Dj/(2*nAccidents);
    return Dj;

}
void StandaloneApplication(int argc, char** argv) {

  clock_t tStart = clock();
  // gRandom = new TRandom3(0);
	// gRandom->SetSeed(0);
  // ==>> here the ROOT macro is called
  vector<vector< double > > A, B;

  ifstream inA, inB;
  inA.open("1000/Matrix A");
  inB.open("1000/Matrix B");
  string line;

  while(getline(inA, line)){
    
    std::istringstream iss(line);
    std::vector<double> tmp;

    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));

    A.push_back(tmp);
  }

  while(getline(inB, line)){
    
    std::istringstream iss(line);
    std::vector<double> tmp;

    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(tmp));

    B.push_back(tmp);
  }

  //numerical estimation of first order sensitivity indices, [8]
  int nAccidents = A.size();
  int nParameters = A.at(1).size();
  cout << "nParameters = " << nParameters << "\tnAccidents = " << nAccidents << endl;
  double mean=0; //.303675 ; //=0.303596; // Matrix A
  double D=0; //.0640519; //=0.0634884; // Matrix A
  
  for (size_t i = 0; i < nAccidents; i++){

    Model dose_model(A.at(i));
    mean += dose_model.getFinal_dose();
    D += dose_model.getFinal_dose()*dose_model.getFinal_dose();
  }
  
  mean = mean/nAccidents;
  D = D/nAccidents-mean*mean;

  cout << "mean: " << mean << "\tVar: " << D << "\n";

  vector<double> S;

  double Dj = sensitivity_index_t(A, B, 0, 1);
  cout << "Par water : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 1, 1);
  cout << "Par activity conc H3 : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 2, 1);
  cout << "Par activity conc C14 : " <<Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 3, 1);
  cout << "Par simga_z : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 4, 1);
  cout << "Par wind_speed : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 5, 1);
  cout << "Par dep vel H3 : " <<Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 6, 1);
  cout << "Par dep vel C14 : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 7, 1);
  cout << "Par lumped par : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 8, 2);
  cout << "Par inh-ing coeff : " << Dj/D << "\n";

  S.push_back(Dj/D);

  Dj = sensitivity_index_t(A, B, 10, 1);
  cout << "Par inh rate : " << Dj/D << "\n";

  S.push_back(Dj/D);
  

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

 