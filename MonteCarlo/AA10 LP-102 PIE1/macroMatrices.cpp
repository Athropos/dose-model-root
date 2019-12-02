#include "headers/Transport.h"
#include "headers/Source.h"
#include "headers/Utility.h"
#include "headers/Output.h"
#include "headers/OuterSource.h"
#include "headers/Dispersion.h"
#include "headers/Dose.h"
#include "headers/Sobol.h"

#include <time.h>
#include <cassert>

#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TH2F.h"
#include "TString.h"
#include "TH1F.h"
// #include "TUnuran.h"
// #include "TUnuranEmpDist.h"
// #include "TUnuranDiscrDist.h"
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
*/


void StandaloneApplication(int argc, char** argv) {

  clock_t tStart = clock();
  gRandom = new TRandom3(0);
	gRandom->SetSeed(0);
  // ==>> here the ROOT macro is called
  int nAccidents = 3000;
//   cout << "Insert number of Accidents" << endl;
//   cin >> nAccidents;

  ofstream out;
  out.open("~/Documenti/Macros/Sensitivity Study/AA10 LP-102 PIE1/10000NoDep/Matrix A");
  
  MyTF1 * water_spilled_PDF = new MyTF1("water_spilled_PDF", "1.0/(875.0-525.0)", 525, 875); //Uniform between +- 25% of 700;

  MyTF1 * activity_conc_H3_PDF  = new MyTF1("activity_conc_H3_PDF" , "gaus(0)", 2.4e+11, 3.5e+11); 
  activity_conc_H3_PDF->SetParameters(1, 2.94e+11, 0.05*2.94e+11);
  MyTF1 * activity_conc_C14_PDF = new MyTF1("activity_conc_C14_PDF", "gaus(0)", 4.93e+8, 6.93e+8); 
  activity_conc_C14_PDF->SetParameters(1, 5.93e+8, 0.05*5.93e+8);

  TFile windDist("windDist.root");
  MyTH1F * wind_PDF = (MyTH1F*)windDist.Get("hWind");

  //deposition velocity for H-3 (uniform) and C-14 (piecewise uniform), ref [7]
  MyTF1 * deposition_vel_H3_PDF  = new MyTF1("deposition_vel_H3_PDF",  "1.0/0.007", 0, 0.007);
  MyTF1 * deposition_vel_C14_PDF = new MyTF1("deposition_vel_C14_PDF", Deposition_vel_C14_PDF, 0.001, 0.015, 0);

  //inhalation and ingestion coefficient for H-3 retrieved from bilognormal distribution
  TFile copula("copula.root");
  MyTH2F * hBilognormal = (MyTH2F*)copula.Get("hBilognormal");

  //inhalation rate (m3/day) and food consuption lognormal distribution from [4]
  MyTF1 * inh_rate_PDF = new  MyTF1("inh_rate_PDF","ROOT::Math::lognormal_pdf(x,[0],[1])",0.1,33.86);
  inh_rate_PDF->SetParameters(2.8147,0.19); 

  //Lumped Translocation Parameter
  MyTF1 * lumped_par_PDF  = new MyTF1("lumped_par_PDF",  "1.0/(0.01-0.0067)", 0.01, 0.0067);

  //stability class discrete PDF, ref [5] 
  //TUnuran stability_random;
  //stability_random.Init(stabilityDiscrDist(), "method=dgt"); //See definition of stabilityDiscrDist in Utility

  ifstream in;
  in.open("Sobol_Sequence/SobolforMatrices10kNoDep.dat");
  string line;

  //Read Sobol numbers from Sobol.dat line by line
  //nAccidents must be equal to the number of lines in Sobol.dat
for(int i=0; i< nAccidents; i++){  
    getline(in, line);  
    std::istringstream iss(line);
    std::vector<double> sobolNumbers;
    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(sobolNumbers));

    std::size_t const half_size = sobolNumbers.size()/2;
    std::vector<double> sobolNumbersA(sobolNumbers.begin(), sobolNumbers.begin() + half_size);
    std::vector<double> sobolNumbersB(sobolNumbers.begin() + half_size, sobolNumbers.end());
    
    //Beginning of the Monte Carlo run LP102-PIE 1
    
      OuterSource outer_source;
  
      double water_spilled = water_spilled_PDF->GetRandom(sobolNumbersA.at(0));
      outer_source.setWater_spilled(water_spilled);
      double release_duration = water_spilled*3600/1.8;
      outer_source.setRelease_duration(release_duration);//reference [1], water pumped at the rate of 1.80 Kg/hours 
  
      double activity_conc_H3 = activity_conc_H3_PDF->GetRandom(sobolNumbersA.at(1));
      double activity_conc_C14 = activity_conc_C14_PDF->GetRandom(sobolNumbersA.at(2));
      vector<double> activity_conc{ activity_conc_H3, activity_conc_C14};
      outer_source.setActivity_conc(activity_conc);

      outer_source.Go();

      Dispersion dispersion(outer_source);
  
      int stability_class = 3;//stability_random.SampleDiscr();
      dispersion.setStability_class(stability_class);
  
      //vertical dispersion coefficient at 300 m lognormal PDF, ref [6]
      MyTF1 sigma_zPDF = Sigma_zPDF(stability_class);//See definition of Sigma_zPDF in Utility
      double sigma_z = sigma_zPDF.GetRandom(sobolNumbersA.at(3));
      dispersion.setSigma_z(sigma_z);
      
      double wind_speed = wind_PDF->GetRandom(sobolNumbersA.at(4));//wind_random.Sample();
      if(wind_speed == 0) wind_speed = 0.5; //min wind speed
      dispersion.setWind_speed(wind_speed);

      //vector<double> deposition_vel{deposition_vel_H3_PDF->GetRandom(sobolNumbersA.at(5)), deposition_vel_C14_PDF->GetRandom(sobolNumbersA.at(6))};
      vector<double> deposition_vel{0, deposition_vel_C14_PDF->GetRandom(sobolNumbersA.at(5))};
      dispersion.setDeposition_vel(deposition_vel);
  
      dispersion.Go();
  
      Dose dose(dispersion);
  
      double lumped_par = lumped_par_PDF->GetRandom(sobolNumbersA.at(6));
      dose.setLumped_translocation(lumped_par);
  
      double inh_coeff, ing_coeff, inh_rate, final_dose;
      hBilognormal->GetRandom2(inh_coeff, ing_coeff,sobolNumbersA.at(7));
      dose.setInhal_coeff(inh_coeff*TMath::Power(10,-11));
      dose.setIng_coeff(ing_coeff*TMath::Power(10,-11));
  
      inh_rate = inh_rate_PDF->GetRandom(sobolNumbersA.at(8));
      dose.setInhal_rate(inh_rate/(3600*24));
  
      dose.Go();
  
      final_dose = 1000*dose.getFinal_dose(); //dose in mSv
      //cout << final_dose << endl;

      //if(i!=0) out << water_spilled << " " << activity_conc.at(0)<< " " << activity_conc.at(1)<< " " << sigma_z<< " " << wind_speed<< " " <<deposition_vel.at(0)<< " " << deposition_vel.at(1)<< " " <<lumped_par<< " " << inh_coeff<< " " << ing_coeff<< " " <<inh_rate<< " " <<"\n";
      if(i!=0) out << water_spilled << " " << activity_conc.at(0)<< " " << activity_conc.at(1)<< " " << sigma_z<< " " << wind_speed<< " " << deposition_vel.at(1)<< " " <<lumped_par<< " " << inh_coeff<< " " << ing_coeff<< " " <<inh_rate<< " " <<"\n";
    }//End of Monte Carlo

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