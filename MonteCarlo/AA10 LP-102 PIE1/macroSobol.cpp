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
#include "TCanvas.h"
#include "TApplication.h"
#include "TMath.h"
#include "THStack.h"


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
  int nAccidents = 5000;
  // cout << "Insert number of Accidents" << endl;
  // cin >> nAccidents;

  MyTF1 * water_spilled_PDF = new MyTF1("water_spilled_PDF", "1.0/(875.0-525.0)", 525, 875); //Uniform between +- 25% of 700;

  MyTF1 * activity_conc_H3_PDF  = new MyTF1("activity_conc_H3_PDF" , "gaus(0)", 2.4e+11, 3.5e+11); 
  activity_conc_H3_PDF->SetParameters(1, 2.94e+11, 0.05*2.94e+11);
  MyTF1 * activity_conc_C14_PDF = new MyTF1("activity_conc_C14_PDF", "gaus(0)", 4.93e+8, 6.93e+8); 
  activity_conc_C14_PDF->SetParameters(1, 5.93e+8, 0.05*5.93e+8);

  TFile windDist("windDist.root");
  MyTH1F * wind_PDF = (MyTH1F*)windDist.Get("hWind");

  //deposition velocity for H-3 (uniform) and C-14 (piecewise uniform), ref [7]
  MyTF1 * deposition_vel_H3_PDF  = new MyTF1("deposition_vel_H3_PDF",  "1.0/0.007", 0.001, 0.007);
  MyTF1 * deposition_vel_C14_PDF = new MyTF1("deposition_vel_C14_PDF", Deposition_vel_C14_PDF, 0.001, 0.013, 0);

  //inhalation and ingestion coefficient for H-3 retrieved from bilognormal distribution
  TFile copula("copula.root");
  MyTH2F * hBilognormal = (MyTH2F*)copula.Get("hBilognormal");

  //inhalation rate (m3/day) and food consuption lognormal distribution from [4]
  MyTF1 * inh_rate_PDF = new  MyTF1("inh_rate_PDF","ROOT::Math::lognormal_pdf(x,[0],[1])",0.1,33.86);
  inh_rate_PDF->SetParameters(2.8147,0.19); 

  //Lumped Translocation Parameter
  MyTF1 * lumped_par_PDF  = new MyTF1("lumped_par_PDF",  "1.0/(0.01-0.0067)", 0.0067, 0.01 );

  //stability class discrete PDF, ref [5] 
  //TUnuran stability_random;
  //stability_random.Init(stabilityDiscrDist(), "method=dgt"); //See definition of stabilityDiscrDist in Utility

  TNtuple *par = new TNtuple("par","par","water_spilled:activity_conc_H3:activity_conc_C14:release_duration:stability_class:sigma_z:wind_speed:deposition_vel_H3:deposition_vel_C14:lumped_par:inh_coeff:ing_coeff:inh_rate:dose");
  TNtuple *doses = new TNtuple("doses", "doses", "total_inhal_dose:total_cs_dose:total_gs_dose:total_ing_dose");
  
  ifstream in;
  in.open("Sobol_Sequence/Sobol5000.dat");
  string line;

  //Read Sobol numbers from Sobol.dat line by line
  //nAccidents must be equal to the number of lines in Sobol.dat
  for(int i=0; i< nAccidents; i++){  
    getline(in, line);
    std::istringstream iss(line);
    std::vector<double> sobolNumbers;
    std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(sobolNumbers));

    //Beginning of the Monte Carlo run LP102-PIE 1
    
      OuterSource outer_source;
  
      double water_spilled = water_spilled_PDF->GetRandom(sobolNumbers.at(0));
      outer_source.setWater_spilled(water_spilled);
      double release_duration = water_spilled*3600/1.8;
      outer_source.setRelease_duration(release_duration);//reference [1], water pumped at the rate of 1.80 Kg/hours 
  
      double activity_conc_H3 = activity_conc_H3_PDF->GetRandom(sobolNumbers.at(1));
      double activity_conc_C14 = activity_conc_C14_PDF->GetRandom(sobolNumbers.at(2));
      vector<double> activity_conc{ activity_conc_H3, activity_conc_C14};
      outer_source.setActivity_conc(activity_conc);
  
      outer_source.Go();

      Dispersion dispersion(outer_source);
  
      int stability_class = 3;//stability_random.SampleDiscr();
      dispersion.setStability_class(stability_class);
  
      //vertical dispersion coefficient at 300 m lognormal PDF, ref [6]
      MyTF1 sigma_zPDF = Sigma_zPDF(stability_class);//See definition of Sigma_zPDF in Utility
      double sigma_z = sigma_zPDF.GetRandom(sobolNumbers.at(3));
      dispersion.setSigma_z(sigma_z);
      
      double wind_speed = wind_PDF->GetRandom(sobolNumbers.at(4));//wind_random.Sample();
      if(wind_speed == 0) wind_speed = 0.5; //min wind speed
      dispersion.setWind_speed(wind_speed);
  
      vector<double> deposition_vel{deposition_vel_H3_PDF->GetRandom(sobolNumbers.at(5)), deposition_vel_C14_PDF->GetRandom(sobolNumbers.at(6))};
      //vector<double> deposition_vel{0, deposition_vel_C14_PDF->GetRandom(sobolNumbers.at(5))};
      dispersion.setDeposition_vel(deposition_vel);

      double stackHeight = 20;
      dispersion.setStack_height(stackHeight);
  
      dispersion.Go();
  
      Dose dose(dispersion);
  
      double lumped_par = lumped_par_PDF->GetRandom(sobolNumbers.at(7));
      dose.setLumped_translocation(lumped_par);
  
      double inh_coeff, ing_coeff, inh_rate, final_dose;
      hBilognormal->GetRandom2(inh_coeff, ing_coeff,sobolNumbers.at(8));
      dose.setInhal_coeff(inh_coeff*TMath::Power(10,-11));
      dose.setIng_coeff(ing_coeff*TMath::Power(10,-11));
  
      inh_rate = inh_rate_PDF->GetRandom(sobolNumbers.at(9));
      dose.setInhal_rate(inh_rate/(3600*24));
  
      dose.Go();
  
      final_dose = 1000*dose.getFinal_dose(); //dose in mSv

      double total_inhal_dose = 1000*dose.getInhal_dose();
      double total_cs_dose    = 1000*dose.getCS_dose();
      double total_gs_dose    = 1000*dose.getGS_dose();
      double total_ing_dose   = 1000*dose.getIng_dose();

    par->Fill(water_spilled, activity_conc_H3,activity_conc_C14,release_duration, stability_class, sigma_z, wind_speed,deposition_vel.at(0), deposition_vel.at(1),lumped_par, inh_coeff, ing_coeff,inh_rate, final_dose);
    doses->Fill(total_inhal_dose, total_cs_dose, total_gs_dose, total_ing_dose);    
    }//End of Monte Carlo

  double max_dose = 3;

  TH1F * hWater = new TH1F("hWater", "Water Spilled PDF", 400, 0, 875);
  TH1F * hActivity_conc_H3 = new TH1F("hActivity_conc_H3", "Activity conc. ^{3}H PDF", 100, 2.4e+11, 3.5e+11);
  TH1F * hActivity_conc_C14 = new TH1F("hActivity_conc_C14", "Activity conc. ^{14}C PDF", 100, 4.93e+8, 6.93e+8);
  TH1F * hTime = new TH1F("hTime", "Time PDF", 150, 0, 4000000);
  TH1F * hStability_class = new TH1F("hStability_class", "Stablity classes PDF", 7, 0, 7);
  TH1F * hSigma_z = new TH1F("hSigma_z", "#sigma_{z} PDF", 1000, 0, 100);
  TH1F * hWind = new TH1F("hWind", "Wind PDF", 150, 0, 15);
  TH1F * hDeposition_vel_H3 = new TH1F("hDeposition_vel_H3", "Deposition velocity ^{3}H PDF", 14, 0, 0.007);
  TH1F * hDeposition_vel_C14 = new TH1F("hDeposition_vel_C14", "Deposition velocity ^{14}C PDF", 24, 0.001, 0.013);
  TH1F * hLumped = new TH1F("hLumped", "Lumped Translocation Parameter PDF", 100, 0.0067, 0.01);
  TH1F * hInh = new TH1F("hInh", "Inh PDF", 150, 0, 5);
  TH1F * hIng = new TH1F("hIng", "Ing PDF", 150, 0, 25);
  TH1F * hInh_rate = new TH1F("hInh_rate", "Inh Rate PDF", 330, 0, 33.86);
  TH2F * WaterVsDose = new TH2F("WaterVsDose", "WaterVsDose", 400, 0, 875, 800, 0, max_dose);
  TH2F * WindVsDose = new TH2F("WindVsDose", "WindVsDose", 150, 0, 15, 800, 0, max_dose);
  TH2F * InhVsIng = new TH2F("InhVsIng", "InhVsIng", 150, 0, 5, 150, 0, 25);
  TH2F * Inh_rateVsDose = new TH2F("Inh_rateVsDose", "Inh_rateVsDose", 330, 0, 33.86, 800, 0, max_dose);
  TH2F * sigma_zVsDose = new TH2F("sigma_zVsDose", "Sigma_zVsDose", 1000, 0, 100, 800, 0, max_dose);
  TH2F * DepositionH3VsDose = new TH2F("DepositionH3VsDose", "DepositionH3VsDose", 140, 0, 0.007, 800, 0 , max_dose);

  TH1F * hDose = new TH1F("hDose", "Dose", 3000, 0, max_dose );
  TH1F * hDoseInhal = new TH1F("hDoseInhal", "Dose", 3000, 0, max_dose );
  TH1F * hDoseIng = new TH1F("hDoseIng", "Dose", 3000, 0, max_dose );
  TH1F * hDoseCS = new TH1F("hDoseCS", "Dose", 3000, 0, max_dose );
  TH1F * hDoseGS = new TH1F("hDoseGS", "Dose", 3000, 0, max_dose );
  

  par->Project("hWater", "water_spilled");
  par->Project("hTime", "release_duration");
  par->Project("hActivity_conc_H3", "activity_conc_H3");  
  par->Project("hActivity_conc_C14", "activity_conc_C14");
  par->Project("hStability_class", "stability_class");
  par->Project("hSigma_z", "sigma_z");
  par->Project("hWind", "wind_speed");
  par->Project("hDeposition_vel_H3", "deposition_vel_H3");
  par->Project("hDeposition_vel_C14", "deposition_vel_C14");
  par->Project("hLumped", "lumped_par");
  par->Project("hInh", "inh_coeff");
  par->Project("hIng", "ing_coeff");
  par->Project("hInh_rate", "inh_rate");
  par->Project("WaterVsDose", "dose:water_spilled");
  par->Project("WindVsDose", "dose:wind_speed");
  par->Project("InhVsIng", "ing_coeff:inh_coeff");
  par->Project("Inh_rateVsDose", "dose:inh_rate");
  par->Project("sigma_zVsDose", "dose:sigma_z");
  par->Project("DepositionH3VsDose", "dose:deposition_vel_H3");

  par->Project("hDose", "dose");

  doses->Project("hDoseInhal", "total_inhal_dose");
  doses->Project("hDoseIng", "total_ing_dose");
  doses->Project("hDoseCS", "total_cs_dose");
  doses->Project("hDoseGS", "total_gs_dose");

  THStack *hs = new THStack("hs","");

  hs->Add(hDoseInhal);
  hs->Add(hDoseIng);
  hs->Add(hDoseCS);
  hs->Add(hDoseGS);

  hs->Draw();


  double quantiles_mark[8] = {0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99};
  double quantiles[8];
  hDose->GetQuantiles(8, quantiles,quantiles_mark );
  
  cout<< "5% :\t"  << quantiles[0] << "\n";
  cout<< "10% :\t" << quantiles[1] << "\n";
  cout<< "25% :\t" << quantiles[2] << "\n";
  cout<< "50% :\t" << quantiles[3] << "\n";
  cout<< "75% :\t" << quantiles[4] << "\n";
  cout<< "90% :\t" << quantiles[5] << "\n";
  cout<< "95% :\t" << quantiles[6] << "\n";
  cout<< "99% :\t" << quantiles[7] << "\n";


  string fTitle = "Results/MC_results_";
  time_t _now =time(NULL );
  struct tm * now_time = localtime ( &_now );
  char* date = asctime(now_time);
  date[strlen(date) - 1] = 0; //remove unwnated \n from asctime output
  fTitle += (string)date;
  fTitle += ".root"; 
  TString rootfTitle = (TString)fTitle;
  TFile *f = new TFile(rootfTitle, "RECREATE");
	f->cd();
  hWater->Write();
  hActivity_conc_H3->Write();
  hActivity_conc_C14->Write();
  hTime->Write();
  hStability_class->Write();
  hSigma_z->Write();
  hWind->Write();
  hDeposition_vel_H3->Write();
  hDeposition_vel_C14->Write();
  hLumped->Write();
  hInh->Write();
  hIng->Write();
  hDose->Write();
  WaterVsDose->Write();
  WindVsDose->Write();
  InhVsIng->Write();
  Inh_rateVsDose->Write();
  sigma_zVsDose->Write();
  DepositionH3VsDose->Write();
  hs->Write();

	f->Close();

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