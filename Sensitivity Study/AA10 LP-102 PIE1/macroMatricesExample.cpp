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
  int nAccidents = 1000;

	MyTF2 *bigaussian = new MyTF2("bigaussian", "ROOT::Math::bigaussian_pdf(x, y, 1, 1, 0.8)", -5,5,-5,5);

	MyTF1* gaussian = new MyTF1("gaussian", "gaus(0)", -5, 5);
	gaussian->SetParameters(1, 0, 1);
	
	ofstream outA;
  outA.open("Example A");

	ofstream outB;
	outB.open("Example B");

	double x1,x2,x3,x4,x5,x6;

  ifstream in;
  in.open("Sobol_Sequence/SobolExample50k.dat");
  string line;

  //Read Sobol numbers from Sobol.dat line by line
  //nAccidents must be equal to the number of lines in Sobol.dat

  	//Read Sobol numbers from Sobol.dat line by line
  	//nAccidents must be equal to the number of lines in Sobol.dat
  	for (size_t i = 0; i < nAccidents; i++)
	  {
    	getline(in, line);
    	std::istringstream iss(line);
    	std::vector<double> sobolNumbers;
    	std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(sobolNumbers));

		std::size_t const half_size = sobolNumbers.size()/2;
    	std::vector<double> sobolNumbersB(sobolNumbers.begin(), sobolNumbers.begin() + half_size);
    	std::vector<double> sobolNumbersB(sobolNumbers.begin() + half_size, sobolNumbers.end());

		x1 = gaussian->GetRandom(sobolNumbersA.at(0));
		x2 = gaussian->GetRandom(sobolNumbersA.at(1));
		bigaussian->GetRandom2(x3, x4, sobolNumbersA.at(2));
		bigaussian->GetRandom2(x5, x6, sobolNumbersA.at(3));

		if(i!=0) {outA << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << endl;}

		x1 = gaussian->GetRandom(sobolNumbersB.at(0));
		x2 = gaussian->GetRandom(sobolNumbersB.at(1));
		bigaussian->GetRandom2(x3, x4, sobolNumbersB.at(2));
		bigaussian->GetRandom2(x5, x6, sobolNumbersB.at(3));

		if(i!=0) {outB << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << endl;}
	}
	
	outA.close();

	outB.close();

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