#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <fstream>
#include "TNtuple.h"
#include "TH2F.h"
#include <vector>
#include <time.h>
#include "TStyle.h"
#include "TGraph.h"
#include "TBox.h"
#include "TFile.h"
#include <time.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <vector>

void generatingExample(){


	MyTF2 *bigaussian = new MyTF2("bigaussian", "ROOT::Math::bigaussian_pdf(x, y, 1, 1, 0.8)", -5,5,-5,5);

	MyTF1* gaussian = new MyTF1("gaussian", "gaus(0)", -5, 5);
	gaussian->SetParameters(1, 0, 1);
	
	ofstream outA;
  	outA.open("Example A");

	ofstream outB;
	outB.open("Example B");

	double x1,x2,x3,x4,x5,x6;

	ifstream in;
  	in.open("Sobol_Sequence/SobolExample5000.dat");
  	string line;

  	//Read Sobol numbers from Sobol.dat line by line
  	//nAccidents must be equal to the number of lines in Sobol.dat
  	while(getline(in, line)){  
    	
    	std::istringstream iss(line);
    	std::vector<double> sobolNumbers;
    	std::copy(std::istream_iterator<double>(iss), std::istream_iterator<double>(), std::back_inserter(sobolNumbers));

		std::size_t const half_size = sobolNumbers.size()/2;
    	std::vector<double> sobolNumbersA(sobolNumbers.begin(), sobolNumbers.begin() + half_size);
    	std::vector<double> sobolNumbersB(sobolNumbers.begin() + half_size, sobolNumbers.end());

		x1 = gaussian->GetRandom(sobolNumbersA.at(0));
		x2 = gaussian->GetRandom(sobolNumbersA.at(1));
		bigaussian->GetRandom2(x3, x4, sobolNumbersA.at(2));
		bigaussian->GetRandom2(x5, x6, sobolNumbersA.at(3));

		outA << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << endl;

		x1 = gaussian->GetRandom(sobolNumbersB.at(0));
		x2 = gaussian->GetRandom(sobolNumbersB.at(1));
		bigaussian->GetRandom2(x3, x4, sobolNumbersB.at(2));
		bigaussian->GetRandom2(x5, x6, sobolNumbersB.at(3));

		outB << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << endl;
	}
	
	outA.close();

	outB.close();

	
}
	

