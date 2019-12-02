#include "TRandom3.h"
#include "TMath.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <fstream>
#include "TNtuple.h"
#include "TH2F.h"
#include <vector>
#include <iostream>
#include <vector>
#include "TStyle.h"
#include "TMatrixD.h"
#include "TGraph.h"
#include "TBox.h"
#include "TFile.h"
#include "TLatex.h"

double sigmazValue(double distance, int stability_class){

	double 		init_ver_broadening = 20; //m Default on LENA

	TMatrixD 	Hosker(6,4);
	//Hosker paramaters for sigma_z, ref[3]
    //a, b, c and d (cols) as a function on Pasquill stability (rows)
    //from ref [4]
    Hosker(0,0) = 0.122;    Hosker(0,1) = 1.06;     Hosker(0,2) = 0.000538; Hosker(0,3) = 0.815;
    Hosker(1,0) = 0.13;     Hosker(1,1) = 0.95;     Hosker(1,2) = 0.000652; Hosker(1,3) = 0.75;
    Hosker(2,0) = 0.112;    Hosker(2,1) = 0.92;     Hosker(2,2) = 0.000905; Hosker(2,3) = 0.718;
    Hosker(3,0) = 0.098;    Hosker(3,1) = 0.889;    Hosker(3,2) = 0.00135;  Hosker(3,3) = 0.688;
    Hosker(4,0) = 0.0609;   Hosker(4,1) = 0.895;    Hosker(4,2) = 0.00196;  Hosker(4,3) = 0.684;
    Hosker(5,0) = 0.0638;   Hosker(5,1) = 0.783;    Hosker(5,2) = 0.00136;  Hosker(5,3) = 0.672;

	double a = Hosker(stability_class, 0);
    double b = Hosker(stability_class, 1);
    double c = Hosker(stability_class, 2);
    double d = Hosker(stability_class, 3);

	double sigma_z = a*TMath::Power(distance,b)/(1+c*TMath::Power(distance,d));
    sigma_z = TMath::Sqrt(TMath::Power(sigma_z, 2) + TMath::Power(init_ver_broadening, 2));

	return sigma_z;



}


void Sigmaz(){

	
	int     	stability_class;
	double		distance;
	vector< vector<double> > x(6), y(6);


	for(stability_class =0; stability_class<6; ++stability_class){

		for(int i=0; i<1000; ++i){

			distance = i*10 + 100;
			x.at(stability_class).push_back(distance);

			double sigma_z = sigmazValue(distance, stability_class);
			y.at(stability_class).push_back(sigma_z);


			}
	}
	
	TGraph * g1 = new TGraph(x.at(0).size(), &x.at(0)[0], &y.at(0)[0]);
	TGraph * g2 = new TGraph(x.at(1).size(), &x.at(1)[0], &y.at(1)[0]);
	TGraph * g3 = new TGraph(x.at(2).size(), &x.at(2)[0], &y.at(2)[0]);
	TGraph * g4 = new TGraph(x.at(3).size(), &x.at(3)[0], &y.at(3)[0]);
	TGraph * g5 = new TGraph(x.at(4).size(), &x.at(4)[0], &y.at(4)[0]);
	TGraph * g6 = new TGraph(x.at(5).size(), &x.at(5)[0], &y.at(5)[0]);

	TCanvas * c1 = new TCanvas("c1", "f", 1000, 800);
	c1->cd();
	c1->SetLogx();
	c1->SetLogy();
	g1->Draw("AL");
	g1->SetTitle("");
	g1->GetYaxis()->SetRangeUser(9,1000);
	g1->GetYaxis()->SetTitle("Vertical dispersion coefficient [m]");
	g1->GetXaxis()->SetTitle("Downwind distance [m]");
	g1->Draw("L");
	g2->Draw("L");
	g3->Draw("L");
	g4->Draw("L");
	g5->Draw("L");
	g6->Draw("L");

	TLatex latex;
   	latex.SetTextSize(0.025);
	latex.SetTextFont(12);
   	latex.DrawLatex(2910,473,"A");
	latex.DrawLatex(2910,220,"B");
	latex.DrawLatex(2910,150,"C");
	latex.DrawLatex(2910,100,"D");
	latex.DrawLatex(2910,60,"E");
	latex.DrawLatex(2910,35,"F");
   	c1->Modified(); c1->Update(); 

	cout << "---------------------------- E5%\tE50%\tE95%\n";
	for(int i=0; i< 6; ++i){
		
		if(i!=0 && i !=5){
			
			cout << "Stability " << i << ", distance 300m: " << y.at(i+1).at(20)  << "\t" << y.at(i).at(20) <<"\t" << y.at(i-1).at(20) <<"\n";
		}
		if(i==0){
			
			cout << "Stability " << i << ", distance 300m: " << y.at(i+1).at(20)  << "\t" << y.at(i).at(20) <<"\t ------ \n";
		}
		if(i==5){
			
			cout << "Stability " << i << ", distance 300m: ------\t" << y.at(i).at(20)<<"\t" << y.at(i-1).at(20) <<"\n";
		}

	}


	
}
	

