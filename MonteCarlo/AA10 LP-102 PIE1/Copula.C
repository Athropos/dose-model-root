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
#include <time.h>
#include "TStyle.h"
#include "TGraph.h"
#include "TBox.h"
#include "TFile.h"


void Copula(){

	TF2 *bigaussian = new TF2("bigaussian", "ROOT::Math::bigaussian_pdf(x, y, 1, 1, 0.5)", -5,5,-5,5);
	

	TCanvas * c1 = new TCanvas("c1", "f", 600, 600);
	c1->SetFrameLineColor(0);
	c1->cd();
	bigaussian->SetNpx(1000);
	bigaussian->SetNpy(1000);
	// bigaussian->Draw("colz");

	TH2F *hBigaussian = (TH2F*)(bigaussian->GetHistogram()->Clone("hBigaussian"));

	TCanvas * c2 = new TCanvas("c2", "h", 600, 600);
	c2->SetFrameLineColor(0);
	c2->cd();
	hBigaussian->Sumw2(0);
	hBigaussian->Draw("COLA");

	TH1F *hBigaussian_px = (TH1F*)hBigaussian->ProjectionX("hBigaussian_px");

	TCanvas * c3 = new TCanvas("c3", "h_px", 600, 600);
	c3->SetFrameLineColor(0);
	c3->cd();
	hBigaussian_px->Sumw2(0);
	hBigaussian_px->Draw();

	TH1F *hBigaussian_py = (TH1F*)hBigaussian->ProjectionY("hBigaussian_py");
	TCanvas * c4 = new TCanvas("c4", "h_py", 600, 600);
	c2->SetFrameLineColor(0);
	c4->cd();
	hBigaussian_py->Sumw2(0);
	hBigaussian_py->Draw();

	TH1F * hLognormal_px = new TH1F("hLognormal_px", "hLognormal_px", 500, 0, 5);
	TH1F * hLognormal_py = new TH1F("hLognormal_py", "hLognormal_py", 2500, 0, 25);

	TH1F * hUniform_px = new TH1F("hUniform_px", "hUniform_px", 100, 0, 1);
	TH1F * hUniform_py = new TH1F("hUniform_py", "hUniform_py", 100, 0, 1);
	TH2F * hBilognormal = new TH2F("hBilognormal", "hBilognormal", 500, 0,5,2500, 0,25);
	TH2F * hCopula = new TH2F("hCopula", "hCopula", 100, 0, 1, 100, 0 , 1);

	for(int i=0; i<10000000; ++i){
		double gaussian_px, gaussian_py;

		hBigaussian->GetRandom2(gaussian_px, gaussian_py);

		double uniform_px = ROOT::Math::normal_cdf(gaussian_px);
		hUniform_px->Fill(uniform_px);
		double lognormal_px = ROOT::Math::lognormal_quantile(uniform_px,0.83,0.21);
		hLognormal_px->Fill(lognormal_px);

		double uniform_py = ROOT::Math::normal_cdf(gaussian_py);
		hUniform_py->Fill(uniform_py);
		double lognormal_py = ROOT::Math::lognormal_quantile(uniform_py,1.72,0.42);
		hLognormal_py->Fill(lognormal_py);

		hCopula->Fill(uniform_px, uniform_py);
		hBilognormal->Fill(lognormal_px,lognormal_py);
	}

	TCanvas * c8 = new TCanvas("c8", "c_uniform_px", 600, 600);
	c2->SetFrameLineColor(0);
	c8->cd();
	hUniform_px->Draw();

	TCanvas * c9 = new TCanvas("c9", "c_uniform_py", 600, 600);
	c2->SetFrameLineColor(0);
	c9->cd();
	hUniform_py->Draw();

	TCanvas * c10 = new TCanvas("c10", "c_copula", 600, 600);
	c2->SetFrameLineColor(0);
	c10->cd();
	hCopula->Draw("COLA");

	TCanvas * c5 = new TCanvas("c5", "h_lognormal_px", 600, 600);
	c2->SetFrameLineColor(0);
	c5->cd();
	hLognormal_px->Draw();

	TCanvas * c6 = new TCanvas("c6", "h_lognormal_py", 600, 600);
	c2->SetFrameLineColor(0);
	c6->cd();
	hLognormal_py->Draw();

	TCanvas * c7 = new TCanvas("c7", "hBilognormal", 600, 600);
	c2->SetFrameLineColor(0);
	c7->cd();
	hBilognormal->Draw("COLA");

	TFile *f = new TFile("copula.root", "RECREATE");
	f->cd();
	hBigaussian->Write();
	hBigaussian_px->Write();
	hBigaussian_py->Write();
	hUniform_px->Write();
	hUniform_py->Write();
	hCopula->Write();
	hLognormal_px->Write();
	hLognormal_py->Write();
	hBilognormal->Write();
	f->Close();


	
}
	

