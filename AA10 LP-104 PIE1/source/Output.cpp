#include "../headers/Output.h"
#include "../headers/Source.h"
#include "../headers/Utility.h"
#include "../headers/Transport.h"
#include "TMath.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TString.h"   
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"

#include <chrono>
#include <ctime>
#include <string>
#include <iostream> 
#include <vector>

Output::Output(): file_root(){}

Output::~Output(){}

void Output::Root(Transport t, TFile *rootf){

    file_root = rootf;

    TMatrixD N = t.getMatrixN();
    TVectorD time_elapsed = t.getTimeElapsed();
    double half_life = t.getSource().getHalfLife();
    string name = t.getSource().getName();

    gROOT->SetBatch(kTRUE);//does not display canvas on screen
    TVectorD N_0 = TMatrixDColumn(N,0);
    TVectorD N_1 = TMatrixDColumn(N,1);
    TVectorD N_2 = TMatrixDColumn(N,2);
    TVectorD N_3 = TMatrixDColumn(N,3);

    TGraph * result0 = new TGraph(time_elapsed, N_0);
        TString gTitle; 
        gTitle.Form("Result %s (%E s)", name.c_str(), half_life);
        result0->SetTitle(gTitle);
        result0->GetYaxis()->SetTitle("Activity in volume");
        result0->GetXaxis()->SetTitle("Time elapsed (s)");
        result0->SetLineColor(kRed);
        result0->SetLineWidth(2);
    TGraph * result1 = new TGraph(time_elapsed, N_1);
        result1->SetLineColor(kYellow);
        result1->SetLineWidth(2);
    TGraph * result2 = new TGraph(time_elapsed, N_2);
        result2->SetLineColor(kGreen);
        result2->SetLineWidth(2);
    TGraph * result3 = new TGraph(time_elapsed, N_3); 
        result3->SetLineColor(kBlue);  
        result3->SetLineWidth(2);
    
    TCanvas * c2 = new TCanvas("c2","Decay Result", 800, 800);
        c2->cd();
        c2->SetLogy();
        result0->Draw("AL");
        result1->Draw("L");
        result2->Draw("L");
        result3->Draw("L");

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Nuclei in Volume #"); 
    legend->AddEntry(result0,"N0","l");
    legend->AddEntry(result1,"N1","l");
    legend->AddEntry(result2,"N2","l");
    legend->AddEntry(result3,"N3 (emission)","l");
    legend->Draw();
    c2->Update();

    file_root->cd();
    c2->Write((TString)name);
    

    delete c2;
    delete result0;
    delete result1;
    delete result2;
    delete result3;
    delete legend;
    
}

void Output::OuterInventory(Transport t , double abs_emi_time, ofstream& f){

        
    f.open("Outer_Inventory", std::ios_base::app);

    TMatrixD N = t.getMatrixN();
    TVectorD time_elapsed = t.getTimeElapsed();

    for(int i=time_elapsed.GetNrows()-1; i>=0; i--){


        if(time_elapsed(i) < abs_emi_time){
            
            f << t.getSource().getName() << "\t" <<  N(i+1, t.getRoom_number() -1) << "\t" << t.getSource().getHalfLife() <<"\n";
            break;
        }
    }

    f.close();

}




