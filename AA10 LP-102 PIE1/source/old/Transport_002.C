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
#include <iterator>
#include <time.h>
#include "TStyle.h"
#include "TFile.h"
#include "TGraph.h"
#include "TBox.h"
#include "TMatrixD.h"
#include "TVectorD.h"

/*
v002: added checkpoints for emission fractions and division in sections
 */

/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] AA10_LP-102-PIE1 Publ UnMit - Transport_calc_a020_out_181014.xlsx
*/

void print_row(TVectorD row, std::ofstream& file){
    Int_t size = row.GetNoElements();

    for(Int_t i=0;i<size; i++){
        file << row(i) << ' ' ;
    }
    file << endl;
    }

TVectorD Mult(TVectorD vec, TMatrixD matrix){
    Int_t size = vec.GetNoElements();
    TVectorD out(size);

    for(Int_t i=0; i<size; i++){
        for(Int_t j=0; j<size; j++){
            out(i)+=vec(j)*matrix(j,i);
        }
    }
    
    return out;
}

void Transport(Int_t nAccidents=100){
//############################# INIT ##########################################

    clock_t tStart = clock();
	
	std::ofstream Results;
	Results.open("Transport_results.txt");

    Double_t time_step = 1;
    Double_t final_time = 40000;
    Int_t steps_number = (Int_t)final_time/time_step;
    Double_t emission_fraction[] {0.01, 0.1, 0.5, 0.9, 0.99, 2};
    const Int_t size = sizeof(emission_fraction)/sizeof(emission_fraction[0]);
    Double_t emission_times[size];
    Int_t time_mark = 0; 
    Int_t room_number = 4;//3 internal + 1 external

//############################# TRANSPORT #####################################
    TVectorD time_elapsed(steps_number+1);
    TVectorD V(room_number);// Init. volume vector, reference [3]
        V(0) = 42;
        V(1) = 812;
        V(2) = 40700;
        V(3) = 1; //EXTERNAL VOLUME NOT RELEVANT

    TMatrixD T(room_number, room_number);   // Init. transport matrix 
    TVectorD S(room_number);                // Init. source vector
    TVectorD D(room_number);                // Init. deposition factors vector

    TMatrixD N(steps_number+1, room_number);  // Init. of matrix with activities in each room (on rows) at different times (on columns) 

    for (Int_t step=0; step<steps_number+1; step++){

        time_elapsed(step)=step*time_step;

        //constant flows, reference [3]
        Double_t f01= 0.03;
        Double_t f12= 0.045167;
        Double_t f23= 2.261167;

        T(0,1) = f01*time_step/V(0);
        T(1,2) = f12*time_step/V(1);
        T(2,3) = f23*time_step/V(2);

        //define the diagonal Tii as 1 - Di - sum_j(Tij)
        TMatrixDDiag diagT(T); diagT = 0;
        TVectorD sumOutflow(room_number); 
        for (Int_t i=0; i<room_number; i++){
            for(Int_t j=0; j<room_number; j++){

                if(j!=i) sumOutflow(i) += T(i,j);
            }
            
            diagT(i)= 1 - D(i) - sumOutflow(i);
        }
        if(step==0) T.Print();

    //############################# SOURCE ########################################
        Double_t source_release = 1; //Activity in Bq
        S(0)= (step==0)*source_release; //instant release

        if(step==0){
            N(0,0) = S(0);
        } else {
            TVectorD row = TMatrixDRow(N,step); 
            TVectorD row_before = TMatrixDRow(N, step-1);
            
            row = Mult(row_before,T);
    
            TMatrixDRow(N,step) = row;
    
            if(time_mark<=size && row(room_number-1)>emission_fraction[time_mark]){
                emission_times[time_mark] = (step-1)*time_step;
                cout << emission_fraction[time_mark]*100 << "%  " <<  emission_times[time_mark] << " s" << endl; 
                time_mark++;
            }
    
            Results<< step << ' '; print_row(row, Results); 

        }

    }
//############################# DECAY #########################################

    Double_t half_life = 1000; //seconds
    Double_t decay_const = 0.693/half_life;

   for (Int_t step=0; step<steps_number+1; step++){

        TVectorD row = TMatrixDRow(N,step); 
        row *= TMath::Exp(-decay_const*step*time_step); 
        TMatrixDRow(N,step) = row;
        
    }

//############################# OUTPUTS #######################################
 
    N.Print();

    TVectorD N_0 = TMatrixDColumn(N,0);
    TVectorD N_1 = TMatrixDColumn(N,1);
    TVectorD N_2 = TMatrixDColumn(N,2);
    TVectorD N_3 = TMatrixDColumn(N,3);

    TGraph * result0 = new TGraph(time_elapsed, N_0);
        result0->GetYaxis()->SetTitle("Activity in volume normalized at 1 Bq");
        result0->GetXaxis()->SetTitle("Time elapsed (s)");
        result0->SetTitle("Results Transport LP-102 PIE1");
        result0->GetYaxis()->SetRangeUser(0.0000001, 1);
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
    
    TCanvas * c1 = new TCanvas("c1","Transport Result", 800, 800);
        c1->cd();
        c1->SetLogy();
        result0->Draw("AL");
        result1->Draw("L");
        result2->Draw("L");
        result3->Draw("L");
    
    cout << "Elapsed Time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

}
