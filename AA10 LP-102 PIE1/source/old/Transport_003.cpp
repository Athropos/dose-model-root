
#include "TRandom3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <time.h>
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TApplication.h"
#include "../../headers/Source.h"
#include "../../headers/Utility.h"



/*
v002: added checkpoints for emission fractions and division in sections
v003 Introduced Source Class and initial functions moved to Utility.cpp. Now the source is an object
and require an half life and a function (optional) at the declaration. Default: t1/2 = 1e9 and constant release from t=60s to t= 1400000s; 
the source is always normalized at 1 Bq.
loader.C introduced.
 */

/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] AA10_LP-102-PIE1 Publ UnMit - Transport_calc_a020_out_181014.xlsx
*/

void Transport(double half_life = 1*TMath::Power(10,9), double time_step = 1000, double final_time = 1600000){
//############################# INIT ##########################################

    clock_t tStart = clock();
	
	std::ofstream Results;
	Results.open("../Transport_results.txt");

    int steps_number = (int)final_time/time_step;
    double emission_fraction[] {0.01, 0.1, 0.5, 0.9, 0.99, 2};
    const int size = sizeof(emission_fraction)/sizeof(emission_fraction[0]);
    double emission_times[size];
    int time_mark = 0; 
    int room_number = 4;//3 internal + 1 external

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

//############################# TRANSPORT #####################################

    for (int step=0; step<steps_number+1; step++){

        time_elapsed(step)=step*time_step;

        //constant flows, reference [3]
        double f01= 0.03;
        double f12= 0.045167;
        double f23= 2.261167;

        T(0,1) = f01*time_step/V(0);
        T(1,2) = f12*time_step/V(1);
        T(2,3) = f23*time_step/V(2);

        //define the diagonal Tii as 1 - Di - sum_j(Tij)
        TMatrixDDiag diagT(T); diagT = 0;
        TVectorD sumOutflow(room_number); 
        for (int i=0; i<room_number; i++){
            for(int j=0; j<room_number; j++){

                if(j!=i) sumOutflow(i) += T(i,j);
            }
            
            diagT(i)= 1 - D(i) - sumOutflow(i);
        }
        if(step==0) T.Print();

    //############################# SOURCE ########################################
        Source S0;
        S(0) = S0.Eval(step*time_step)*time_step;

        if(step==0){
            TVectorD row = S;
            TMatrixDRow(N,step) = row;
        } else {
            TVectorD row = TMatrixDRow(N, step);
            TVectorD row_before = TMatrixDRow(N, step-1);
            
            row = S + Mult(row_before,T);
    
            TMatrixDRow(N,step) = row;
    
            if(time_mark<=size && row(room_number-1)>emission_fraction[time_mark]){
                emission_times[time_mark] = (step-1)*time_step;
                cout << emission_fraction[time_mark]*100 << "%  " <<  emission_times[time_mark] << " s" << endl; 
                time_mark++;
            }
    
            Results<< step << ' '; print_row(row, Results); 

        }

    }
    
    Results.close();
//############################# DECAY #########################################

    double decay_const = 0.693/half_life;

   for (int step=0; step<steps_number+1; step++){

        TVectorD row = TMatrixDRow(N,step); 
        row *= TMath::Exp(-decay_const*step*time_step); 
        TMatrixDRow(N,step) = row;
        
    } 

//############################# OUTPUTS #######################################
 
    //N.Print();

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

    auto legend = new TLegend(0.1,0.7,0.48,0.9);
    legend->SetHeader("Nuclei in Volume #"); 
    legend->AddEntry(result0,"N0","l");
    legend->AddEntry(result1,"N1","l");
    legend->AddEntry(result2,"N2","l");
    legend->AddEntry(result3,"N3 (emission)","l");
    legend->Draw();
    
    cout << "Elapsed Time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;
}

