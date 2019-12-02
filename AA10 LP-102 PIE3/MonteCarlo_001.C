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
#include "TFile.h"
#include "TGraph.h"
#include "TBox.h"
#include "TMatrixD.h"
#include "TVectorD.h"

/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
*/

//Triangular Probability Density Function
Double_t TriangularPDF(Double_t *x, Double_t *par){
    //Parameters
    Double_t xx= x[0];
    Double_t min = par[0];
    Double_t max = par[1];
    Double_t mode = par[2];

    //Cases definition
    if(xx < min) return 0;
    if(xx >= min && xx < mode) return 2*(xx - min)/((max-min)*(mode-min));
    if(xx >= mode && xx <= max) return 2*(max - xx)/((max-min)*(max - mode));
    if(xx > max) return 0;
    return 0;
}

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

void MonteCarlo(Int_t nAccidents=100){

    ////#############################   INIT   #####################################
     clock_t tStart = clock();
 
	 std::ofstream Results;
	 Results.open("Results.txt");
 
     Double_t time_step = 100;
     Double_t final_time = 400000;
     Int_t steps_number = (Int_t)final_time/time_step;
     Int_t room_number = 4;//3 internal + 1 external

     //Triangular Distribution

        Double_t a = 500;  //min
        Double_t b = 1640; //max equal to MAR, see [1]
        Double_t c = 700; //mode 

        TF1 * triangularPDF = new TF1("triangularPDF", TriangularPDF, a, b, 3);
        triangularPDF->SetParameters(a, b, c);  
    
     TH1F * hTriangular = new TH1F("hTriangular", "Triangular PDF", 800, 0, b);
     TString title_emitted = Form("Normalized Activity Emitted in %.0f s", final_time);
     TH1F * hE0 = new TH1F("hE0", title_emitted, 1000, 0, 1);

    //Beginning of the Monte Carlo run LP102-PIE 1
    for(Int_t i=0; i< nAccidents; i++){
     //########################## RANDOM VARIABLES ################################# 
        Double_t activity_conc = 5.13*TMath::Power(10,11); //Bq/kg reference [1]

        Double_t water_spilled = triangularPDF->GetRandom();
        hTriangular->Fill(water_spilled);
    
        Double_t max_activity_released = 1640 * activity_conc;
        Double_t activity_released = water_spilled * activity_conc;

     //############################# TRANSPORT #####################################
        TVectorD time_elapsed(steps_number+1);
        TVectorD V(room_number);                // Init. volume vector
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

            //constant flows
            Double_t f01= 0.03;
            Double_t f12= 0.0452;
            Double_t f23= 2.261;

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

            //############################# SOURCE #####################################
            S(0)=(step==0)*activity_released; //instant release

            if(step==0){
                N(0,0) = S(0);
            } else {
            TVectorD row = TMatrixDRow(N,step); 
            TVectorD row_before = TMatrixDRow(N, step-1);

            row = Mult(row_before,T);

            TMatrixDRow(N,step) = row;
            }

        }   
     //#############################   DECAY   #####################################

        Double_t half_life = 1000000000; //seconds
        Double_t decay_const = 0.693/half_life;

        for (Int_t step=0; step<steps_number+1; step++){

            TVectorD row = TMatrixDRow(N,step); 
            row *= TMath::Exp(-decay_const*step*time_step); 
            TMatrixDRow(N,step) = row;
        
        }

        hE0->Fill(N(steps_number,3)/max_activity_released);

    }
    //#############################   OUTPUTS   ####################################
        TCanvas * c1 = new TCanvas("c1", "Water Spilled", 600, 600);
        c1->cd();
        hTriangular->Draw();

        TCanvas * c2 = new TCanvas("c2", "Activity Emitted" , 600, 600);
        c2->cd();
        hE0->Draw();
    
        cout << "Elapsed Time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << endl;

}
