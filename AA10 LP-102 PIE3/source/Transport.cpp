#include "../headers/Transport.h"
#include "TMath.h"
#include <time.h>
#include "TApplication.h"
#include "../headers/Source.h"
#include "../headers/Utility.h"


/* LIST OF REFERENCES
[1] AA10 - Accident Analysis Report - Public: Water leakage in Monolith Vessel, ESS-0052199
[2] ESS - Activity transport and dose calculation models and tools used in safety analyses at ESS, ESS-0092033
[3] AA10_LP-102-PIE1 Publ UnMit - Transport_calc_a020_out_181014.xlsx
*/

Transport::Transport(Source S0, double time_step = 100, double final_time = 100000):pTime_step(time_step),
                                                                pFinal_time(final_time),
                                                                room_number(4),
                                                                steps_number ((int)final_time/time_step),
                                                                emission_fraction(),
                                                                emission_times(),
                                                                src(S0),
                                                                V(room_number),// Init. volume vector
                                                                time_elapsed(steps_number+1),
                                                                T(room_number, room_number),   // Init. transport matrix 
                                                                S(room_number),// Init. source vector
                                                                D(room_number),// Init. deposition factors vector
                                                                N(steps_number+1, room_number)// Init. of matrix with activities in each room (on rows) at different times (on columns) 
                                                                {
                                                                    this->Init();
                                                                }


Transport::~Transport(){}
//############################# INIT ##########################################

void Transport::Init(){    // Init. volume vector, reference [3]

    V(0) = 42;
    V(1) = 812;
    V(2) = 40700;
    V(3) = 1; //EXTERNAL VOLUME NOT RELEVANT

    //constant flows, reference [3]
    double f01= 0.03;
    double f12= 0.045167;
    double f23= 2.261167;
    T(0,1) = f01*pTime_step/V(0);
    T(1,2) = f12*pTime_step/V(1);
    T(2,3) = f23*pTime_step/V(2);
    //define the diagonal Tii as 1 - Di - sum_j(Tij)
    TMatrixDDiag diagT(T); diagT = 0;
    TVectorD sumOutflow(room_number); 
    for (int i=0; i<room_number; i++){
        for(int j=0; j<room_number; j++){
            if(j!=i) sumOutflow(i) += T(i,j);
        }
        
        diagT(i)= 1 - D(i) - sumOutflow(i);
    }

}
void Transport::Go(){
//############################# TRANSPORT #####################################
    double time_step = pTime_step;
    steps_number = (int) pFinal_time/time_step;
    int time_mark = 0;
    emission_fraction = {0.01, 0.1, 0.5, 0.9, 0.99};

    for (int step=0; step<steps_number+1; step++){

        time_elapsed(step)=step*time_step;

        S(0) = src.Eval(step*time_step)*time_step;

        if(step==0){
            TVectorD row = S;
            TMatrixDRow(N,step) = row;
        } else {
            TVectorD row = TMatrixDRow(N, step);
            TVectorD row_before = TMatrixDRow(N, step-1);
            
            row = S + Mult(row_before,T);
            row *= TMath::Exp(-src.getDecayConst()*time_step); 
    
            TMatrixDRow(N,step) = row;
            
            if(time_mark<emission_fraction.size() && row(room_number-1)>emission_fraction.at(time_mark)){
                emission_times.push_back( (step-1)*time_step );
                time_mark++;
            }
    
        }

        

    }
    
}
    void Transport::setSource( Source s){src=s;}; 
    void Transport::setTime_step( double time_step ){pTime_step=time_step;this->Init();};
    void Transport::setFinal_time( double final_time){pFinal_time=final_time;this->Init();}
    void Transport::setMatrixN( TMatrixD newN){N = newN;}

    Source Transport::getSource() const {return src;}
    TMatrixD Transport::getMatrixN() const {return N;}
    TVectorD Transport::getTimeElapsed() const {return time_elapsed;}
    TMatrixD Transport::getMatrixT() const {return T;}
    double Transport::getTime_step() const {return pTime_step;}
    int Transport::getRoom_number() const {return room_number;}
    vector<double> Transport::getEmission_fraction() {return emission_fraction;}
    vector<double> Transport::getEmission_times(){return emission_times;} 