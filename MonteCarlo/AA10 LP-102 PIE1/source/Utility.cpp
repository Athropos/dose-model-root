#include "../headers/Utility.h"


#include <fstream>


void print_row(TVectorD row, std::ofstream& file){
    int size = row.GetNoElements();

    for(int i=0;i<size; i++){
        file << row(i) << ' ' ;
    }
    file << endl;
    }

TVectorD Mult(TVectorD  vec, TMatrixD  matrix){
    int size = vec.GetNoElements();
    TVectorD out(size);

    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            out(i) += vec(j) * matrix(j,i);
        }
    }
    
    return out;
}

//Triangular Probability Density Function
Double_t WaterPDF(Double_t *x, Double_t *par){
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



// TUnuranEmpDist windEmpDist(){

//     ifstream in;
//     in.open("wind_speed");
//     double wind_speed;
//     vector<double> data;

//     while(in >> wind_speed){
//   	    data.push_back(wind_speed);
//     }
//     in.close();
//     TUnuranEmpDist dist( data.begin(), data.end());

//     return dist;
// }

// TUnuranDiscrDist stabilityDiscrDist(){

//     double stability[6] = {0, 0.04, 0.09, 0.67, 0.09, 0.11};
//     TUnuranDiscrDist dist( stability, stability+6);
    
//     return dist;
// }

const char* stability_conv( int num){

    switch (num){
        case 0:
            return "A";
        break;

        case 1:
            return "B";
        break;

        case 2:
            return "C";
        break;

        case 3:
            return "D";
        break;

        case 4:
            return "E";
        break;

        case 5:
            return "F";
        break;

        default:
            cout << "Non valid number! Insert only numbers from 0 to 5\n"; 
            return "0";
        break;
    } 
}  

MyTF1 Sigma_zPDF(int stability_class){
    MyTF1 A("A","ROOT::Math::lognormal_pdf(x,3.97,0.13)",30,90);
    MyTF1 B("B","ROOT::Math::lognormal_pdf(x,3.54,0.09)",10,90);
    MyTF1 C("C","ROOT::Math::lognormal_pdf(x,3.35,0.05)",10,90);
    MyTF1 D("D","ROOT::Math::lognormal_pdf(x,3.21,0.04)",10,90);
    MyTF1 E("E","ROOT::Math::lognormal_pdf(x,3.09,0.03)",10,90);
    MyTF1 F("F","ROOT::Math::lognormal_pdf(x,3.03,0.02)",10,90);

    switch (stability_class)
    {
    case 0:
        return A;
        break;
    case 1:
        return B;
        break;
    case 2:
        return C;
        break;
    case 3:
        return D;
        break;
    case 4:
        return E;
        break;
    case 5:
        return F;
        break;
                        
    default:
        cout << "Invalid stability class\n";
        break;
    }
    

}

Double_t Deposition_vel_C14_PDF(Double_t *x, Double_t *par){
    //Parameters
    Double_t xx= x[0];

    if(xx>=0.001 && xx<0.0015) {return (2./3.)/0.0005;}        //winter prob.
    else if(xx>=0.0015 && xx<=0.013) {return (1./3.)/0.0115;}  //summer prob.
    else{return 0.;}
}

