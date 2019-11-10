#include "../headers/Utility.h"

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

