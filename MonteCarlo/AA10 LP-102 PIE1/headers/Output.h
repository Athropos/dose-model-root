#ifndef OUTPUT_H
#define OUTPUT_H

#include "Transport.h"
#include "TFile.h"

#include <fstream>

class Output
{
    public: 
    Output(); 
    virtual ~Output();

    void Root(Transport,TFile*);
    void OuterInventory(Transport,double, ofstream& );
    private:
    TFile* file_root;

};
#endif


