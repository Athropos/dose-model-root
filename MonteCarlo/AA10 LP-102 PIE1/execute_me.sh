#! /bin/bash

g++ -o accident macroSobol.cpp `root-config --cflags --libs` source/Source.cpp source/Sobol.cpp source/Utility.cpp source/Transport.cpp source/OuterSource.cpp source/Dispersion.cpp source/Dose.cpp 
./accident

