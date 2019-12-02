#! /bin/bash

`root-config --cxx --cflags` -o accident macroBCI.cpp  source/Source.cpp source/Utility.cpp source/Sobol.cpp source/Transport.cpp source/OuterSource.cpp source/Dispersion.cpp source/Dose.cpp source/Model.cpp `root-config --libs`
./accident


