#! /bin/bash

`root-config --cxx --cflags` -o accident macro.cpp source/Source.cpp source/Utility.cpp source/Transport.cpp source/Output.cpp source/OuterSource.cpp source/Dispersion.cpp source/Dose.cpp `root-config --libs`
./accident
