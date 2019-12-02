#! /bin/bash

g++ -o accident macro.cpp `root-config --cflags --libs` source/Source.cpp source/Utility.cpp source/Transport.cpp source/Output.cpp source/OuterSource.cpp source/Dispersion.cpp source/Dose.cpp
./accident

