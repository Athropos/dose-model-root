#! /bin/bash

`root-config --cxx --cflags` -o example example.cpp  source/Source.cpp `root-config --libs`
./example

