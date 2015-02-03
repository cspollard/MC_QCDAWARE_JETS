#!/usr/bin/bash

mkfifo $1.hepmc

run-pythia -s -i $1.cmnd -o $1.hepmc 2>&1 > $1.log

rivet --pwd -a MC_QCDAWARE_JETS $1.hepmc -H $1.yoda

rm $1.hepmc
