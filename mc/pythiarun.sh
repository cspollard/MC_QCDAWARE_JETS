#!/bin/bash

CMNDFILE=$1
HEPMCFILE=`mktemp`.hepmc
PYTHIALOG=${1/cmnd/pythia.log}
RIVETLOG=${1/cmnd/rivet.log}
YODAFILE=${1/cmnd/yoda}

NEVT=$2

mkfifo $HEPMCFILE

run-pythia -s -n $NEVT -e 13000 -i $CMNDFILE \
    -o $HEPMCFILE > $PYTHIALOG 2>&1 &

rivet --pwd -a MC_QCDAWARE_JETS $HEPMCFILE \
    -H $YODAFILE > $RIVETLOG 2>&1

rm $HEPMCFILE
