#!/bin/bash

CMNDFILE=$1
HEPMCFILE=${1/cmnd/hepmc}
PYTHIALOG=${1/cmnd/pythia.log}

NEVT=$2

run-pythia -s -n $NEVT -e 13000 -i $CMNDFILE \
    -o $HEPMCFILE > $PYTHIALOG 2>&1 &
