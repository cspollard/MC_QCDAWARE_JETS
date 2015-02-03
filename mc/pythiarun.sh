#!/usr/bin/bash

CMNDFILE=$1
HEPMCFILE=${1/cmnd/hepmc}
LOGFILE=${1/cmnd/log}
YODAFILE=${1/cmnd/yoda}

mkfifo $HEPMC

run-pythia -s -i $CMNDFILE -o $HEPMCFILE 2>&1 > $LOGFILE

rivet --pwd -a MC_QCDAWARE_JETS $HEPMCFILE -H $YODAFILE 2>&1 > $LOGFILE

rm $HEPMCFILE
