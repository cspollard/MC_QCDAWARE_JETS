#!/bin/bash

CMNDFILE=$1
HEPMCFILE=`mktemp`.hepmc
LOGFILE=${1/cmnd/log}
YODAFILE=${1/cmnd/yoda}

NEVT=$2

mkfifo $HEPMCFILE

run-pythia -s -n $NEVT -e 13000 -i $CMNDFILE \
    -o $HEPMCFILE > $LOGFILE 2>&1 &

rivet --pwd -a MC_QCDAWARE_JETS $HEPMCFILE \
    -H $YODAFILE > $LOGFILE 2>&1

rm $HEPMCFILE
