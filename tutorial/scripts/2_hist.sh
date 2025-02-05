#!/bin/bash

# Path
WORKDIR=.
SCRIPT=${WORKDIR}/hist.C
INPUT=${WORKDIR}/tutorial_detect.root
OUTPUT=${WORKDIR}/tutorial_hist.root

# Execute
time root -l -b -q ${SCRIPT}\(\"${INPUT}\",\"${OUTPUT}\"\)
