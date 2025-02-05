#!/bin/bash

# Path
WORKDIR=.
SCRIPT=${WORKDIR}/detect.C
INPUT=${WORKDIR}/../tutorial.root
OUTPUT=${WORKDIR}/../tutorial_detect.root

# Execute
time root -l -b -q ${SCRIPT}\(\"${INPUT}\",\"${OUTPUT}\"\)
