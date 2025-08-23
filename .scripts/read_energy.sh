#!/bin/bash

DIR=$HOME/C5O-Kinetics/calc/$1
cd $DIR

grep -R --include "REVDSD.log" "Zero point energy" 
grep -R --include "CCSDT.log" "FINAL SINGLE POINT ENERGY" 