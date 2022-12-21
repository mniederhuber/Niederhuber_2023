#!/bin/bash


peakfile=$1
fn=$(basename -- "${peakfile}")
fn=${fn%.*}

echo $fn
awk '$9 >= 10 {print}' $peakfile > "Peaks/${fn}-q10.narrowPeak"
