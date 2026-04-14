#!/bin/bash

for j in $(seq 0 29); do
    Rscript simulation_1_tps.R 250 $j 1.00
done  


