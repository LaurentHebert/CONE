#! /bin/bash

python3 decomposition $1

g++ -std=c++11 -O3 -o tevol_source_sis_cone ./tevol_source_sis_cone.cpp $(gsl-config --cflags) $(gsl-config --libs)

rm ./results/results.dat

for beta in {00..99}; do ./tevol_source_sis_cone '0.'$beta >>'./results/results.dat'; done

python3 plot_figure_localization.py