#!/bin/bash
# Script with commands needed to compile and run code semi-automatically.

# Erase previous output files.
rm output_*.csv
rm output_*.xyz

# Compilation.
gfortran -pg -c gram_joint_inversion_host.f08 -o gram_joint_inversion_host.o -llapack -lblas
gfortran -pg -c parameters.f08 -o parameters.o
gfortran -Wall -g -pg *.o -llapack -lblas -Ofast

# Execution.
mv a.out run_inversion.exe
./run_inversion.exe

# Change some model results to .xyz format for plotting.
cp output_jointGram_gv_minv_iter\ 7.csv output_jointGram_gv_minv_iter\ 7.xyz 
cp output_jointGram_mg_minv_iter\ 7.csv output_jointGram_mg_minv_iter\ 7.xyz
cp output_jointGram_gv_minv_iter30.csv output_jointGram_gv_minv_iter30.xyz 
cp output_jointGram_mg_minv_iter30.csv output_jointGram_mg_minv_iter30.xyz

#Profiling.
#gprof ./run_inversion.exe

# Clean up temporary files.
rm -f ./run_inversion.exe
rm -f *.o
rm -f *.mod
#rm -f fort.6