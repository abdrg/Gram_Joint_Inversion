#!/bin/bash
# Script with commands needed to compile and run code semi-automatically.

# Erase previous output files.
rm output_*.csv
rm output_*.xyz

# Compilation.
#gfortran -pg -c gram_joint_inversion_CPU.f08 -o gram_joint_inversion_CPU.o -llapack -lblas
gfortran -pg -c auxiliar.f08 -o auxiliar.o
gfortran -pg -c parameters.f08 -o parameters.o
gfortran -Wall -g -pg gram_joint_inversion_CPU.o auxiliar.o parameters.o -llapack -lblas -Ofast

# Execution.
mv a.out run_inversion.exe
./run_inversion.exe

#Profiling.
#gprof ./run_inversion.exe

# Clean up temporary files.
rm -f ./run_inversion.exe
rm -f gmon.out


# Change some model results to .xyz format for plotting.
#cp output_jointGram_gv_minv_iter\ 7.csv output_jointGram_gv_minv_iter\ 7.xyz 
#cp output_jointGram_mg_minv_iter\ 7.csv output_jointGram_mg_minv_iter\ 7.xyz
#cp output_jointGram_gv_minv_iter30.csv output_jointGram_gv_minv_iter30.xyz 
#cp output_jointGram_mg_minv_iter30.csv output_jointGram_mg_minv_iter30.xyz