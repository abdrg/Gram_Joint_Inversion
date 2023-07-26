#!/bin/bash
# Script con comandos necesarios para compilar y correr el código forward_potential.

# Borrar resultados experimentos previos
rm output_*.csv
rm output_*.xyz

# Compilacion
nvfortran -cuda -gpu=rdc -c gram_joint_inversion.cuf -o gram_joint_inversion.o -cudalib=cublas -lblas
nvfortran -cuda -gpu=rdc -c parameters.f08 -o parameters.o
nvfortran *.o -cuda -gpu=rdc -cudalib=cublas -lblas -llapack -Mlarge_arrays #-cuda=fastmath #-g -O3 -Mfprelaxed  # -Mcuda=ptxinfo

# Ejecucion
mv a.out run_inversion.exe
./run_inversion.exe

#gedit output_misfit_jointGram.txt
cp output_jointGram_gv_minv_iter\ 7.csv output_jointGram_gv_minv_iter\ 7.xyz 
cp output_jointGram_mg_minv_iter\ 7.csv output_jointGram_mg_minv_iter\ 7.xyz
cp output_jointGram_gv_minv_iter30.csv output_jointGram_gv_minv_iter30.xyz 
cp output_jointGram_mg_minv_iter30.csv output_jointGram_mg_minv_iter30.xyz
cp output_jointGram_gv_minv_iter15.csv output_jointGram_gv_minv_iter15.xyz 
cp output_jointGram_mg_minv_iter15.csv output_jointGram_mg_minv_iter15.xyz

#Profiling
#cuda-memcheck ./run_inversion.exe
#cuda-gdb ./run_inversion.exe
#gprof ./run_inversion.exe

rm -f ./run_inversion.exe
rm -f *.o
rm -f *.mod
#rm -f fort.6
