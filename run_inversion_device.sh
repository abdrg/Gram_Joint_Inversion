#!/bin/bash
# Script con comandos necesarios para compilar y correr el código forward_potential.

# Erase previous output files.
rm output_*.csv
rm output_*.xyz

rm -f output_misfit_jointGram_c1.txt
rm -f output_parameters_inv.txt
rm -f output_time.txt
rm -f run_inversion.exe
rm -f *.o
rm -f *.mod


# Al no tener funcion main la bandera -c indica que se salte esa advertencia.
#gfortran -pg -c auxiliar.f08 -o auxiliar.o
#gfortran -pg -c operators.f08 -o operators.o
#gfortran -pg -c model_ini.f08 -o model_ini.o
#gfortran -pg -c forward_modeling.f08 -o forward_modeling.o
#gfortran -pg -c inversion_serial.f08 -o inversion_serial.o -llapack -lblas
#gfortran -pg -c forward_script.f08 -o forward_script.o
#gfortran -pg -c inversion_script.f08 -o inversion_script.o

# Compilation.
nvfortran -cuda -gpu=rdc -c auxiliar.f08 -o auxiliar.o
nvfortran -cuda -gpu=rdc -c operators.f08 -o operators.o
nvfortran -cuda -gpu=rdc -c model_ini.f08 -o model_ini.o
nvfortran -cuda -gpu=rdc -c forward_modeling.f08 -o forward_modeling.o
#nvfortran -cuda -gpu=rdc -c inversion_serial.f08 -o inversion_serial.o -llapack -lblas
nvfortran -cuda -gpu=rdc -c inversion_parallel.cuf -o inversion_parallel.o -cudalib=cublas -lblas -llapack
nvfortran -cuda -gpu=rdc -c inversion_script.f08 -o inversion_script.o

# Compilar el programa completo usando todos los modulos y objetos compilados anteriormente.
nvfortran *.o -cuda -gpu=rdc -cudalib=cublas -lblas -llapack #-Mlarge_arrays #-cuda=fastmath # -g -O3 -Mfprelaxed  # -Mcuda=ptxinfo

# Execution.
mv a.out run_inversion.exe
./run_inversion.exe

#Profiling
#cuda-memcheck ./run_inversion.exe
#cuda-gdb ./run_inversion.exe
#gprof ./run_inversion.exe

# Clean up temporary files.
rm -f run_inversion.exe
rm -f *.o
rm -f *.mod
#rm -f fort.6
#rm -f gmon.out


# Renombrar archivos de salida modelado sintetico para volver archivos de entrada de script inversion
#cp output_gv_anomaly_m_ini_noised.csv input_d_gv.csv
#cp output_mg_anomaly_m_ini_noised.csv input_d_mg.csv
#cp output_gv_m_ini.csv input_m_gv.csv
#cp output_mg_m_ini.csv input_m_mg.csv

# Cambiar extension modelos salida para poder graficar directo en 3D/voxel_convert de oasis montaj
#cp output_jointGram_gv_minv_c1_iter9.csv output_jointGram_gv_minv_c1_iter9.xyz
#cp output_jointGram_mg_minv_c1_iter9.csv output_jointGram_mg_minv_c1_iter9.xyz

# Para cuando el tamaño de I es de dos digitos (I2 en vez de I1)
cp output_jointGram_gv_minv_c1_iter\ 1.csv output_jointGram_gv_minv_c1_iter\ 1.xyz
cp output_jointGram_mg_minv_c1_iter\ 1.csv output_jointGram_mg_minv_c1_iter\ 1.xyz
cp output_jointGram_gv_minv_c1_iter\ 5.csv output_jointGram_gv_minv_c1_iter\ 5.xyz
cp output_jointGram_mg_minv_c1_iter\ 5.csv output_jointGram_mg_minv_c1_iter\ 5.xyz
cp output_jointGram_gv_minv_c1_iter10.csv output_jointGram_gv_minv_c1_iter10.xyz
cp output_jointGram_mg_minv_c1_iter10.csv output_jointGram_mg_minv_c1_iter10.xyz
cp output_jointGram_gv_minv_c1_iter15.csv output_jointGram_gv_minv_c1_iter15.xyz
cp output_jointGram_mg_minv_c1_iter15.csv output_jointGram_mg_minv_c1_iter15.xyz
cp output_jointGram_gv_minv_c1_iter20.csv output_jointGram_gv_minv_c1_iter20.xyz
cp output_jointGram_mg_minv_c1_iter20.csv output_jointGram_mg_minv_c1_iter20.xyz
cp output_jointGram_gv_minv_c1_iter25.csv output_jointGram_gv_minv_c1_iter25.xyz
cp output_jointGram_mg_minv_c1_iter25.csv output_jointGram_mg_minv_c1_iter25.xyz
#cp output_jointGram_gv_minv_c1_iter30.csv output_jointGram_gv_minv_c1_iter30.xyz
#cp output_jointGram_mg_minv_c1_iter30.csv output_jointGram_mg_minv_c1_iter30.xyz
#cp output_jointGram_gv_minv_c1_iter35.csv output_jointGram_gv_minv_c1_iter35.xyz
#cp output_jointGram_mg_minv_c1_iter35.csv output_jointGram_mg_minv_c1_iter35.xyz

###################################################################################################################################
# -c flag para no hacer vinculo, permite compilar modulos por separado, generando modulos .mod y objetos .o
# -o flag indica que nombre al archivo de salida.
# -Wall flag genera reporte de errores.
# -llapack flag incluye la libreria lapack.
# -lblas flag incluye la libreria blas.
# -cuda flag indica que se usara cuda.
# -gpu=rdc flag indica Relocatable Device Code para poder hacer uso de modulos de CudaFor en otro archivo.
# -g -G flag indica que se generara informacion de depuracion de host y de device para cuda. 
# -traceback flag para imprimir la pila de llamadas cuando se genere un error.
# -Ofast  This option enables aggressive optimizations at the cost of reduced floating point precision.
# -O3 -Mfprelaxed  This option enables aggressive optimizations at the cost of reduced floating point precision.
