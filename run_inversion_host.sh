#!/bin/bash
# Script con comandos necesarios para compilar y correr el código forward_potential.

# Erase previous output files.
rm output_*.csv
rm output_*.xyz

#rm output_d_mg.csv
#cp output_d_gv_noised.csv input_d_gv.csv
#cp output_d_mg_noised.csv input_d_mg.csv
#cp output_m_gv_ini.csv input_m_gv.csv
#cp output_m_mg_ini.csv input_m_mg.csv
#cp output_a_gv.dat input_a_gv.dat
#cp output_a_mg.dat input_a_mg.dat
#cp output_dx.dat input_dx.dat
#cp output_dy.dat input_dy.dat
#cp output_dz.dat input_dz.dat
#cp output_Lapl.dat input_Lapl.dat

rm -f output_misfit_jointGram_c1.txt
rm -f output_time.txt
rm -f ./run_inversion.exe
rm -f *.o
rm -f *.mod


# Al no tener funcion main la bandera -c indica que se salte esa advertencia.
gfortran -pg -c auxiliar.f08 -o auxiliar.o
gfortran -pg -c operators.f08 -o operators.o
gfortran -pg -c model_ini.f08 -o model_ini.o
gfortran -pg -c forward_modeling.f08 -o forward_modeling.o
gfortran -pg -c inversion_serial.f08 -o inversion_serial.o -llapack -lblas
gfortran -pg -c forward_script.f08 -o forward_script.o
#gfortran -pg -c inversion_script.f08 -o inversion_script.o

# Compilar el programa completo usando todos los modulos y objetos compilados anteriormente.
gfortran -Wall -g -pg *.o -llapack -lblas -Ofast

# Renombra el ejecutable .out y corre el nuevo ejecutable .exe
mv a.out run_inversion.exe
./run_inversion.exe

#Profiling
#gprof ./run_inversion.exe

rm -f ./run_inversion.exe
rm -f *.o
rm -f *.mod
rm -f fort.6
rm -f gmon.out

# Renombrar archivos de salida modelado sintetico para volver archivos de entrada de script inversion
#cp output_gv_anomaly_m_ini_noised.csv input_d_gv.csv
#cp output_mg_anomaly_m_ini_noised.csv input_d_mg.csv
#cp output_gv_m_ini.csv input_m_gv_test.csv
#cp output_mg_m_ini.csv input_m_mg_test.csv
#rm output_gv_anomaly_m_ini.csv
#rm output_mg_anomaly_m_ini.csv
#rm output_gv_anomaly_m_ini_noised.csv
#rm output_mg_anomaly_m_ini_noised.csv
#cp output_gv_m_ini.csv input_m_gv.csv
#cp output_mg_m_ini.csv input_m_mg.csv
#rm output_gv_m_ini.csv
#rm output_mg_m_ini.csv

# Cambiar extension modelos salida para poder graficar directo en 3D/voxel_convert de oasis montaj
#cp output_jointGram_gv_minv_c1_iter4.csv output_jointGram_gv_minv_c1_iter4.xyz
#cp output_jointGram_mg_minv_c1_iter4.csv output_jointGram_mg_minv_c1_iter4.xyz

# Para cuando el tamaño de I es de dos digitos (I2 en vez de I1)
#cp output_jointGram_gv_minv_c\ 1_iter\ 5.csv output_jointGram_gv_minv_c\ 1_iter\ 5.xyz
#cp output_jointGram_mg_minv_c\ 1_iter\ 5.csv output_jointGram_mg_minv_c\ 1_iter\ 5.xyz

###################################################################################################################################
# -c flag para no hacer vinculo, permite compilar modulos por separado, generando modulos .mod y objetos .o
# -o flag indica que nombre al archivo de salida.
# -g flag genera informacion de depuracion.
# -Wall flag genera reporte de errores.
# -pg flag para hacer profiling con gprof.
# -llapack flag incluye la libreria lapack.
# -lblas flag incluye la libreria blas.
# -cuda flag indica que se usara cuda.
# -gpu=rdc flag indica Relocatable Device Code para poder hacer uso de modulos de CudaFor en otro archivo.
# -Ofast muchas optimizaciones del compilador para hacerlo mas rapido. Puede generar errores de presicion.
