# Gram_Joint_Inversion
Fortran Parallel Code to perform Joint Inversion of Gravity and Magnetic with structural similarity approach based on Gramian constraints.

This repository contains all you need to replicate a synthetic experiment proposed by Li & Oldenburg (1998) and try to obtain better results by using a modern robust joint inversion technique.

Main program 'parameters.f08' allows you to setup all inversion parameters, meanwhile 'gram_joint_inversion.cuf' and 'gram_joint_inversion_host.f08' contains different code implementations to calculate solutions and some auxiliar subroutines. Scripts 'run_inversion_XXX.sh' allow you to compile and execute the experiment easily.

To run the code only with CPU you should use the script 'run_inversion_host.sh'. This option is the slowest but it can be run in any computer with gfortran compiler preinstalled. 
To run the code on GPU you must have a PC with dedicated graphics card that supports CUDA, preinstall the NVIDIA HPC SDK (https://developer.nvidia.com/hpc-sdk-downloads) and use the script 'run_inversion_device.sh'. This option allows to test CPU, GPU and a CPU-GPU hybrid versions of the code. 

You can download binaries of matrices A_gv, A_mg, dx, dy, dy and Lapl for this specific experiment in the following link (https://drive.google.com/drive/folders/16lP8HbDLCub-WlDDTDghrTXKrTE3ZNza?usp=sharing), which has been precalculated to isolate their computational cost from the main program in order to evaluate the improvement on performance with different code implementations.


References
Li, Y., & Oldenburg, D. W. (1998). 3-D inversion of gravity data. GEOPHYSICS, 63(1), 109–119. http://segdl.org/
