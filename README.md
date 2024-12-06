# YHPIC
YHPIC is a highly scalable and efficient parallel solver for PIC plasma simulations on the MT-3000 processor

### YHPIC-CPU
YHPIC-CPU is a program that can only run on the CPU, and the work in the paper involves porting several core functions of this program to the MT-3000 processor.

1. MPI Parallel Computing Environment
The YHPIC-CPU program has been successfully tested in the MPICH 3.0 or later environment. To install the MPI environment on Ubuntu Linux, you can run the following command:

  $ sudo apt-get install mpich libmpichdev

2. HDF5 Scientific Data Repository
The YHPIC-CPU program uses HDF5 for large-scale scientific data storage. To install the HDF5 environment on Ubuntu Linux, run the following command:

  $ sudo apt-get install libhdf5-parallel-dev

3. Compilation and Installation
Once the MPI and HDF5 libraries and environment are installed correctly, unzip the YHPIC-CPU.tar.gz file. A YHPIC-CPU directory will be created in the current directory, which includes a Makefile src/Makefile. By simply configuring the installation directories for MPI and the HDF5 libraries in the Makefile, the program can be compiled directly by running the command:

  $ make

By default, the program is compiled with the -O3 optimization option, and the executable program  will be generated in the /run/ directory. Finally, for ease of visualizing and analyzing computational data, users should also install the Paraview and Matplotlib programs and libraries.

4. Running
Navigate to the run subdirectory and run the following command:
  $ mpirun -np N ./lared-p   # N represents the number of parallel CPU cores.
