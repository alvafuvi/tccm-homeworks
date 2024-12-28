    Default compilation
=============================

For compilation, one should use the following command considering that the installation of TREXIO was successful: 

gcc -I/usr/local/include -L/usr/local/lib mp2.c -o <__name_of_executable__> -ltrexio

If compiling was successful, there should be an executable binary file in the working directory.


    Compiling with CMake
==============================

For CMake compiling, locate yourself in the build directory and jus execute the following commands

    ~$ cmake ..
    ~$ make

The executable file should be created in the source directory.
