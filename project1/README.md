In this program, users can calculate the Hartree-Fock and MP2 energy for different molecules.

It is written completely in C and uses the external library TREXIO (https://github.com/TREX-CoE/trexio). This library depends on another one called HDF5, so in order for this program to work, it is necessary to install both.

Input file examples hve been added in the test directory. They can be generated using TREXIO writing functions as: trexio_write_{group}_{attribute}.

After executing, the program will ask in the prompt to choose a TREXIO input file (it has to be in the same directory as the executable). By doing this, some data will be printed in the terminal, including the Hartree-Fock and MP2 energy in a.u.


