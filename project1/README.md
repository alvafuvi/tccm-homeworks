In this program, users can calculate the Hartree-Fock and MP2 energy for different molecules.

It is written completely in C and uses the external library TREXIO (https://github.com/TREX-CoE/trexio). This library depends on another one called HDF5, so in order for this program to work, it is necessary to install both.

Input file examples have been added in the _test directory. They can be generated using TREXIO writing functions as: ```trexio_write_{group}_{attribute}```.

Regarding execution, after compiling the source code (more info about this in ```INSTALL.md``` ) the program must be executed as:

```
./mp2 [file.h5]
```

The argument here has to be the path to the TREXIO input file (if the input file is in the same directory only the name is needed). By doing this, some data will be printed in the terminal, including the Hartree-Fock and MP2 energy in a.u.


