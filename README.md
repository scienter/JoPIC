# JoPIC

1. Installation

I recommand installing in LINUX system.
Present version of codes is composed of 'cartesian' and 'cylindrical' codes.

In each directory, you may find 'makefile'.

Denpending on your computer system, specify 'CFLAGS' and 'LDFLAGS' which may include pathes of gs library(gsl), fftw, and hdf5.

If you don't have gsl, fftw, and hdf5 in you computer, please install three libraries.

After defining pathes, type 'make'.



2. Running codes

Type below,

mpirun -np [number of mpi_cores] [path of code]/jopic [input file]




3. Input file

... in preparation.
