# JoPIC

# 1. Installation

I recommand installing 'JoPIC' code in LINUX system.

Present version is composed of 'cartesian' and 'cylindrical' codes and each code is located in different directory.

In directories, you may find 'makefile' which will install the code.
Denpending on your computer system, you have to specify '**CFLAGS**' and '**LDFLAGS**' which will define pathes of gs library(gsl), fftw, and hdf5.

If you don't have gsl, fftw, and hdf5 in you computer, please install three libraries.

After defining pathes in '**makefile**', type '**make**'.


# 2. Running codes

Type below,

**mpirun -np [number of mpi_cores] [path of code] / jopic [input file] [step for restart]**

**[step for restart]** : iteration step only when dump file (titled with '**dumpField[step].h5**' and '**dumpParticle[step].h5**' exist.

*example :*

* When you want to rerun starting 5000 time-step with existing '**dumpField5000.h5**' and '**dumpParticle5000.h5**'.

  **mpirun -np 8 /home/scienter/Git/JoPIC/cartesian/jopic lwfa.inp 5000**

* When you run from beginning

  **mpirun -np 8 /home/scienter/Git/JoPIC/cartesian/jopic lwfa.inp**



# 3. Input file

Some examples are prepared in '**input**' folder. Input file can be any txt file. Examples have '.inp' extension, but there is no rule to define file name.
