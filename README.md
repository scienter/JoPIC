# JoPIC
Send E-mail if you need help in using JoPIC.

**mh0309@postech.ac.kr**

# 1. Installation

I recommand installing 'JoPIC' code in LINUX system.

Present version is composed of 'cartesian' and 'cylindrical' codes and each code is located in different directory.

In directories, you may find 'makefile' which will install the code.
Denpending on your computer system, you have to specify '**CFLAGS**' and '**LDFLAGS**' which will define pathes of gs library(gsl), fftw, and hdf5.

If you don't have gsl, fftw, and hdf5 in you computer, please install three libraries.

After defining pathes in '**makefile**', 

type '**make**'.


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

Some examples are prepared in '**input**' folder. Input file can be any txt file. Examples have '.inp' extension, but there is no rule to define file name and extension.

You may find 
In writing....

# 4. Postprocessing

Output files are written in hdf5 file format. Output files can be converted or explored using python or c programming.

In '**diag**' folder, you may find some postprocessing codes, those convert hdf5 output to txt file. 

The converted txt file will write and arrange data to be able to plot in 'gnuplot' program which is well-known free software and available for windows and linux system.

To use diag programs, you have to compile codes. In **diag** directory, you may find **compile** file. After specifying library pathes in **compile**, simply type

**sh compile** or **./compile**

Once finishing compilation, you will have list of programs.

# *Postprocessing field data.*

When you type './hdf_field', it will show message as below.

**./hdf_field [mode] [initial] [final] [timestep] [stX] [stY]

mode(0): [fileType] [angle] [minY] [maxY] [shift]

mode(1): [fileType] [angleDivision] [cellNum]

mode(2): [minY] [maxY] [angle] [numS] [species1] [species2] ...

mode(3 : summing): [minY] [maxY] [angle] [numS] [species1] [species2] ...**

There are 4 modes, but mode 0 will be enough for obtaining basic data.

**[stX]** **[stY]** : defining stripe step. 
If **[stX]**=1, **[stY]**=1, all field grids are included. 
If **[stX]**=1, **[stY]**=2, y grid number will be half of origin data. 

**[fieldType]** : file name except of time step. 
If field file has name of **fieldSplit4000.h5**, **[fieldType]** is **fieldSplit**.
In case of density file like **0density5000.h5**, **[fieldType]** is **0density**.

**[angle]** : azimuthal angle in degree. ex) 0 or 90

**[minY] [maxY]** : plot range in transverse

**[shift]** : just type '0'.

In writing....

In writing....


Other postprocessing program written in python will be preparing in the future.



