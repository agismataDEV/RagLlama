by Agisilaos Matalliotakis
24 January 2020, Medical Imaging Group, TU Delft , Delft, The Netherlands

INCS codes uses multiple libraries in order to run in a right and an efficient way.
The libraries that are used are static libraries and not shared. The difference in how to separate 
them is mainly that static start with "lib" like lidhdf5_fortan while the respective shared library would be hdf5_fortran
These libraries, bin and include files should be included when the code runs.
It is possible that I have included more files but with just these the code works in a fine way.
The libraries that are used are:
1. FFTW3 library because of the version 3.3.8 . Any fftw library should be ok, even the future versions. You can add them with -lfftw3 and -lfftw3_threads

2. HDF5 library v1.8.7 or previous can be used because parallel is enabled for CentOS v6.10 (Current Linux distribution)
   HDF5 v1.8.21 was tested but could not be compiled with --enable-parallel configuration. Release Notes for every version
   write exactly in which platforms the library was used and for what mode. You can add them with -lhdf5_fortran and -lhdf5.
   Others libraries that hdf5 has a dependence on are zlib and szlib so if they are built they can be used. Otherwise you should build
   them before trying to build hdf5. In case hdf5 goes to old or it is impossible to build them with parallel mode then another output library
   should be used.

   In order for HDF5  library to work, it should be compiled by the proper compiler. IF the compiler used, is the Intel Compiler(x64) 
   then something like that should be used to build the libraries manually in a linux distribution : 

	module load openmpi
	export CC=mpiicc
	export FC=mpiifort
	export F9X=mpiifort
	export F90=mpiifort
	export CXX=mpiicpc
	export MPICC=mpiicc
	export FFLAGS='-fpp -DDEC$=DEC_ -DMS$=MS_'
	source /opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin/mpivars.sh
	./configure --prefix=/home/agismatalliota/hdf5_ifort --with-zlib --with-szlib --enable-fortran --enable-parallel --enable-production

3. MPI library for parallel computing . In this case , because the compiler is Intel Compiler, for Fortran90, mpiifort is used to do the compilation 
   and mpiexec to run the file with multiple cores.

Makefile contains all the needed commands for compilation but some explanation follow to clarify the process:

i. The objects that will be used for compilation should be written in a serial way and their position in the file is really important. This is because
   in order to compile a file, all the dependencies should be compiled first . So the dependency order should be satisfied.

ii. The commands for the compilation that were used are :

	for the libraries LIB and the include INC path were used from Liberario Demi's work and hdf5 libs and inc for the reasons mentioned before.
	LIB_PATH_INCS = -L ~/INCS_Agis/LIB -L ~/hdf5_ifort_v1.8.7/lib #-L ~/hdf5_ifort_v1.8.21/lib
	INC_PATH_INCS = -I ~/INCS_Agis/INCLUDE -I ~/hdf5_ifort_v1.8.7/include #-I ~/hdf5_ifort_v1.8.21/include
	
	The libraries used are static as mentioned for mpi parallel programming, hdf5, zlib , fftw3  and a standalone C mathematical library libm -lm
	LIBS = -lmpifort -lhdf5_fortran -lhdf5 -lz -lhdf5_hl -lhdf5hl_fortran -lfftw3 -lfftw3_threads -lm  -static-intel

	This flag was used for a warning in some files Bad # preprocessor line #endif -^
	CFLAGS = -fpp 



!!!In case errors appear for MPI , load openmpi and export PATH and Libraries like :

%% openmpi should not be loaded
export PATH=$PATH:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/bin/intel64:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/bin:/opt/intel/debugger_2017/gdb/intel64_mic/bin:/opt/openlava-3.0/bin:/bin:/usr/bin:/usr/sbin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/openlava-3.0/bin:/home/agismatalliota/bin:/opt/openlava-3.0/bin:/opt/openlava-3.0/bin:/home/agismatalliota/bin:/opt/openlava-3.0/bin:${MKLROOT}/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/mpi/intel64/lib:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/mpi/mic/lib:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/ipp/lib/intel64:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/compiler/lib/intel64_lin:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64_lin:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/tbb/lib/intel64/gcc4.4:/opt/intel/debugger_2017/iga/lib:/opt/intel/debugger_2017/libipt/intel64/lib:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/daal/lib/intel64_lin:/opt/insy/intel/2017u4/compilers_and_libraries_2017.4.196/linux/daal/../tbb/lib/intel64_lin/gcc4.4:${MKLROOT}/lib/intel64:${MKLROOT}/include
