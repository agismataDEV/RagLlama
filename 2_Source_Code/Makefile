## Make file version 2
## Should load first ifort with module load intel/2017 and openmpi with module load openmpi-x86_x64
objects = fftw3.o Types.o Constants.o Cisi.o FDWeight.o ParnacTransformDef.o ParnacDataDef.o Stopwatch.o ParnacGeneral.o ParnacParamDef.o ParnacIdentifierDef.o ParnacTransformInit.o ParnacDataSpaceInit.o ParnacOutput.o ParnacDataMap.o  Specialfun.o ContrastGenFilt.o ParnacDataRedist.o f2kcli.o ParnacParamInit.o  ParnacTransformFilter.o ParnacDerivative.o RP.o RKF45.o ODEPACK.o ORDERPACK.o ParnacContrastFunctions.o ParnacBubbleContrast.o ParnacDataSourceInit.o ParnacConvolution.o ParnacIterStepMain.o ParnacSchemes.o ParnacMain.o

LIB_PATH_HDF5 = -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/hdf5_ifort_v1.8.7/lib #-L ~/hdf5_ifort_v1.8.21/lib
INC_PATH_HDF5 = -I /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/hdf5_ifort_v1.8.7/include #-I ~/hdf5_ifort_v1.8.21/include
LIB_PATH_MKL = -L ${MKLROOT}/lib/intel64/
INC_PATH_MKL = -I ${MKLROOT}/include
BIN_PATH_FFTW = -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/bin
LIB_PATH_FFTW = -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/lib
INC_PATH_FFTW = -I /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/include

LIBS_HDF5 = -lhdf5_fortran -lhdf5 -lhdf5_hl -lhdf5hl_fortran  
LIBS_MPI = -lmpifort -lmpi 
LIBS_FFTW3 = $(LIB_PATH_FFTW) $(INC_PATH_FFTW) $(BIN_PATH_FFTW) -lfftw3 -lfftw3_mpi -lfftw3_threads ##-L ~/FFTW_3.3.9/lib/libfftw3.a -L ~/FFTW_3.3.9/libfftw3_mpi.a -lm 
LIBS =  $(LIBS_HDF5) $(LIBS_MPI) $(LIBS_FFTW3) -lz -lm -shared-intel

CFLAGS_MKL = -mkl -liomp5 -lpthread -lm -qopenmp $(INC_PATH_MKL)
CFLAGS_OPT = -O3 -funroll-all-loops #-xSSE4.2 -axAVX,CORE-AVX2 # -fast -ipo
CFLAGS_DEBUG = -g -traceback -debug extended #-check all -fpe0 -warn 
CFLAGS =  -fpp -heap-arrays $(CFLAGS_OPT)  
# For debugging use -g -check all -fpe0 -warn -traceback -debug extended
FFTW_WRAPPERS  = -I fftw_wrappers.inc

## O3 optimization method 
INCS_Agis : ${objects} 
	mpiifort -o INCS_Agis ${objects} $(INC_PATH_HDF5) $(LIB_PATH_HDF5) $(INC_PATH_FFTW) $(LIB_PATH_FFTW) $(LIBS) $(CFLAGS) $(CFLAGS_MKL) 

## External modules
Specialfun.o : Specialfun.f90 
	mpiifort -c Specialfun.f90 $(CFLAGS) 
Constants.o : Constants.f90
	mpiifort -c Constants.f90 $(CFLAGS)
Stopwatch.o : Stopwatch.f90
	mpiifort -c Stopwatch.f90 $(CFLAGS)
Types.o : Types.f90
	mpiifort -c Types.f90 $(CFLAGS)
FDWeight.o : FDWeight.f90
	mpiifort -c FDWeight.f90 $(CFLAGS)
f2kcli.o : f2kcli.f90 
	mpiifort -c f2kcli.f90 $(CFLAGS)
Cisi.o : Cisi.f90
	mpiifort -c Cisi.f90 $(CFLAGS)
fftw3.o : fftw3.f90
	mpiifort -c fftw3.f90 $(CFLAGS)

##Parnac Files
##Definitions
ParnacDataDef.o : ParnacDataDef.f90 
	mpiifort -c ParnacDataDef.f90 $(CFLAGS)
ParnacIdentifierDef.o : ParnacIdentifierDef.f90 
	mpiifort -c ParnacIdentifierDef.f90 $(CFLAGS)
ParnacTransformDef.o : ParnacTransformDef.f90
	mpiifort -c ParnacTransformDef.f90 $(CFLAGS) $(LIBS_FFTW3)
ParnacParamDef.o :ParnacParamDef.f90
	mpiifort -c ParnacParamDef.f90 $(CFLAGS)

## Initializations
ParnacTransformInit.o : ParnacTransformInit.f90
	mpiifort -c ParnacTransformInit.f90 $(CFLAGS) $(LIBS_FFTW3)
ParnacParamInit.o : ParnacParamInit.f90
	mpiifort -c ParnacParamInit.f90 $(CFLAGS) 
ParnacDataSpaceInit.o : ParnacDataSpaceInit.f90
	mpiifort -c ParnacDataSpaceInit.f90 $(CFLAGS)
ParnacDataSourceInit.o : ParnacDataSourceInit.f90
	mpiifort -c ParnacDataSourceInit.f90 $(CFLAGS)

## Main Files
## Data
ParnacDataMap.o : ParnacDataMap.f90
	mpiifort -c ParnacDataMap.f90 $(CFLAGS)
ParnacDataRedist.o : ParnacDataRedist.f90
	mpiifort -c ParnacDataRedist.f90 $(CFLAGS) $(LIBS_FFTW3)

##Convolution
ParnacConvolution.o : ParnacConvolution.f90
	mpiifort -c ParnacConvolution.f90 $(CFLAGS)

##Contrast
ParnacContrastFunctions.o : ParnacContrastFunctions.f90
	mpiifort -c ParnacContrastFunctions.f90 $(CFLAGS)
ParnacBubbleContrast.o : ParnacBubbleContrast.f90
	mpiifort -c ParnacBubbleContrast.f90 $(CFLAGS) $(CFLAGS_MKL)  $(LIBS) 

##Iterative Schemes
ParnacIterStepMain.o  : ParnacIterStepMain.f90
	mpiifort -c ParnacIterStepMain.f90 $(CFLAGS)
ParnacSchemes.o :ParnacSchemes.f90
	mpiifort -c ParnacSchemes.f90 $(CFLAGS)

##Transform and Filter
ParnacTransformFilter.o : ParnacTransformFilter.f90
	mpiifort -c ParnacTransformFilter.f90 $(CFLAGS) $(LIBS_FFTW3)
ContrastGenFilt.o : ContrastGenFilt.f90
	mpiifort -c ContrastGenFilt.f90 $(CFLAGS) $(FFTW_WRAPPERS)

##Miscellaneous
ParnacDerivative.o : ParnacDerivative.f90
	mpiifort -c ParnacDerivative.f90 $(CFLAGS)
ParnacGeneral.o : ParnacGeneral.f90
	mpiifort -c ParnacGeneral.f90 $(CFLAGS)
RP.o : RP.f90
	mpiifort -c RP.f90 $(CFLAGS)
RKF45.o : RKF45.f90
	mpiifort -c RKF45.f90 $(CFLAGS)
ODEPACK.o : ODEPACK.f
	mpiifort -c ODEPACK.f $(CFLAGS)
ORDERPACK.o : ORDERPACK.f90
	mpiifort -c ORDERPACK.f90 $(CFLAGS)

##Output
ParnacOutput.o : ParnacOutput.f90
	mpiifort -c ParnacOutput.f90 -fpp -O3 $(INC_PATH_HDF5)

##Main 
ParnacMain.o :ParnacMain.f90
	mpiifort -c ParnacMain.f90 $(CFLAGS) $(LIBS_FFTW3)
ParnacMain_Compact.o :ParnacMain_Compact.f90
	mpiifort -c ParnacMain_Compact.f90 $(CFLAGS)

move : 
	mv *.o ../3_Compiled_Files/
	mv *.mod ../3_Compiled_Files/
	mv INCS_Agis ../3_Compiled_Files/INCS_Agis_${OUTPUT_FOLDER}
	
move-clean:
	rm ../3_Compiled_Files/*.o 
	rm ../3_Compiled_Files/*.mod 
	rm ../3_Compiled_Files/INCS_Agis_${OUTPUT_FOLDER}
		
clean:
	rm ${objects}
	rm *.mod
	
distclean:
	rm INCS_Agis
	
copy : 
	cp *.o ../3_Compiled_Files/
	cp *.mod ../3_Compiled_Files/
	cp INCS_Agis ../3_Compiled_Files/INCS_Agis_${OUTPUT_FOLDER}
	
	




