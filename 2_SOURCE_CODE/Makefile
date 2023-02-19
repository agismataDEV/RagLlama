## Make file version 2
## Should load first ifort with module load intel/2017 and openmpi with module load openmpi-x86_x64

OBJDIR =  ../3_COMPILED_FILES

objects = 	 $(addprefix $(OBJDIR)/, fftw3.o Types.o Constants.o Cisi.o FDWeight.o ParnacTransformDef.o ParnacDataDef.o \
								  Stopwatch.o ParnacGeneral.o ParnacParamDef.o ParnacIdentifierDef.o ParnacTransformInit.o \
								  ParnacDataSpaceInit.o ParnacOutput.o ParnacDataMap.o  Specialfun.o ContrastGenFilt.o \
								  ParnacDataRedist.o f2kcli.o ParnacParamInit.o  ParnacTransformFilter.o ParnacDerivative.o \
			 					  ORDERPACK.o ODEPACK.o RP.o RKF45.o ParnacContrastFunctions.o ParnacPointSourceCloud.o \
								  ParnacDataSourceInit.o ParnacConvolution.o ParnacIterStepMain.o ParnacMain.o\
								)  
 
LIB_PATH_HDF5 	= -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/hdf5_ifort_v1.8.7/lib #-L ~/hdf5_ifort_v1.8.21/lib
INC_PATH_HDF5 	= -I /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/hdf5_ifort_v1.8.7/include #-I ~/hdf5_ifort_v1.8.21/include
LIB_PATH_MKL 	= -L ${MKLROOT}/lib/intel64/
INC_PATH_MKL 	= -I ${MKLROOT}/include
BIN_PATH_FFTW 	= -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/bin
LIB_PATH_FFTW 	= -L /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/lib
INC_PATH_FFTW 	= -I /tudelft.net/staff-bulk/tnw/IST/AK/hpc/agismatalliota/FFTW_3.3.9/include

LIBS_HDF5 	= 	-lhdf5_fortran -lhdf5 -lhdf5_hl -lhdf5hl_fortran  
LIBS_MPI 	= 	-lmpifort -lmpi
LIBS_FFTW3 	= 	$(LIB_PATH_FFTW) $(INC_PATH_FFTW) $(BIN_PATH_FFTW) -lfftw3 ##-lfftw3_mpi -lfftw3_threads ##-L ~/FFTW_3.3.9/lib/libfftw3.a -L ~/FFTW_3.3.9/libfftw3_mpi.a -lm 
LIBS		=  	$(LIBS_HDF5) $(LIBS_MPI) $(LIBS_FFTW3) -lz -lm -shared-intel

CFLAGS_MKL 		= 	-mkl -liomp5 -lpthread -lm -qopenmp $(INC_PATH_MKL)
CFLAGS_OPT 		= 	-O3 -funroll-all-loops #-xSSE4.2 -axAVX,CORE-AVX2 # -fast -ipo
CFLAGS_DEBUG 	= 	-g -traceback -debug extended #-check all -fpe0 -warn 
CFLAGS 			=  	-fpp -heap-arrays $(CFLAGS_OPT) 

# For debugging use -g -check all -fpe0 -warn -traceback -debug extended
FFTW_WRAPPERS   = 	-I fftw_wrappers.inc

## O3 optimization method 
$(OBJDIR)/INCS_Agis : ${objects} 
	mpiifort -o $(OBJDIR)/INCS_Agis ${objects} $(INC_PATH_HDF5) $(LIB_PATH_HDF5) $(INC_PATH_FFTW) $(LIB_PATH_FFTW) $(LIBS) $(CFLAGS) $(CFLAGS_MKL) 

## External modules
$(OBJDIR)/Specialfun.o : Specialfun.f90 
	mpiifort -c Specialfun.f90 $(CFLAGS)  -o $(OBJDIR)/Specialfun.o -module $(OBJDIR)
$(OBJDIR)/Constants.o : Constants.f90
	mpiifort -c Constants.f90 $(CFLAGS) -o $(OBJDIR)/Constants.o -module $(OBJDIR)
$(OBJDIR)/Stopwatch.o : Stopwatch.f90
	mpiifort -c Stopwatch.f90 $(CFLAGS) -o $(OBJDIR)/Stopwatch.o -module $(OBJDIR)
$(OBJDIR)/Types.o : Types.f90
	mpiifort -c Types.f90 $(CFLAGS) -o $(OBJDIR)/Types.o -module $(OBJDIR)
$(OBJDIR)/FDWeight.o : FDWeight.f90
	mpiifort -c FDWeight.f90 $(CFLAGS) -o $(OBJDIR)/FDWeight.o -module $(OBJDIR)
$(OBJDIR)/f2kcli.o : f2kcli.f90 
	mpiifort -c f2kcli.f90 $(CFLAGS) -o $(OBJDIR)/f2kcli.o -module $(OBJDIR)
$(OBJDIR)/Cisi.o : Cisi.f90
	mpiifort -c Cisi.f90 $(CFLAGS) -o $(OBJDIR)/Cisi.o -module $(OBJDIR)
$(OBJDIR)/fftw3.o : fftw3.f90
	mpiifort -c fftw3.f90 $(CFLAGS) -o $(OBJDIR)/fftw3.o -module $(OBJDIR)

##Parnac Files
##Definitions
$(OBJDIR)/ParnacDataDef.o : ParnacDataDef.f90 
	mpiifort -c ParnacDataDef.f90 $(CFLAGS) -o $(OBJDIR)/ParnacDataDef.o -module $(OBJDIR)
$(OBJDIR)/ParnacIdentifierDef.o : ParnacIdentifierDef.f90 
	mpiifort -c ParnacIdentifierDef.f90 $(CFLAGS) -o $(OBJDIR)/ParnacIdentifierDef.o -module $(OBJDIR)
$(OBJDIR)/ParnacTransformDef.o : ParnacTransformDef.f90
	mpiifort -c ParnacTransformDef.f90 $(CFLAGS) $(LIBS_FFTW3) -o $(OBJDIR)/ParnacTransformDef.o -module $(OBJDIR)
$(OBJDIR)/ParnacParamDef.o :ParnacParamDef.f90
	mpiifort -c ParnacParamDef.f90 $(CFLAGS) -o $(OBJDIR)/ParnacParamDef.o -module $(OBJDIR)

## Initializations
$(OBJDIR)/ParnacTransformInit.o : ParnacTransformInit.f90
	mpiifort -c ParnacTransformInit.f90 $(CFLAGS) $(LIBS_FFTW3) -o $(OBJDIR)/ParnacTransformInit.o -module $(OBJDIR)
$(OBJDIR)/ParnacParamInit.o : ParnacParamInit.f90
	mpiifort -c ParnacParamInit.f90 $(CFLAGS)  -o $(OBJDIR)/ParnacParamInit.o -module $(OBJDIR)
$(OBJDIR)/ParnacDataSpaceInit.o : ParnacDataSpaceInit.f90
	mpiifort -c ParnacDataSpaceInit.f90 $(CFLAGS) -o $(OBJDIR)/ParnacDataSpaceInit.o -module $(OBJDIR)
$(OBJDIR)/ParnacDataSourceInit.o : ParnacDataSourceInit.f90
	mpiifort -c ParnacDataSourceInit.f90 $(CFLAGS) -o $(OBJDIR)/ParnacDataSourceInit.o -module $(OBJDIR)

## Main Files
## Data
$(OBJDIR)/ParnacDataMap.o : ParnacDataMap.f90
	mpiifort -c ParnacDataMap.f90 $(CFLAGS) -o $(OBJDIR)/ParnacDataMap.o -module $(OBJDIR)
$(OBJDIR)/ParnacDataRedist.o : ParnacDataRedist.f90
	mpiifort -c ParnacDataRedist.f90 $(CFLAGS) $(LIBS_FFTW3) -o $(OBJDIR)/ParnacDataRedist.o -module $(OBJDIR)

##Convolution
$(OBJDIR)/ParnacConvolution.o : ParnacConvolution.f90
	mpiifort -c ParnacConvolution.f90 $(CFLAGS) -o $(OBJDIR)/ParnacConvolution.o -module $(OBJDIR)

##Contrast
$(OBJDIR)/ParnacContrastFunctions.o : ParnacContrastFunctions.f90
	mpiifort -c ParnacContrastFunctions.f90 $(CFLAGS) -o $(OBJDIR)/ParnacContrastFunctions.o -module $(OBJDIR)
$(OBJDIR)/ParnacPointSourceCloud.o : ParnacPointSourceCloud.f90
	mpiifort -c ParnacPointSourceCloud.f90 $(CFLAGS) $(CFLAGS_MKL)  $(LIBS)  -o $(OBJDIR)/ParnacPointSourceCloud.o -module $(OBJDIR)

##Iterative Schemes
$(OBJDIR)/ParnacIterStepMain.o  : ParnacIterStepMain.f90
	mpiifort -c ParnacIterStepMain.f90 $(CFLAGS) -o $(OBJDIR)/ParnacIterStepMain.o -module $(OBJDIR)
# $(OBJDIR)/ParnacSchemes.o :ParnacSchemes.f90
# 	mpiifort -c ParnacSchemes.f90 $(CFLAGS) -o $(OBJDIR)/ParnacSchemes.o -module $(OBJDIR)

##Transform and Filter
$(OBJDIR)/ParnacTransformFilter.o : ParnacTransformFilter.f90
	mpiifort -c ParnacTransformFilter.f90 $(CFLAGS) $(LIBS_FFTW3) -o $(OBJDIR)/ParnacTransformFilter.o -module $(OBJDIR)
$(OBJDIR)/ContrastGenFilt.o : ContrastGenFilt.f90
	mpiifort -c ContrastGenFilt.f90 $(CFLAGS) $(FFTW_WRAPPERS) -o $(OBJDIR)/ContrastGenFilt.o -module $(OBJDIR)

##Miscellaneous
$(OBJDIR)/ParnacDerivative.o : ParnacDerivative.f90
	mpiifort -c ParnacDerivative.f90 $(CFLAGS) -o $(OBJDIR)/ParnacDerivative.o -module $(OBJDIR)
$(OBJDIR)/ParnacGeneral.o : ParnacGeneral.f90
	mpiifort -c ParnacGeneral.f90 $(CFLAGS) -o $(OBJDIR)/ParnacGeneral.o -module $(OBJDIR)
$(OBJDIR)/RP.o : RP.f90
	mpiifort -c RP.f90 $(CFLAGS) -o $(OBJDIR)/RP.o -module $(OBJDIR)
$(OBJDIR)/RP_SP.o : RP_SP.f90
	mpiifort -c RP_SP.f90 $(CFLAGS) -o $(OBJDIR)/RP_SP.o -module $(OBJDIR)
$(OBJDIR)/RKF45.o : RKF45.f90
	mpiifort -c RKF45.f90 $(CFLAGS) -o $(OBJDIR)/RKF45.o -module $(OBJDIR)
$(OBJDIR)/ODEPACK.o : ODEPACK.f
	mpiifort -c ODEPACK.f $(CFLAGS) -o $(OBJDIR)/ODEPACK.o -module $(OBJDIR)
$(OBJDIR)/ODEPACK_SP.o : ODEPACK_SP.f
	mpiifort -c ODEPACK_SP.f $(CFLAGS) -o $(OBJDIR)/ODEPACK_SP.o -module $(OBJDIR)
$(OBJDIR)/ORDERPACK.o : ORDERPACK.f90
	mpiifort -c ORDERPACK.f90 $(CFLAGS) -o $(OBJDIR)/ORDERPACK.o -module $(OBJDIR)

##Output
$(OBJDIR)/ParnacOutput.o : ParnacOutput.f90
	mpiifort -c ParnacOutput.f90 -fpp -O3 $(INC_PATH_HDF5) -o $(OBJDIR)/ParnacOutput.o -module $(OBJDIR)

##Main 
$(OBJDIR)/ParnacMain.o :ParnacMain.f90
	mpiifort -c ParnacMain.f90 $(CFLAGS) $(LIBS_FFTW3) -o $(OBJDIR)/ParnacMain.o -module $(OBJDIR)
# $(OBJDIR)/ParnacMain_Compact.o :ParnacMain_Compact.f90
# 	mpiifort -c ParnacMain_Compact.f90 $(CFLAGS) -o $(OBJDIR)/ParnacMain_Compact.o -module $(OBJDIR)
		
clean:
	rm $(OBJDIR)/*.o
	rm $(OBJDIR)/*.mod

remove:
    @read -p "Enter File Name: " OBJRM; 
	rm $(OBJDIR)/$OBJRM.o 
	rm $(OBJDIR)/$OBJRM.mod
	
distclean:
	rm $(OBJDIR)/INCS_Agis
	
	




