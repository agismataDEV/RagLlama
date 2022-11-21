    include "mkl_vsl.f90"
	MODULE ParnacParamInit

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen / Jasper de Koning
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    !   Copyright (C)   Laboratory of Electromagnetic Research
    !                   Delft University of Technology, Delft, The Netherlands
    !
    !   Permission to copy or distribute this software or documentation
    !   in hard copy or soft copy granted only by written license
    !   obtained from the Laboratory of Electromagnetic Research.
    !   All rights reserved. No parts of this publication may be
    !   reproduced, stored in a retrieval system (e.g. in memory, disk,
    !   or core) or be transmitted by any means, electronic, mechanical,
    !   photocopy, recording, or otherwise, without written permission
    !   from the publisher.
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The module ParnacParamInit contains functions and subroutines employed in
    !   setting the program's i/o directory and file names, and the source, model
    !   and domain parameters.
    !
    ! *****************************************************************************
    !
    !   MODULES USED
    !
	!include "mkl.fi"
    USE Types
    USE Constants
    USE ParnacGeneral
    USE ParnacParamDef
    USE ParnacDataDef
    USE ParnacIdentifierDef

    USE f2kcli, ONLY : COMMAND_ARGUMENT_COUNT, GET_COMMAND_ARGUMENT
    USE MKL_VSL
	USE MKL_VSL_TYPE
    ! *****************************************************************************
    !
    !   GLOBAL DECLARATIONS
    !
    !   none
    !
    IMPLICIT NONE

    ! *****************************************************************************
    !
    !   CONTAINED SUBROUTINES AND FUNCTIONS
    !
    !   ReadCommandLineAndEnvironment   sub   checks and evaluates the command line
    !                                         arguments and the environment
    !                                         variables
    !   ReadConfigFile   sub   Read the configuration file, set all parameters
    !   ConfigFileInfo   sub   Echo the parameters loaded from the config file
    !                            to screen
    !   SpecialChecksConfigParameters   sub   Run a number of checks and
    !                                         corrections on the parameters
    !   NormalizeConfigParameters  sub   Normalize all parameters to the normalized
    !                                    dimensions as employed in the program
    !   BINARY_read1Darray         sub   Read a 1D array from a binary file
    !   BINARY_read1Darraysize     fun   Read the size of a 1D array stored in a
    !                                    binary file
    !   BINARY_read3Darray         sub   Read a 3D array from a binary file
    !   BINARY_read3Darraysize     fun   Read the sizes of a 3D array stored in a
    !                                    binary file
    !   BINARY_read3DarrayTdelay   fun   Read the maximum time delay of the delay
    !                                    array stored in a binary file
    !   BINARY_save3Darray_test    sub   Write two numbers to a file - test for
    !                                    determining the Fortran binary format
    !   ASCII_read1Darray          sub   Read a 1D array from an ASCII file
    !   ASCII_read1Darraysize      fun   Read the size of a 1D array stored in an
    !                                    ASCII file
    !   ASCII_read3Darray          sub   Read a 3D array from an ASCII file
    !   ASCII_read3Darraysize      fun   Read the sizes of a 3D array stored in an
    !                                    ASCII file
    !   remove_CR                  sub   replaces all carriage returns in a string
    !                                    by spaces
    !
    ! =============================================================================
    !
    CONTAINS

    SUBROUTINE ReadCommandLineAndEnvironment(iProcID)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen / Jasper de Koning
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine ReadCommandLineAndEnvironment reads the arguments from the
    !   command line call to ParnacMain, and it sets the global variables for the
    !   inputfilename sInputName, the root of the output filenames sOutputRoot,
    !   and the input/output directories sInputDir and sOutputDir either to the
    !   requested or to default values. For the temporary directory name in
    !   sTempDir, the environment variable TMPDIR is checked. The configuration
    !   file and the directories are checked for existence.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   iProcID   i   i8b   ID number of the current processor
    !
    integer(i8b), intent(in) :: iProcID

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   i        i8b    temporary variable
    !   err      i8b    error number
    !   error    lgt    error flag
    !   acTemp   char   Temporary char array for output log messages
    !   CMD      char   arguments of the command line
    !   NARG     int    number of arguments
    !
    integer(i8b) ::                 i,err
    logical(lgt) ::					error
    character(LEN = 2048)::                acTemp;
    CHARACTER(LEN = 2048)  :: CMD
    INTEGER            :: NARG


    ! *****************************************************************************
    !
    !   I/O
    !
    !   Command line and environment variable input, checking for the existence of
    !   the configuration file and directories
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   COMMAND_ARGUMENT_COUNT
    !   GET_COMMAND_ARGUMENT
    !   getenv
    !
    ! =============================================================================

    CMD = ""
    NARG = COMMAND_ARGUMENT_COUNT();
    if (NARG /=0) then
        call GET_COMMAND_ARGUMENT(1,CMD)
    end if
    if ((NARG ==0) .or. &
        (trim(CMD) .eq. "-h") .or. (trim(CMD) .eq. "-help")  .or. &
        (trim(CMD) .eq. "-?") .or. (trim(CMD) .eq. "--help") .or. &
        (trim(CMD) .eq. "/?") .or. (trim(CMD) .eq. "/h")     .or. &
        (trim(CMD) .eq. "/help") ) then
    ! Print help information and exit program
    print *, " "
    print *, "Iterative solver for full-wave, nonlinear, pulsed pressure fields in a"
    print *, " three-dimensional domain from sources with an arbitrary plane aperture."
    print *, " "
    print *, "Usage: ParnacMain configfile [outputroot] [inputdir] [outputdir]"
    print *, " "
    print *, "  configfile   name of the configuration file that contains all settings"
    print *, "  outputroot   root name of the output files. If not supplied, then the"
    print *, "               inputfilename without extension is used."
    print *, "  inputdir     directory where all input files are located. If not supplied,"
    print *, "               then a standard directory is used."
    print *, "  outputdir    directory where all output files are stored. If not supplied,"
    print *, "               then a standard directory is used."
    print *, " "
    print *, "ParnacMain expects the environment variable TMPDIR to contain the name of the"
    print *, "temporary directory. If not supplied, then the current directory is used."
    stop
    else
        ! First argument is inputfile
        sInputName = CMD
    end if
    
    if (NARG > 1) then
        call GET_COMMAND_ARGUMENT(2,CMD)
        sOutputRoot = CMD
    else
        i=index(sInputName,'.')
        if (i==0) then
            sOutputRoot = trim(sInputName)
        else
            sOutputRoot = sInputName(1:i-1)
        end if
    end if

    if (NARG > 2) then
        call GET_COMMAND_ARGUMENT(3,CMD)
        sInputDir = CMD
    else
        sInputDir = sDefInputDir
    end if
    ! If necessary, add slash to input dir
    i=len_trim(sInputDir)
    if(sInputDir(i:i)/=sSystemSlash) then
        sInputDir = trim(sInputDir) // sSystemSlash
    end if

    if (NARG > 3) then
        call GET_COMMAND_ARGUMENT(4,CMD)
        sOutputDir = CMD
    else
        sOutputDir = sDefOutputDir
    end if
    ! If necessary, add slash to output dir
    i=len_trim(sOutputDir)
    if(sOutputDir(i:i)/=sSystemSlash) then
        sOutputDir = trim(sOutputDir) // sSystemSlash
    end if

    ! Get the environment variable TMPDIR that should contain the name of the
    ! temporary file directory. If not present, sTempDir becomes sDefTempDir
    call getenv("TMPDIR",sTempDir)
    if(trim(sTempDir)=="") then
        sTempDir = sDefTempDir
    end if
    ! If necessary, add slash to temp dir
    i=len_trim(sTempDir)
    if(sTempDir(i:i)/=sSystemSlash) then
        sTempDir = trim(sTempDir) // sSystemSlash
    end if

    if (iProcID==0) then
        !Check existence of input, output and temp directories
        open(unit=iImportUNIT,file=trim(sInputDir)//"test.tst",iostat=err)
        if (err/=0) then
            write(*, "('Error, input directory ',A, ' does not exist')") trim(sInputDir)
            stop 'Error, input directory does not exist'
        else
            close(unit=iImportUNIT,status="DELETE",iostat = err)
        end if
        open(unit=iImportUNIT,file=trim(sOutputDir)//"test.tst",iostat=err)
        if (err/=0) then
            write(*, "('Error, output directory ',A, ' does not exist')") trim(sOutputDir)
            stop 'Error, output directory does not exist'
        else
            close(unit=iImportUNIT,status="DELETE",iostat = err)
        end if
        ! Check the existence of the Temporary directory
        open(unit=iImportUNIT,file=trim(sTempDir)//"test.tst",iostat=err)
        if (err/=0) then
            write(*, "('Error, temp directory ',A, ' does not exist')") trim(sTempDir)
            stop 'Error, temp directory does not exist'
        else
            close(unit=iImportUNIT,status="DELETE",iostat = err)
        end if

        ! Check the existence of the inputfile in the input dir
        inquire(FILE=trim(sInputDir)//trim(sInputName),EXIST=error)
        if(error .EQV. .FALSE.) then
            write(*, "('Error, cannot find inputfile',A)") trim(sInputDir)//trim(sInputName)
            stop 'Error, cannot find inputfile'
        end if
    end if
    
    END SUBROUTINE ReadCommandLineAndEnvironment

    SUBROUTINE ReadConfigFile(inputfileversion,errstatus)
    ! =============================================================================
    !
    !   Programmer: Koos Huijssen / Jasper de Koning
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine ReadConfigFile reads the configuration file and sets the
    !   parameters in the parameter structures that define the configuration,
    !   source, contrast, domain and that include various model settings.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   inputfileversion   io   dp   Version of the configuration file that is
    !                                needed for the current version of the program.
    !                                In the RunEstimator program, it returns the
    !                                real version of the config file.
    !   errstatus          o   i8b   error status, is zero if all is okay nonzero
    !                                if there are one or more errors.
    !
    real(dp), intent(inout) :: inputfileversion
    integer(i8b), intent(out) :: errstatus

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   readstatus   i8b    status of each separate read command
    !   i            i8b    temporary variable
    !   fileversion  i8b    version of the config file
    !   acTemp       char   Temporary char array for output log messages
    !
    integer(i8b) ::                 readstatus
    integer(i8b) ::                 i
    real(dp) ::						fileversion
    character(LEN = 1024)::         acTemp;
	
    TYPE (VSL_STREAM_STATE) :: STREAM2

    ! *****************************************************************************
    !
    !   I/O
    !
    !   Reads the configuration file
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   remove_CR
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=trim(sInputDir)//trim(sInputName),form="FORMATTED",status="OLD",iostat=errstatus)

    !check file version
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) fileversion
    errstatus = errstatus + readstatus
    if(abs(fileversion-inputfileversion)>1e-3) then
#ifndef RUNESTIMATOR
        write(*,"('Config file version number is',F5.2,', should be',F5.2)") fileversion,inputfileversion
        stop
        !#else
        !				errstatus=1
        !				inputfileversion=fileversion
        !		        close(unit=iImportUNIT)
        !				return
#endif
    end if

    !read medium parameters
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.rho0
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.c0
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams%mu
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams%P0
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.fc0
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.beta
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.a_att
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cMediumParams.b_att
    errstatus = errstatus + readstatus

    !read domain parameters
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Lx
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Ly
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Lt
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Sz
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Lz
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Thetax
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Thetay
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Thetat
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.freq0
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Fnyq
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams%PPW 
    errstatus = errstatus + readstatus
    
    !read other model parameters
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) acTemp
    errstatus = errstatus + readstatus
    call remove_CR(acTemp)
    if(trim(acTemp)=='qsource') then
        cModelParams.PrimarySourceType=iEI_QSOURCE	!Q source
    else if(trim(acTemp)=='fsource') then
        cModelParams.PrimarySourceType=iEI_FSOURCE	!F source
    else if(trim(acTemp)=='dtqsource') then
        cModelParams.PrimarySourceType=iEI_DTQSOURCE	!dQ/dt source from measured signature or in FieldII fashion
    else if(trim(acTemp)=='planewave') then
        cModelParams.PrimarySourceType=iEI_PLANEWAVE	!Plane wave primary field solutionelse if(trim(acTemp)=='planewave') then
    else if(trim(acTemp)=='pointsourcecloud') then
        cModelParams.PrimarySourceType=iEI_POINTSOURCECLOUD	!Point Source Cloud
        read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%N
		ALLOCATE(PointSourceCloudParams%R0(PointSourceCloudParams%N))
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%PointSourceAmplitude
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%Dist_Amplitude
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%R0
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%ClusterDimsRatio(:,1)
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%ClusterDimsRatio(:,2)
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%ClusterDimsRatio(:,3)
		read(iImportUNIT,*,IOSTAT = readstatus) PointSourceCloudParams%MinInBetweenDist
    else if(trim(acTemp)=='filename') then
		read(iImportUNIT,*)
        read(iImportUNIT,*,IOSTAT = readstatus) acTemp
        errstatus = errstatus + readstatus
        call remove_CR(acTemp)
        cModelParams.PrimarySourceType=iEI_LOADFIELD	
		cModelParams.FieldFilename = acTemp
		do i = 1,3
            read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.xyzfielddim(i,:)
		enddo
    else
        cModelParams.PrimarySourceType=0 ! Doesn't exist
    end if
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) acTemp
    errstatus = errstatus + readstatus
    call remove_CR(acTemp)
    if (trim(acTemp)=='nonlin') then
        cModelParams.ContrastSourceType=iCI_NONLIN
    else if (trim(acTemp)=='complexcontrast') then
        cModelParams.ContrastSourceType=iCI_COMPLEXCONTRAST
    else if (trim(acTemp)=='luneberg') then
        cModelParams.ContrastSourceType=iCI_LUNEBERG
    else if (trim(acTemp)=='sphere') then
        cModelParams.ContrastSourceType=iCI_SPHERE
    else if (trim(acTemp)=='blob') then
        cModelParams.ContrastSourceType=iCI_BLOB
    else if (trim(acTemp)=='scatterer') then
        cModelParams.ContrastSourceType=iCI_SCATTERER
    else
        cModelParams.ContrastSourceType=0 ! Doesn't exist
    end if
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.FDTorder
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.FDXorder
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) iDebugLvl
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) i
    if(i==0) then
        cModelParams.UseAntiAliasing = .false.
    else
        cModelParams.UseAntiAliasing = .true.
    end if
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) i
    if(i==0) then
        cModelParams.UseSupportTapering=.false.
    else
        cModelParams.UseSupportTapering=.true.
    end if
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) i
    if(i==0) then
        cModelParams.UseFreqTapering=.false.
    else
        cModelParams.UseFreqTapering=.true.
    end if
    errstatus = errstatus + readstatus
    read(iImportUNIT,*,IOSTAT = readstatus) i 
    if (i == 0) then
    	cModelParams.UseYSymmetry = .false.
    else
    	cModelParams.UseYSymmetry = .true.
    endif   
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.Numbeams
    errstatus = errstatus + readstatus
    if(cModelParams%Numbeams>0) then
        allocate(cModelParams%numiterations(1:cModelParams%Numbeams))
        do i=1,cModelParams%Numbeams
            read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.numiterations(i)
            errstatus = errstatus + readstatus
        end do
    end if
    read(iImportUNIT,*)
    read(iImportUNIT,*,IOSTAT = readstatus) acTemp
    errstatus = errstatus + readstatus
    call remove_CR(acTemp)
    if (trim(acTemp)=='firstandlast') then
        cModelParams%Slicesavespecifier=iAI_FIRSTANDLAST	!Save first and last iterates
    else if (trim(acTemp)=='last') then
        cModelParams%Slicesavespecifier=iAI_LAST	!Save last iterate only
    else if (trim(acTemp)=='all') then
        cModelParams%Slicesavespecifier=iAI_ALL	!Save all iterates
    else
        cModelParams%Slicesavespecifier=0	!Doesn't exist
    end if
    read(iImportUNIT,*,IOSTAT = readstatus) cModelParams%Numslices
    errstatus = errstatus + readstatus
    if(cModelParams%Numslices>0) then
        allocate(cModelParams%xyzslicedim(1:cModelParams%Numslices))
        allocate(cModelParams%xyzslicepos(1:cModelParams%Numslices))
        allocate(cModelParams%xyzslicebeam(1:cModelParams%Numslices))
        allocate(cModelParams%xyzsliceindex(1:cModelParams%Numslices))
        do i=1,cModelParams%Numslices
            read(iImportUNIT,*,IOSTAT = readstatus) cModelParams.xyzslicedim(i),&
                cModelParams.xyzslicepos(i)
            errstatus = errstatus + readstatus
        end do
    end if
	
	if (cModelParams.PrimarySourceType/=iEI_LOADFIELD) then
		!read source input filename
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.srcfilename
		call remove_CR(cSourceParams.srcfilename)
		errstatus = errstatus + readstatus

		if(trim(cSourceParams.srcfilename)=='none') then
			! read source signature parameters
			read(iImportUNIT,*)
			read(iImportUNIT,*,IOSTAT = readstatus) acTemp
			call remove_CR(acTemp)
			errstatus = errstatus + readstatus
			if (trim(acTemp)=='gaussian') then
				cSourceParams.srcsigntype = iGI_GAUSSIAN

				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Pstart
				errstatus = errstatus + readstatus
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Tpulse
				errstatus = errstatus + readstatus
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Twidth
				errstatus = errstatus + readstatus
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Tdelay
				errstatus = errstatus + readstatus
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Power
				errstatus = errstatus + readstatus
			else if(trim(acTemp) == 'blackman') then
				cSourceParams.srcsigntype = iGI_BLACKMAN

				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Pstart
				errstatus = errstatus + readstatus
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Tpulse
				errstatus = errstatus + readstatus
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.Tdelay
				errstatus = errstatus + readstatus
			else
				! otherwise, assume that we have a filename
				cSourceParams.srcsigntype = iGI_FILE
				cSourceParams.srcsignfilename = acTemp
			end if

			if (cModelParams.PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams.PrimarySourceType/=iEI_POINTSOURCECLOUD) then
				!read source shape parameters
				read(iImportUNIT,*)
				read(iImportUNIT,*,IOSTAT = readstatus) acTemp
				call remove_CR(acTemp)
				errstatus = errstatus + readstatus
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				read(iImportUNIT,*)
				if (trim(acTemp)=='cylindrical') then
					cSourceParams%srcshapetype = iHI_CYLINDRICAL
					!print *, "Acciaaaaaaa"
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.radius
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.radfocus
					errstatus = errstatus + readstatus
				else if (trim(acTemp)=='rectangular') then
					cSourceParams%srcshapetype = iHI_RECTANGULAR

					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.rectwidth
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.rectheight
					errstatus = errstatus + readstatus
				else if (trim(acTemp)=='triangular') then
					cSourceParams%srcshapetype = iHI_TRIANGULAR

					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.triangedge
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.triangxcenter
					errstatus = errstatus + readstatus
				else if (trim(acTemp)=='pointsource') then
					cSourceParams%srcshapetype = iHI_POINTSOURCE

				else if (trim(acTemp)=='phasedarray') then
					cSourceParams%srcshapetype = iHI_PHASEDARRAY

					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.numel
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.elwidth
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.elheight
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.kerf
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.focusx
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.focusz
					errstatus = errstatus + readstatus
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.elevationfocusz
					errstatus = errstatus + readstatus
					read(iImportUNIT,*)
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.td_filename
					errstatus = errstatus + readstatus
					call remove_CR(cSourceParams.td_filename)
					read(iImportUNIT,*)
					read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.phaseapodfilename
					errstatus = errstatus + readstatus
					call remove_CR(cSourceParams.phaseapodfilename)
					!KH					else if (trim(acTemp)=='matrixarray') then
					!KH						cSourceParams%srcshapetype = iHI_MATRIXARRAY
					!KH
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matnumelx
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matnumely
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matelwidth
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matelheight
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matkerfx
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matkerfy
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matfocusx
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matfocusy
					!KH						errstatus = errstatus + readstatus
					!KH						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.matfocusz
					!KH						errstatus = errstatus + readstatus
					!KH!				        read(iImportUNIT,*)
					!KH!						read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.phaseapodfilename
					!KH!						errstatus = errstatus + readstatus
					!KH!						call remove_CR(cSourceParams.phaseapodfilename)
				else
					cSourceParams%srcshapetype = 0	!Doesn't exist
				end if
			end if
		end if 
		! L.d. 09-05-2012 Extension of the input informations
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.NKC
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.LCSM
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BiCGSTAB
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.CG
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.LINEARCASEIMPLEMENTATION
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.DFM
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SHIFT
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.alpha1
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.b1
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SOS1
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.KAPPA1
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BETA1
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.alpha2
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.b2
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SOS2
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.KAPPA2
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BETA2
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.alpha3
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.b3
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SOS3
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.KAPPA3
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BETA3
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.alpha4
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.b4
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SOS4
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.KAPPA4
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BETA4
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.alpha5
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.b5
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.SOS5
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.KAPPA5
		read(iImportUNIT,*,IOSTAT = readstatus) cSourceParams.BETA5
		if (cModelParams.ContrastSourceType==iCI_SCATTERER) then
            read(iImportUNIT,*)
            read(iImportUNIT,*)
            read(iImportUNIT,*)
		
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ScType
			call remove_CR(ScattererParams%ScType)
			if (trim(ScattererParams%ScType) == 'microbubble') then
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%kappa_s
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%sigma_w
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%sigma_R0
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%gama
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%chi
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%Solver_Method
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%Solver_Normalize
			elseif (trim(ScattererParams%ScType) == 'linear') then
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%rho1
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%c1
			endif	
	    	
			read(iImportUNIT,*)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%N
			ALLOCATE(ScattererParams%R0(ScattererParams%N))
			read(iImportUNIT,*)
			read(iImportUNIT,*)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%Distribution
			call remove_CR(ScattererParams%Distribution)
			if (trim(ScattererParams%Distribution) == 'monodisperse') then
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%R0(1)
				ScattererParams%R0(:) = ScattererParams%R0(1)
			elseif (trim(ScattererParams%Distribution) == 'polydisperse') then
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%PDRange
				!CALL RANDOM_NUMBER(ScattererParams%R0)
				!ScattererParams%R0 = ScattererParams%R0*ScattererParams%PDRange(2) + (1-ScattererParams%R0)*ScattererParams%PDRange(1)
				errstatus = vslnewstream( STREAM2, VSL_BRNG_MCG31, 1 )
				errstatus = vdrnggamma( VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, STREAM2, INT(ScattererParams%N), ScattererParams%R0 , 1.2D0, 0.0D0, 1.0D0 )
				ScattererParams%R0 = (ScattererParams%R0-MINVAL(ScattererParams%R0))
				ScattererParams%R0 = ScattererParams%R0/MAXVAL(ScattererParams%R0)*(ScattererParams%PDRange(2)-ScattererParams%PDRange(1))+ScattererParams%PDRange(1)
			endif
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ClusterDimsRatio(:,1)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ClusterDimsRatio(:,2)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ClusterDimsRatio(:,3)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%MinInBetweenDist
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ClusterSlicesN
				ALLOCATE(ScattererParams%ClusterSliceDim(ScattererParams%ClusterSlicesN))
			if (ScattererParams%ClusterSlicesN > 0) then
				read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%ClusterSliceDim
			endif
			read(iImportUNIT,*)
			read(iImportUNIT,*,IOSTAT = readstatus) ScattererParams%GridPointsPressure
			call remove_CR(ScattererParams%GridPointsPressure)
	 	endif
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*)
		read(iImportUNIT,*,IOSTAT = readstatus) acTemp
		call remove_CR(acTemp)
		errstatus = errstatus + readstatus
		cSourceParams.residualdirectory = trim(sOutputDir)//acTemp !'\\srv256\NASSHARE\ldemi\DATA\OFDM Second Round\MixingComponents\'

    end if

    close(unit=iImportUNIT)

    END SUBROUTINE ReadConfigFile

    SUBROUTINE ConfigFileInfo(fileversion)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine ConfigFileInfo echoes the configuration parameters, as
    !   obtained from the configuration file, to the screen.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   fileversion   i   dp   Version of the configuration file.
    !
    real(dp), intent(in) :: fileversion

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   i            i8b    temporary variable
    !
    integer(i8b) :: i

    ! *****************************************************************************
    !
    !   I/O
    !
    !   Echoes the configuration file parameters to screen
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    write ( *, '("Configuration input file: ",A)') trim(sInputDir)//trim(sInputName)
    write ( *, '(" ")')
    write ( *, '("---------------------------------------------------------------------------")')
    write ( *, '("Input file, version ",F10.2)') fileversion
    write ( *, '("---------------------------------------------------------------------------")')
    write ( *, '(" ")')
    write ( *, '("rho0:  ",E10.3)') cMediumParams.rho0
    write ( *, '("c0:    ",E10.3)') cMediumParams.c0
    write ( *, '("mu:    ",E10.3)') cMediumParams.mu
    write ( *, '("P0:    ",E10.3)') cMediumParams.P0
    write ( *, '("fc0:   ",E10.3)') cMediumParams.fc0
    write ( *, '("beta:  ",E10.3)') cMediumParams.beta
    write ( *, '("a_att: ",E10.3)') cMediumParams.a_att
    write ( *, '("b_att: ",E10.3)') cMediumParams.b_att
    write ( *, '(" ")')
    write ( *, '("Lx:     ",E10.3)') cModelParams.Lx
    write ( *, '("Ly:     ",E10.3)') cModelParams.Ly
    write ( *, '("Lt:     ",E10.3)') cModelParams.Lt
    write ( *, '("Sz:     ",E10.3)') cModelParams.Sz
    write ( *, '("Lz:     ",E10.3)') cModelParams.Lz
    write ( *, '("ThetaX: ",E10.3)') cModelParams.Thetax
    write ( *, '("ThetaY: ",E10.3)') cModelParams.Thetay
    write ( *, '("ThetaT: ",E10.3)') cModelParams.Thetat
    write ( *, '("freq0:  ",E10.3)') cModelParams.freq0
    write ( *, '("Fnyq:   ",E10.3)') cModelParams.Fnyq
    write ( *, '(" ")')
    write ( *, '("PrimSType:  ",I4)') cModelParams.PrimarySourceType
    write ( *, '("ContSType:  ",I4)') cModelParams.ContrastSourceType
    write ( *, '("FDTorder:   ",I4)') cModelParams.FDTorder
    write ( *, '("FDXorder:   ",I4)') cModelParams.FDXorder
    write ( *, '("Debuglvl:   ",I4)') iDebugLvl
    write ( *, '("A-Aliasing: ",L4)') cModelParams.UseAntiAliasing
    write ( *, '("SupTaper:   ",L4)') cModelParams.UseSupportTapering
    write ( *, '("FreqTaper:  ",L4)') cModelParams.UseFreqTapering
    write ( *, '("Numbeams:   ",I4)') cModelParams.Numbeams
    do i=1,cModelParams%Numbeams
        write ( *, '(" NumIt(",I2,"):   ",I4)') i,cModelParams.numiterations(i)
    end do
    write ( *, '("Slicesavespecifier: ",I5)') cModelParams.slicesavespecifier
    write ( *, '("Numslices:          ",I5)') cModelParams.Numslices
    do i=1,cModelParams%Numslices
        write ( *, '(" xyzslicedim(",I2,"):    ",A1)') i,cModelParams.xyzslicedim(i)
        write ( *, '(" xyzslicepos(",I2,"):    ",E10.3)') i,cModelParams.xyzslicepos(i)
    end do
    write ( *, '(" ")')
    write ( *, '("Srcfilename: ",A)') cSourceParams.srcfilename
    if(trim(cSourceParams.srcfilename)=='none') then
        ! source signature parameters
        write ( *, '(" SrcSignType: ",I2)') cSourceParams.srcsigntype
        if (cSourceParams.srcsigntype==iGI_GAUSSIAN) then
            write ( *, '(" Pstart: ",E10.3)') cSourceParams.Pstart
            write ( *, '(" Tpulse: ",E10.3)') cSourceParams.Tpulse
            write ( *, '(" Twidth: ",E10.3)') cSourceParams.Twidth
            write ( *, '(" Tdelay: ",E10.3)') cSourceParams.Tdelay
            write ( *, '(" Power:  ",I3)') cSourceParams.Power
        else if (cSourceParams.srcsigntype==iGI_BLACKMAN) then
            write ( *, '(" Pstart: ",E10.3)') cSourceParams.Pstart
            write ( *, '(" Tpulse: ",E10.3)') cSourceParams.Tpulse
            write ( *, '(" Tdelay: ",E10.3)') cSourceParams.Tdelay
        else if (cSourceParams.srcsigntype==iGI_FILE) then
            write (*, "('Source signature file ',A)") trim(cSourceParams.srcsignfilename)
        end if
        write ( *, '(" ")')
        if (cModelParams.PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams.PrimarySourceType/=iEI_POINTSOURCECLOUD) then
            ! source shape parameters
            write ( *, '(" SrcShapeType: ",I2)') cSourceParams.srcshapetype
            if (cSourceParams.srcshapetype == iHI_CYLINDRICAL) then
                write ( *, '(" Radius:   ",E10.3)') cSourceParams.radius
                write ( *, '(" Radfocus: ",E10.3)') cSourceParams.radfocus
            else if (cSourceParams.srcshapetype == iHI_RECTANGULAR) then
                write ( *, '(" Rectwidth:  ",E10.3)') cSourceParams.rectwidth
                write ( *, '(" Rectheight: ",E10.3)') cSourceParams.rectheight
            else if (cSourceParams.srcshapetype == iHI_TRIANGULAR) then
                write ( *, '(" Triangedge:    ",E10.3)') cSourceParams.triangedge
                write ( *, '(" Triangxcenter: ",E10.3)') cSourceParams.triangxcenter
            else if (cSourceParams.srcshapetype == iHI_POINTSOURCE) then
            else if (cSourceParams.srcshapetype == iHI_PHASEDARRAY) then
                write ( *, '(" Numel:        ",I3)') cSourceParams.numel
                write ( *, '(" Elwidth:      ",E10.3)') cSourceParams.elwidth
                write ( *, '(" Elheight:     ",E10.3)') cSourceParams.elheight
                write ( *, '(" Kerf:         ",E10.3)') cSourceParams.kerf
                write ( *, '(" Focusx:       ",E10.3)') cSourceParams.focusx
                write ( *, '(" Focusz:       ",E10.3)') cSourceParams.focusz
                write ( *, '(" Elfocus:      ",E10.3)') cSourceParams.elevationfocusz
                write ( *, '(" TD file:       ",A)') cSourceParams.td_filename
                write ( *, '(" Apodfile:      ",A)') cSourceParams.phaseapodfilename
                !KH		else if (cSourceParams.srcshapetype == iHI_MATRIXARRAY) then
                !KH				write ( *, '(" Numelx:    ",I3)') cSourceParams.matnumelx
                !KH				write ( *, '(" Numely:    ",I3)') cSourceParams.matnumely
                !KH				write ( *, '(" Elwidth:  ",E10.3)') cSourceParams.matelwidth
                !KH				write ( *, '(" Elheight: ",E10.3)') cSourceParams.matelheight
                !KH				write ( *, '(" Kerfx:     ",E10.3)') cSourceParams.matkerfx
                !KH				write ( *, '(" Kerfy:     ",E10.3)') cSourceParams.matkerfy
                !KH				write ( *, '(" Focusx:   ",E10.3)') cSourceParams.matfocusx
                !KH				write ( *, '(" Focusy:   ",E10.3)') cSourceParams.matfocusx
                !KH				write ( *, '(" Focusz:   ",E10.3)') cSourceParams.matfocusz
                !KH!				write ( *, '(" Apodfile: ",A)') cSourceParams.phaseapodfilename
                write ( *, '(" ")')
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '(" ")')
                write ( *, '(" Symmetry in the Y axis:             ",I3)') abs(cModelParams%UseYSymmetry)
                write ( *, '(" Nonlinear Compressibility Contrast: ",I3)') cSourceParams.NKC
                write ( *, '(" Linearized Contrast Source Method:  ",I3)') cSourceParams.LCSM
                write ( *, '(" Bi-CGSTAB:                          ",I3)') cSourceParams.BiCGSTAB
                write ( *, '(" Steepest Descent:                   ",I3)') cSourceParams.CG
                write ( *, '(" Linear case Steepest Descent:       ",I3)') cSourceParams.LINEARCASEIMPLEMENTATION
                write ( *, '(" Dual Frequency Modality:            ",I3)') cSourceParams.DFM
                write ( *, '(" Shift for Dual Frequency Modality:  ",I3)') cSourceParams.SHIFT
                write ( *, '(" ")')
                write ( *, '("---------------------------------------------------------------------------")')

                write ( *, '(" ")')
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '("Parameters Section ")')
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '(" ")')
                write ( *, '(" alpha1: ",E10.3)') cSourceParams.alpha1
                write ( *, '(" b1:     ",E10.3)') cSourceParams.b1
                write ( *, '(" SOS1:   ",E10.3)') cSourceParams.SOS1
                write ( *, '(" KAPPA1: ",E10.3)') cSourceParams.KAPPA1
                write ( *, '(" BETA1:  ",E10.3)') cSourceParams.BETA1
                write ( *, '(" ")')
                write ( *, '(" alpha2: ",E10.3)') cSourceParams.alpha2
                write ( *, '(" b2:     ",E10.3)') cSourceParams.b2
                write ( *, '(" SOS2:   ",E10.3)') cSourceParams.SOS2
                write ( *, '(" KAPPA2: ",E10.3)') cSourceParams.KAPPA2
                write ( *, '(" BETA2:  ",E10.3)') cSourceParams.BETA2
                write ( *, '(" ")')
                write ( *, '(" alpha3: ",E10.3)') cSourceParams.alpha3
                write ( *, '(" b3:     ",E10.3)') cSourceParams.b3
                write ( *, '(" SOS3:   ",E10.3)') cSourceParams.SOS3
                write ( *, '(" KAPPA3: ",E10.3)') cSourceParams.KAPPA3
                write ( *, '(" BETA3:  ",E10.3)') cSourceParams.BETA3
                write ( *, '(" ")')
                write ( *, '(" alpha4: ",E10.3)') cSourceParams.alpha4
                write ( *, '(" b4:     ",E10.3)') cSourceParams.b4
                write ( *, '(" SOS4:   ",E10.3)') cSourceParams.SOS4
                write ( *, '(" KAPPA4: ",E10.3)') cSourceParams.KAPPA4
                write ( *, '(" BETA4:  ",E10.3)') cSourceParams.BETA4
                write ( *, '(" ")')
                write ( *, '(" alpha5: ",E10.3)') cSourceParams.alpha5
                write ( *, '(" b5:     ",E10.3)') cSourceParams.b5
                write ( *, '(" SOS5:   ",E10.3)') cSourceParams.SOS5
                write ( *, '(" KAPPA5: ",E10.3)') cSourceParams.KAPPA5
                write ( *, '(" BETA5:  ",E10.3)') cSourceParams.BETA5
                write ( *, '(" ")')
                write ( *, '("--------------------------------------------------------------------------")')
                write ( *, '("Residual :",A)') trim(cSourceParams.residualdirectory)
                write ( *, '("--------------------------------------------------------------------------")')
            
            
            end if
            
        end if
        if (cModelParams.ContrastSourceType==iCI_SCATTERER) then
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '("Parameters for scatterer")')
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '(" Type: 	     ",A)'	  		) trim(ScattererParams%ScType)
				if (  trim(ScattererParams%ScType) == 'microbubble') then 
					write ( *, '(" kappa_s:   ",E10.3)') ScattererParams%kappa_s
					write ( *, '(" sigma_w:   ",E10.3)') ScattererParams%sigma_w
					write ( *, '(" sigma_R0:  ",E10.3)') ScattererParams%sigma_R0
					write ( *, '(" gamma:     ",E10.3)') ScattererParams%gama
					write ( *, '(" chi:       ",E10.3)') ScattererParams%chi
					write ( *, '(" Marmottant Solver:       ",A)') trim(ScattererParams%Solver_Method)
					write ( *, '(" Normalize in:       ",A)') trim(ScattererParams%Solver_Normalize)
					write ( *, '(" ")')
				elseif (  trim(ScattererParams%ScType) == 'linear') then 
					write ( *, '(" rho1:   ",E10.3)') ScattererParams%rho1
					write ( *, '(" c1:   ",E10.3)') ScattererParams%c1
				endif
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '("Cluster Parameters")')
                write ( *, '("---------------------------------------------------------------------------")')
                write ( *, '(" ")')
                write ( *, '(" Distribution: 	     ",A)'	  		) trim(ScattererParams%Distribution)
            if (trim(ScattererParams%Distribution) == 'monodisperse') then
	    		write ( *, '(" Radius:		 	     ",E10.3)'	  	) ScattererParams%R0(1)
            elseif (trim(ScattererParams%Distribution) == 'polydisperse') then
	    		write ( *, '(" Radius Range:	     ",E10.3,E10.3)') MINVAL(ScattererParams%R0),MAXVAL(ScattererParams%R0)
                endif
                write ( *, '(" No. of Bubbles:       ",I<int(log10(real(ScattererParams%N,dp))+1)>)'  ) ScattererParams%N
                write ( *, '(" [X_Bmin X_Bmax]/ Lx:  ",E10.3,E10.3)') ScattererParams%ClusterDimsRatio(:,1)
                write ( *, '(" [Y_Bmin Y_Bmax]/ Ly:  ",E10.3,E10.3)') ScattererParams%ClusterDimsRatio(:,2)
                write ( *, '(" [Z_Bmin Z_Bmax]/ Lz:  ",E10.3,E10.3)') ScattererParams%ClusterDimsRatio(:,3)
                write ( *, '(" Minimum Bubble Dist:  ",E10.3)') ScattererParams%MinInBetweenDist
                write ( *, '(" Calculate Pressure to ",A," gridpoints")') trim(ScattererParams%GridPointsPressure)
        endif
    end if
    write ( *, '(" ")')
     write ( *, '("---------------------------------------------------------------------------")')

    END SUBROUTINE ConfigFileInfo

SUBROUTINE SpecialChecksConfigParameters

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine SpecialChecksConfigParameters sets and checks the following
    !   regarding the configuration parameters
    !			- Checks a number of invalid parameter values
    !			- Are the FD orders even? if not, correct
    !			- Check and build input/output files of accompanying files
    !			- Is there symmetry in the y-dimension? set UseYSymmetry flag and
    !               divide ly by 2.
    !			- If the first position is sz = 0, reset it to \delta z, and reset
    !               lz as well.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   none
    !
    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   i        i8b    temporary variable
    !   err      i8b    error number
    !   error    lgt    error flag
    !   acTemp   char   Temporary char array for output log messages
    !
    integer(i8b) ::                 i,err
    logical(lgt) ::					error
    character(LEN = 1024)::         acTemp;

    ! *****************************************************************************
    !
    !   I/O
    !
    !   Echoes the error messages to screen, checks the existence of input files
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

#ifndef RUNESTIMATOR
    if(cModelParams.primarysourcetype == 0) then
        write(*,"('Error, primary source type is not recognized')")
        stop
    end if

    if(cModelParams.contrastsourcetype == 0) then
        write(*,"('Error, contrast source type is not recognized')")
        stop
    end if

    if(cModelParams.slicesavespecifier == 0) then
        write(*,"('Error, slice save specifier is not recognized')")
        stop
    end if

    if(trim(cSourceParams.srcfilename)=='none') then
        if(cSourceParams.srcsigntype == 0) then
            write(*,"('Error, source signature type is not recognized')")
            stop
        end if
        if (cModelParams.PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams.PrimarySourceType/=iEI_POINTSOURCECLOUD) then
            if(cSourceParams.srcshapetype == 0) then
                write(*,"('Error, source shape type is not recognized')")
                stop
            end if
        end if
    end if

    if (cModelParams.Numbeams<=0) then
        write(*,"('Error, number of beams should be larger than zero')")
        stop
    end if

    if (any(cModelParams.NumIterations==0)) then
        write(*,"('Error, iterations for all beams should be larger than zero')")
        stop
    end if

    !Check Numiterations
    if (cModelParams%Numbeams>1) then
        do i=2,cModelparams%Numbeams
            if (cModelParams%Numiterations(i)>cModelParams%Numiterations(i-1)) then
                write(*,"('Error: Numiterations should not increase for increasing beams')")
                stop
            end if
        end do
    end if

    if (cModelParams.Numslices<=0) then
        write(*,"('Error, number of slices should be larger than zero')")
        stop
    end if

    if ((cModelParams.PrimarySourceType==iEI_PLANEWAVE).and.(trim(cSourceParams.srcfilename)/='none')) then
        write(*,"('Error, plane wave solution cannot be combined with a source inputfile')")
        stop
    end if
#endif

    !srcfilename
    if (trim(cSourceParams.srcfilename)=="standard") then
        i=index(sInputName,'.')
        cSourceParams.srcfilename=sInputName(1:i-1)//'_src.bin'
    end if

    if (trim(cSourceParams%srcfilename)/="none") then
#ifndef RUNESTIMATOR
        inquire(FILE=trim(sInputDir)//cSourceParams.srcfilename,EXIST=error)
        if(error .EQV. .FALSE.) then
            write(*,"('Error, source inputfile ',A,' cannot be found.')") trim(sInputDir)//trim(cSourceParams%srcfilename)
            stop
        end if
#endif
    else
        if (cModelParams.PrimarySourceType==iEI_PLANEWAVE) then
            ! If the linear field solution is the plane wave field, then
            ! treat source configuration as though it is a point source
            ! further definition of the source shape type is ignored
            cSourceParams.srcshapetype= iHI_POINTSOURCE
        end if
        !srcsigntype
        if (cSourceParams%srcsigntype == iGI_FILE) then
            if (trim(cSourceParams.srcsignfilename)=="standard") then
                i=index(sInputName,'.')
                cSourceParams.srcsignfilename=sInputName(1:i-1)//'_sig.bin'
            end if
#ifndef RUNESTIMATOR
            inquire(FILE=trim(sInputDir)//cSourceParams.srcsignfilename,EXIST=error)
            if(error .EQV. .FALSE.) then
                write(*,"('Error, signature inputfile ',A,' cannot be found.')") trim(sInputDir)//trim(cSourceParams%srcsignfilename)
                stop
            end if
#endif
        end if

		!td_filename
        if (cSourceParams%srcshapetype == iHI_PHASEDARRAY) then
            if (trim(cSourceParams.td_filename)=="standard") then
                i=index(sInputName,'.')
                cSourceParams.td_filename=sInputName(1:i-1)//'_td.bin'
            end if
            if     (trim(cSourceParams.td_filename) == 'none') then
            else
#ifndef RUNESTIMATOR
                inquire(FILE=trim(sInputDir)//trim(cSourceParams.td_filename),EXIST=error)
                if(error .EQV. .FALSE.) then
                    write(*,"('Error, time delay inputfile ',A,' cannot be found.')") trim(cSourceParams.td_filename)
                    stop
                end if
#endif
            end if
			
            if (trim(cSourceParams.phaseapodfilename)=="standard") then
                i=index(sInputName,'.')
                cSourceParams.phaseapodfilename=sInputName(1:i-1)//'_phapo.bin'
            end if
            if     (trim(cSourceParams.phaseapodfilename) == 'none') then
            else
#ifndef RUNESTIMATOR
                inquire(FILE=trim(sInputDir)//trim(cSourceParams.phaseapodfilename),EXIST=error)
                if(error .EQV. .FALSE.) then
                    write(*,"('Error, apodization inputfile ',A,' cannot be found.')") trim(cSourceParams.phaseapodfilename)
                    stop
                end if
#endif
            end if
        end if

    end if

    !Check FD orders: should be even
    if (mod(int(cModelParams.FDTorder),2)==1) then
#ifndef RUNESTIMATOR
        write(*,"('Warning: FDTorder should be even, correcting it to one order higher')")
#endif
        cModelParams.FDTorder = cModelParams.FDTorder+1
    end if
    if (mod(int(cModelParams.FDXorder),2)==1) then
#ifndef RUNESTIMATOR
        write(*,"('Warning: FDXorder should be even, correcting it to one order higher')")
#endif
        cModelParams.FDXorder = cModelParams.FDXorder+1
    end if

    !Use Y symmetry: only if Thetay=pi/2 and and we use nonlinear contrast
    if ((abs(cModelParams.Thetay-half_pi)<1e-5) .AND. cModelParams.ContrastSourceType==iCI_NONLIN .AND. cModelParams.UseYSymmetry) then
        cModelParams.UseYSymmetry=.true.
        cModelParams.Ly=cModelParams.Ly/2.0_dp
    else
        cModelParams.UseYSymmetry=.false.
    end if

    !If Sz was zero, set to dz - one stepsize away from the source - in mm
    !Also reset Lz
    if(abs(cModelParams.Sz)<1e-6) then
        cModelParams.Sz = cMediumParams.c0 / (2.0_dp*cModelParams.freq0*cModelParams.FNyq) * 1.0d3
        cModelParams.Lz = cModelParams.Lz - cModelParams.Sz
    end if

    END SUBROUTINE SpecialChecksConfigParameters

    SUBROUTINE NormalizeConfigParameters

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NormalizeConfigParameters normalize the constants to the
    !   normalization as applied within the program: all temporal dimensions with
    !   respect to the fundamental period 1/freq0, all spatial dimensions with
    !   respect to the fundamental wavelength c0/freq0. Some derived medium
    !   parameters are also set.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   none
    !
    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   dLambdaMM         dp   Fundamental wavelength in mm
    !   dMaxPhaseDelayX   dp   maximum phase delay of the phased array in the x-
    !                          dimension, in fundamental periods.
    !   dMaxPhaseDelayY   dp   maximum phase delay of the phased array in the y-
    !                          dimension, in fundamental periods.
    !   eltdelay          dp   delay of the elements of the phased array
    !   elx               dp   x-position of the elements of the phased array
    !   i                i8b   temporary variable
    !
    real(dp) :: dLambdaMM
    real(dp) :: dMaxPhaseDelayX,dMaxPhaseDelayY
    real(dp),allocatable :: eltdelay(:),elx(:),sourcetdtemp(:)
    integer(i8b) :: i

    ! *****************************************************************************
    !
    !   I/O
    !
    !   none
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    !Normalization for all spatial parameters
    dLambdaMM = cMediumParams.c0/cModelParams.freq0 * 1.0d3;

    !extra medium parameters
    cMediumParams.kappa0     = 1.0_dp/(cMediumParams.rho0 * cMediumParams.c0**2)
    cMediumParams.Z0         = cMediumParams.rho0*cMediumParams.c0
    cMediumParams.alpha0_att = cMediumParams.a_att*100.0_dp/(TWO_PI*1.0E6_dp)**cMediumParams.b_att

    !normalize domain parameters
    cModelParams.Lx = cModelParams.Lx / dLambdaMM;
    cModelParams.Ly = cModelParams.Ly / dLambdaMM;
    cModelParams.Lz = cModelParams.Lz / dLambdaMM;
    cModelParams.Sz = cModelParams.Sz / dLambdaMM;
    cModelParams.Lt = cModelParams.Lt
    !Lt is already in time periods...

    !normalize slices
    do i=1,cModelParams%Numslices
        cModelParams%xyzslicepos(i) = cModelParams%xyzslicepos(i) / dLambdaMM;
    end do

    !normalize source parameters

    if (trim(cSourceParams%srcfilename) == 'none') then
        if (cSourceParams.srcshapetype==iHI_CYLINDRICAL) then
            cSourceParams.radius=cSourceParams.radius / dLambdaMM;
            cSourceParams.radfocus=cSourceParams.radfocus / dLambdaMM;
            if(abs(cSourceParams.radfocus)>1e-5) then
                cSourceParams.maxphasedelay=sqrt(cSourceParams.radius**2+cSourceParams.radfocus**2)-abs(cSourceParams.radfocus)
            else
                cSourceParams.maxphasedelay=0
            end if
        else if (cSourceParams.srcshapetype==iHI_RECTANGULAR) then
            cSourceParams.rectwidth=cSourceParams.rectwidth / dLambdaMM;
            cSourceParams.rectheight=cSourceParams.rectheight / dLambdaMM;
            cSourceParams.maxphasedelay=0
        else if (cSourceParams.srcshapetype==iHI_TRIANGULAR) then
            cSourceParams.triangedge=cSourceParams.triangedge / dLambdaMM;
            cSourceParams.triangxcenter=cSourceParams.triangxcenter / dLambdaMM;
            cSourceParams.maxphasedelay=0
        else if (cSourceParams.srcshapetype==iHI_POINTSOURCE) then
            cSourceParams.maxphasedelay=0
        else if (cSourceParams.srcshapetype==iHI_PHASEDARRAY) then
            cSourceParams.elwidth=cSourceParams.elwidth / dLambdaMM;
            cSourceParams.elheight=cSourceParams.elheight / dLambdaMM;
            cSourceParams.kerf=cSourceParams.kerf / dLambdaMM;
            cSourceParams.focusx=cSourceParams.focusx / dLambdaMM;
            cSourceParams.focusz=cSourceParams.focusz / dLambdaMM;
            cSourceParams.elevationfocusz=cSourceParams.elevationfocusz / dLambdaMM;

            !derived parameters arraywidth and maxphasedelay
            cSourceParams.arraywidth=real(cSourceParams.numel,dp)*(cSourceParams.elwidth+cSourceParams.kerf)-cSourceParams.kerf

            !max phase delay in the x-direction due to element focusing
            allocate(eltdelay(1:cSourceParams%numel),						&
                elx(1:cSourceParams%numel))

            elx=-cSourceParams%arraywidth/2.0_dp+cSourceParams%elwidth/2.0_dp+ &
                (/ (i,i=0,cSourceParams%numel-1) /)*(cSourceParams%elwidth+cSourceParams%Kerf)
			if (trim(cSourceParams%td_filename)=='none') then
				eltdelay=-sqrt((cSourceParams%focusx-elx)**2+cSourceParams%focusz**2)*SIGN(1.0D0,cSourceParams%focusz)
			else
				i=index(cSourceParams%td_filename,'.')
				allocate(sourcetdtemp(1:cSourceParams%numel))
				if (cSourceParams%td_filename(i+1:i+3)=='bin') then
					if (BINARY_read1Darraysize(trim(sInputDir)//trim(cSourceParams%td_filename))/=cSourceParams%numel) then
						call PrintToLog("Error: number of items in apodization file is not equal to number of array elements",1)
						stop
					end if
					call BINARY_read1Darray(trim(sInputDir)//trim(cSourceParams%td_filename),sourcetdtemp)
				else
					if (ASCII_read1Darraysize(trim(sInputDir)//trim(cSourceParams%td_filename))/=cSourceParams%numel) then
						call PrintToLog("Error: number of items in apodization file is not equal to number of array elements",1)
						stop
					end if
					call ASCII_read1Darray(trim(sInputDir)//trim(cSourceParams%td_filename),sourcetdtemp)
				end if
				eltdelay=sourcetdtemp*cModelParams.freq0*1E-6; ! Normalize microsec 
				deallocate(sourcetdtemp)
			end if
			dMaxPhaseDelayX=maxval(eltdelay-minval(eltdelay))

            deallocate(eltdelay,elx)

!			!alternative approach to determine dMaxPhaseDelayX
!			dMaxPhaseDelayX = sqrt((abs(cSourceParams.focusx)+cSourceParams.arraywidth/2.0)**2 &
!									+cSourceParams.focusz**2)-abs(cSourceParams.focusz);
!			if (abs(cSourceParams.focusx)>cSourceParams.arraywidth/2.0) then
!				!subtract the delay of the element closest to focusx
!				! if |focusx|>arrwidth/2
!				dMaxPhaseDelayX = dMaxPhaseDelayX - &
!							(sqrt((abs(cSourceParams.focusx)-cSourceParams.arraywidth/2.0)**2 &
!									+cSourceParams.focusz**2)-abs(cSourceParams.focusz));
!			end if

            !phase delay in the y-direction due to elevation focusing
            if(cSourceParams.elevationfocusz>1.0e-5_dp) then
                dMaxPhaseDelayY = sqrt((cSourceParams.elheight/2.0)**2 &
                    +cSourceParams.elevationfocusz**2) - abs(cSourceParams.elevationfocusz);
            else
                dMaxPhaseDelayY = 0
            end if

            cSourceParams.maxphasedelay=dMaxPhaseDelayX + dMaxPhaseDelayY

            !KH			else if (cSourceParams.srcshapetype==iHI_MATRIXARRAY) then

            !KH         !normalize parameters and determine maximum phase delay

        end if
    end if

    END SUBROUTINE NormalizeConfigParameters

    SUBROUTINE BINARY_read1Darray(filename,darr)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine BINARY_read1Darray reads a fortran binary file holding a
    !   1D array of real*8 data. The binary file is an unformatted file that
    !   contains the following variables:
    !
    ! 		integer*8    arrdim(1:2)              size of the dimensions of the
    !                                             array, either (1,Numel) or
    !                                             (Numel,1)
    ! 		real*8    darr(arrdim(1),arrdim(2))   the array
    !
    !   The standard fortran format to store data in an unformatted file is for
    !   each variable to store an integer*4 number containing the number of
    !   bytes used by the variable before and after the variable.
    !
    !   The size of the array that is input to the subroutine needs to be equal to
    !   the total size of the  array in the file. Read it with the function
    !   BINARY_read1Darraysize and allocate the array before calling
    !   BINARY_read1Darray.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename   i   char   name of the file
    !   darr       o    dp    data array as read from the file
    !
    character(*), intent(in) :: filename
    real(dp), intent(out) :: darr(:)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err   i8b   error number
    !   da    i8b   2-element vector containing the size of the array
    !
    integer(i8b) :: err, da(2)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read a binary file containing a 1D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=trim(filename),status="OLD",form="UNFORMATTED")

    read(iImportUNIT) da

    read(iImportUNIT) darr

    close(iImportUNIT)

    END SUBROUTINE BINARY_read1Darray

    FUNCTION BINARY_read1Darraysize(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The function BINARY_read1Darraysize reads the number of elements stored
    !   in a binary file with a 1D array. See the header of BINARY_read1Darray for
    !   more information.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                i   char   name of the file
    !   BINARY_read1Darraysize  o   i8b    size of the array
    !
    character(*), intent(in) :: filename
    integer(i8b) :: BINARY_read1Darraysize

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err   i8b   error number
    !   da    i8b   2-element vector containing the size of the array
    !
    integer(i8b) :: err, da(2)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read number of elements from a binary file containing a 1D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=trim(filename),status="OLD",form="UNFORMATTED")

    read(iImportUNIT) da

    close(iImportUNIT)

    BINARY_read1Darraysize = da(1)*da(2)

    END FUNCTION BINARY_read1Darraysize

    SUBROUTINE BINARY_read_sourcedef(filename,darr,dTdelay)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine BINARY_read_sourcedef reads a fortran binary file holding a
    !   source definition in a 3D array of real*8 data. The file also includes the
    !   space-dependent time delays. This subroutine is NOT a generic reading of a
    !   3D array from a binary file!!
    !   The binary file is an unformatted file that holds the following variables:
    !
    ! 		integer*8  arrdim(1:3)            size of the dimensions of the array
    !       integer*8  iTdelay                number of points to take extra in the
    !                                         temporal dimension of the source to
    !                                         be able to hold the delayed source
    !                                         pulse at all positions
    !       real*8     dTdelay(arrdim(2),arrdim(3))  time delays as function of
    !                                                x and y, in milliseconds
    ! 		real*8     darr(arrdim(1),arrdim(2),arrdim(3))   the array
    !
    !   The standard fortran format to store data in an unformatted file is for
    !   each variable to store an integer*4 number containing the number of
    !   bytes used by the variable before and after the variable.
    !
    !   The size of the array that is input to the subroutine needs to be equal to
    !   the total size of the  array in the file. Read it with the function
    !   BINARY_read3Darraysize and allocate the array before calling
    !   BINARY_read_sourcedef.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename   i   char   name of the file
    !   darr       o    dp    data array as read from the file
    !   dTdelay    o    dp    space-dependent time delays, in milliseconds
    !
    character(*), intent(in) :: filename
    real(dp), intent(out) :: darr(:,:,:)
    real(dp), intent(out) :: dTdelay(:,:)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   err         i8b   error number
    !   da          i8b   3-element vector containing the size of the array
    !   dMaxTdelay  dp    maximum time delay in milliseconds
    !
    integer(i8b) :: err, da(3)
    real(dp)   :: dMaxTdelay

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read a binary file containing a primary source definition in the form of
    !   a 3D array and a time delay array.
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=filename,status="OLD",form="UNFORMATTED",iostat=err)
    if (err/=0) then
        write(*,"('Error reading file ',A,', error code ',I3)") trim(filename),err
        stop
    end if

    read(iImportUNIT) da
    read(iImportUNIT) dMaxTdelay
    read(iImportUNIT) dTdelay
    read(iImportUNIT) darr

    close(unit=iImportUNIT)

    END SUBROUTINE BINARY_read_sourcedef

    FUNCTION BINARY_read3Darraysize(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The function BINARY_read3Darraysize reads the dimensions of the array
    !   stored in a binary file with a 3D array for the primary source definition.
    !   See the header of BINARY_read_sourcedef for more information.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                i   char   name of the file
    !   BINARY_read3Darraysize  o   i8b    dimensions of the array
    !
    character(*), intent(in) :: filename
    integer(i8b) :: BINARY_read3Darraysize(3)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err   i8b   error number
    !   da    i8b   3-element vector containing the size of the array
    !
    integer(i8b) :: err, da(3)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read the array dimensions from a binary file containing a 3D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=trim(filename),status="OLD",form="UNFORMATTED")

    read(iImportUNIT) da

    close(iImportUNIT)

    BINARY_read3Darraysize=da

    END FUNCTION BINARY_read3Darraysize

    FUNCTION BINARY_read3DarrayTdelay(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The function BINARY_read3DarrayTdelay reads the maximum time delay (in
    !   milliseconds) used in a binary source definition file. See the header of
    !   BINARY_read_sourcedef for more information.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                  i   char  name of the file
    !   BINARY_read3DarrayTdelay  o   dp    maximum time delay in milliseconds
    !
    character(*), intent(in) :: filename
    real(dp) :: BINARY_read3DarrayTdelay

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err     i8b   error number
    !   da      i8b   3-element vector containing the size of the array
    !   dTdelay dp    maximum time delay in milliseconds
    !
    integer(i8b) :: err, da(3)
    real(dp) :: dTdelay

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read the maximum time delay from a binary file containing a 3D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=trim(filename),status="OLD",form="UNFORMATTED")

    read(iImportUNIT) da
    read(iImportUNIT) dTdelay

    close(iImportUNIT)

    BINARY_read3DarrayTdelay = dTdelay

    END FUNCTION BINARY_read3DarrayTdelay

    SUBROUTINE BINARY_save3Darray_test(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine BINARY_save3Darray_test writes a test binary file to check
    !   the format in which Fortran writes its unformatted files.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                  i   char  name of the file
    !
    character(*), intent(in) :: filename

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err     i8b   error number
    !   da      i8b   3-element i8b test vector
    !   dar     dp    3-element dp test vector
    !
    integer(i8b) :: err, da(3)
    real(dp) :: dar(3)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   writes a test binary file
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    da = (/ 3, 4, 5 /)
    dar = real(da,dp)

    open(unit=iExportUNIT,file=trim(filename),status="NEW",form="UNFORMATTED")

    write(iExportUNIT) da

    write(iExportUNIT) dar

    close(iExportUNIT)

    END SUBROUTINE BINARY_save3Darray_test

    SUBROUTINE ASCII_read1Darray(filename,darr)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine BINARY_read1Darray reads an ASCII file holding a
    !   1D array of double precision data. The correct format is obtained when
    !   the file is created with the following Matlab command line:
    !
    !   sizearr=size(array); save <filename> -ASCII -DOUBLE sizearr array
    !
    !   It holds the following variables:
    !
    ! 		integer*8    arrdim(1:2)              size of the dimensions of the
    !                                             array, either (1,Numel) or
    !                                             (Numel,1)
    ! 		real*8    darr(arrdim(1)*arrdim(2))   the array, stored either
    !                                             columnwise or rowwise
    !
    !   The size of the array that is input to the subroutine needs to be equal to
    !   the total size of the  array in the file. Read it with the function
    !   ASCII_read1Darraysize and allocate the array before calling
    !   ASCII_read1Darray.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename   i   char   name of the file
    !   darr       o    dp    data array as read from the file
    !
    character(*), intent(in) :: filename
    real(dp), intent(out) :: darr(:)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   ia       i8b   temporary variable for determining the array format string
    !                   for reading one row
    !   ib       i8b   idem
    !   ic       i8b   idem
    !   ipower   i8b   maximum order of 10 in the number of columns
    !   err      i8b   error number
    !   arrnumrows   i8b   number of rows in the array
    !   arrnumcols   i8b   number of columns in the array
    !   da       dp    2-element vector containing the size of the array
    !   cformat  char  format string for reading one row
    !
    integer(i8b) :: ia,ib,ic, ipower, err, arrnumrows,arrnumcols
    real(dp) :: da(2)
    character(1024) :: cformat

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read an ASCII file containing a 1D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=filename,status="OLD",form="FORMATTED",iostat=err)
    if (err/=0) then
        write(*,"('Error reading file ',A,', error code ',I3)") trim(filename),err
        stop
    end if

    !First two numbers in the file are the sizes of the array: numrows, numcols
    read(unit=iImportUNIT,iostat=err,fmt="(2ES25.16E3)") da

    arrnumrows = da(1)
    arrnumcols = da(2)

    !Build format string
    cformat="("
    ic=arrnumcols
    ipower=log10(real(ic))
    do ia = ipower,1,-1
        ib = ic/10**ia
        cformat=trim(cformat)//char(48+ib)
        ic = ic-ib*10**ia
    end do
    cformat=trim(cformat)//char(48+ic)//"ES25.16E3)"

    !Read array
    do ia = 0,arrnumrows-1
        read(unit=iImportUNIT,fmt=cformat) darr(ia*arrnumcols+1:ia*arrnumcols+arrnumcols)
    end do

    close(unit=iImportUNIT)

    END SUBROUTINE ASCII_read1Darray

    FUNCTION ASCII_read1Darraysize(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The function ASCII_read1Darraysize reads the number of elements stored
    !   in an ASCII file with a 1D array. See the header of ASCII_read1Darray for
    !   more information.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                i   char   name of the file
    !   ASCII_read1Darraysize   o   i8b    size of the array
    !
    character(*), intent(in) :: filename
    integer(i8b) :: ASCII_read1Darraysize

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err   i8b   error number
    !   da    i8b   2-element vector containing the size of the array
    !
    integer(i8b) :: err
    real(dp) :: da(2)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read number of elements from an ASCII file containing a 1D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=filename,status="OLD",form="FORMATTED",iostat=err)
    if (err/=0) then
        write(*,"('Error reading temporary file ',A,', error code ',I3)") trim(filename),err
        stop
    end if

    read(unit=iImportUNIT,iostat=err,fmt="(2ES25.16E3)") da

    ASCII_read1Darraysize = da(1)*da(2)

    close(unit=iImportUNIT)

    END FUNCTION ASCII_read1Darraysize

    SUBROUTINE ASCII_read3Darray(filename,darr)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine BINARY_read3Darray reads an ASCII file holding a
    !   3D array of double precision data. The correct format is obtained when
    !   the file is created with the following Matlab command line:
    !
    !   sizearr=size(array); array=array(:); save <filename> -ASCII -DOUBLE sizearr array
    !
    !   Mark the flattening of the 3D array before it is saved.
    !
    !   The file holds the following variables:
    !
    ! 		integer*8    arrdim(1:2)              size of the dimensions of the
    !                                             array, either (1,Numel) or
    !                                             (Numel,1)
    ! 		real*8    darr(arrdim(1)*arrdim(2)*arrdim(3))   the array
    !
    !   The size of the array that is input to the subroutine needs to be equal to
    !   the total size of the  array in the file. Read it with the function
    !   ASCII_read3Darraysize and allocate the array before calling
    !   ASCII_read3Darray.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename   i   char   name of the file
    !   darr       o    dp    data array as read from the file
    !
    character(*), intent(in) :: filename
    real(dp), intent(out) :: darr(:,:,:)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   ia         i8b   temporary variable for determining the array format string
    !                    for reading one row
    !   ib         i8b   idem
    !   ic         i8b   idem
    !   arrdim     i8b   dimensions of the array
    !   ipower     i8b   maximum order of 10 in the number of columns
    !   err        i8b   error number
    !   da         dp    3-element vector containing the size of the array
    !   cformat    char  format string for reading one row
    !
    integer(i8b) :: ia,ib,ic, err, arrdim(3), ipower
    real(dp) :: da(3)
    character(1024) :: cformat

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read an ASCII file containing a 3D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=filename,status="OLD",form="FORMATTED",iostat=err)
    if (err/=0) then
        write(*,"('Error reading file ',A,', error code ',I3)") trim(filename),err
        stop
    end if

    !First three numbers in the file are the sizes of the array
    read(unit=iImportUNIT,iostat=err,fmt="(3ES25.16E3)") da

    arrdim = da

    !Build format string
    cformat="(1ES25.16E3)"

    !Read array, in column-major order (used in Matlab as well as in Fortran)
    do ic = 0,arrdim(3)-1
        do ib = 0,arrdim(2)-1
            do ia = 0,arrdim(1)-1
                read(unit=iImportUNIT,fmt=cformat) darr(ia+1,ib+1,ic+1)
            end do
        end do
    end do

    close(unit=iImportUNIT)

    END SUBROUTINE ASCII_read3Darray

    FUNCTION ASCII_read3Darraysize(filename)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The function ASCII_read3Darraysize reads the dimensions of a 3D array
    !   stored as an ASCII file. See the header of ASCII_read3Darray for
    !   more information.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   filename                i   char   name of the file
    !   ASCII_read3Darraysize   o   i8b    dimensions of the array
    !
    character(*), intent(in) :: filename
    integer(i8b) :: ASCII_read3Darraysize(3)

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   err   i8b   error number
    !   da    i8b   3-element vector containing the dimensions of the array
    !
    integer(i8b) :: err
    real(dp) :: da(3)

    ! *****************************************************************************
    !
    !   I/O
    !
    !   read dimensions of an array from an ASCII file containing a 3D array
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !
    ! =============================================================================

    open(unit=iImportUNIT,file=filename,status="OLD",form="FORMATTED",iostat=err)
    if (err/=0) then
        write(*,"('Error reading temporary file ',A,', error code ',I3)") trim(filename),err
        stop
    end if

    read(unit=iImportUNIT,iostat=err,fmt="(3ES25.16E3)") da

    ASCII_read3Darraysize = da

    close(unit=iImportUNIT)

    END FUNCTION ASCII_read3Darraysize

    SUBROUTINE remove_CR(string)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     090505  Original code (KH)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine remove_CR replaces any Carriage Return symbols (ASCII code
    !   13) by a space. This is necessary to prevent a mismatch between Fortran
    !   and C (or Fortran and Windows) line ending types.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   string   io   char   string to be processed
    !
    character(LEN=*), intent(inout) :: string

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   lenstring   i4b   length of the string
    !   i           i4b   temporary variable
    !
    integer(i4b) :: lenstring,i

    ! *****************************************************************************
    !
    !   I/O
    !
    !   none
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   none
    !

    lenstring=LEN_TRIM(string)

    do i=1,lenstring
        if (ichar(string(i:i))==13) then
            string(i:i)=' '
        end if
    end do

    END SUBROUTINE remove_CR

    SUBROUTINE ScattererInit()
    
    ScattererParams%sBubbleDir = 'Bubbles/'
    
	ScattererParams%time_norm = 1.0D0
	ScattererParams%rad_norm  = 1.0D0
	
    !------------ Experimental Values
    ScattererParams.coeff_fit =(/ -326138851.8775966167449951171875Q0, 3329551876.01073455810546875Q0, -15287331213.0629024505615234375Q0, 41570311455.0132293701171875Q0, -74140077882.8638153076171875Q0, 90617932686.12164306640625Q0, -76870588777.044830322265625Q0, 44688537547.7557373046875Q0, -17039124392.086650848388671875Q0, 3847686270.230488300323486328125Q0, -390758718.15616351366043090820312Q0 /)
    ScattererParams.A0c  = 0.0Q0; 
    END SUBROUTINE ScattererInit

    END MODULE ParnacParamInit
