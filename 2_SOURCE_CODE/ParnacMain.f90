    PROGRAM Parnac_Main



    !*******************************************************************************
    ! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     030221  Original code (KH - LD - AM)
    !
    !
    !   Copyright (C)   Laboratory of Electromagnetic Research
    !                   Delft University of Technology, Delft, The Netherlands
    !                   Laboratory of Wavefield Imaging
    !                   Delft University of Technology, Delft, The Netherlands
    !
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The program Parnac_Main contains the core of the INCS program. With this
    !   program, the acoustic pressure field from a plane source with an arbitrary
    !   aperture may be obtained through a nonlinear and inhomogeneous (in the attenuation
    !   in the coefficient of nonlinearity and in the speed of sound) medium.
    !   The code is written for a multi-processor system, where processor
    !   communication is performed through MPI, which requires the MPI library to
    !   be linked with the object files. The output files are in the HDF5 format,
    !   which requires the HDF5 library. For the Fast Fourier Transforms (FFT's),
    !   employed on more than one occasion in the code, we utilize the routines in
    !   the FFTW library, which is therefore also required.
    !
    !   Method: We employ a contrast source description and an iterative solution
    !   in the form of a Neumann, Bi-CGSTAB or Steepest Descent iterative scheme.
    !   Each scheme is has pro and cons. In general Neumann suits for modeling
    !   nonlinear propagation and moderate losses, Bi-CGSTAB allows for stronger
    !   contrasts and Steepest descent allows for inhomogeneous speed of sound.
    !
    !
    !   More information:
    !   - J. Huijssen, Modeling of Nonlinear Medical Diagnostic Ultrasound, PhD.
    !       thesis. Delft: Delft University of Technology, 2008. available from
    !       depository.tudelft.nl
    !   - J. Huijssen, Parnac - User manual.
    !   - L. Demi, Modeling nonlinear ultrasound through inhomogeneous biomedical
    !       tissue, PhD. thesis.
    !
    ! *****************************************************************************
    !
    !   MODULES USED
    	
	USE Types
    USE ParnacGeneral
    USE ParnacDataDef
    USE ParnacParamDef
    USE ParnacIdentifierDef

    USE ParnacDataSpaceInit, ONLY : &
        InitSpaceFromPhysicalDims, InitSpace, InitGrid, DestructSpace, SpaceInfo, &
        GridDistr0CreateEmpty, GridDistr0Deallocate, &
        iBeamOffsetT, iBeamOffsetX, iBeamOffsetY, &
        InitRRMSNorm, DestroyRRMSNorm, &
        InitRRMSPreviousCorrection, DestroyRRMSPreviousCorrection, &
        InitSaveSlices
    USE ParnacParamInit, ONLY : &
        ReadCommandLineAndEnvironment, ReadConfigFile, &
        ConfigFileInfo, NormalizeConfigParameters, &
        SpecialChecksConfigParameters
    USE ParnacIterStepMain, ONLY : &
        PrimarySourceToField, PlaneWaveField, &
        ContrastSourceToField, ContrastSourcetoFieldConj, FieldToContrastSource, FieldtoContrastSourceLin, & !! LSCM L.D. 12-07-2010
    CFPlaneSourceToField, &
        LoadField, StoreField, AddStoredToField, AddStoredtoFieldDFM, SubtractData, SumData, Copydata, SubtractStoredToFieldBis, SubtractStoredToFieldTris, sumsquare, sumandmultiply, sumandmultiplyCG, sumandmultiplyCGbis, sumandmultiplyCGtris, TwoRealValueCalculation, AbsoluteValueCalculation,  &
        LoadPlane, StorePlane, ZeroFirstPlane, &
        StoreSlices, UpdateRRMSError, ExportRRMSError, &
        test_isnan,PointSourceCloudField
    USE ParnacContrastFunctions, ONLY : &
        InitContrastSpace
    USE ParnacOutput, ONLY : &
        ExportSlice,Print_Error


    !! *****************************************************************************
    !!
    !!   INPUT/OUTPUT PARAMETERS
    !!
    !!   From the command line, the ParnacMain program should be invoked with
    !!
    !!   mpiexec -n 48 INCS.exe Inputfilename.in
    !!
    !!   in this example 48 cores have been used
    !! *****************************************************************************
    !!
    !!   LOCAL PARAMETERS - Examples
    !!
    !!   codeversion       dp    Current version of the code
    !!   inputfileversion  dp    Current version of the inpu config file
    !!   cRefBeam        type(Space)   Reference beam space
    !!   cBeam           type(Space)   Current beam space
    !!   cBeamCS         type(Space)   Previous beam space acting as contrast
    !!                                 source to current beam
    !!   cBeamST         type(Space)   Static Term beam space, for linearized contrast source method
    !!   cPlaneCF        type(Space)   Beam space containing the contrast field
    !!                                 plane summing up all contributions
    !!                                 of the previous beams except the last one
    !!   cInhomContrast  type(Space)   Beam space containing the inhomogeneity
    !!                                 contrast
    !!   cRRMSNorm       type(RRMSNormtype) structure containing the RRMS errors
    !!                                 of all beams and all iterations
    !!   iBeam             i8b         loop counter for the current beam
    !!   iIter             i8b         loop counter for the current iteration
    !!   i                 i8b         general loop counter
    !!   error             i8b         error number
    !!   iStartT           i8b         starting index of the beam in T
    !!   iStartX           i8b         starting index of the beam in X
    !!   iStartY           i8b         starting index of the beam in Y
    !!   iStartZ           i8b         starting index of the beam in Z
    !!   iProcN            i8b         number of processors
    !!   iProcID           i8b         identification number of the current proc
    !!   iProcN4           i4b         same as iProcN, but now in i4b form
    !!   iProcID4          i4b         same as iProcID, but now in i4b form
    !!   iErr              i4b         error number, in i4b form
    !!   acTemp            char        temporary char array for output log messages
    IMPLICIT NONE



    real(dp) ::				codeversion=1.3, inputfileversion=1.3
    real(8) ::              deltaepsilonNegative,deltaepsilon, ThirOrderPolinomialValueDeltaNegative, ThirOrderPolinomialValueDelta, ThirOrderPolinomialValue
    real(8) ::              AlphaCG, OmegaCG, RhoOld, RhoNew, BetaCG, RandomNumber,AlphaCoeffCGscheme,Phi2,Phi3,Phi4,TempPhi,z,q,acoefficient1,acoefficient2,singacoefficient1,singacoefficient2,TempResidualOutput1,TempResidualOutput2,TempResidualOutput3
    REAL(8), ALLOCATABLE :: TempResidualOutputVectorNum(:,:), TempResidualOutputVectorDen(:,:)
    
	type(Space)::           cRefBeam, cBeam, cBeamDFM, cBeamCS, cBeamST, cBeamTot, cResidual, cResidualNKCDX,cResidualNKCDY,cResidualNKCDZ, cPCG, cTemp, cTemp1,cTemp2,cTemp3,cTemp4, cP, cXCG, cV, cPlaneCF, cT, cInhomContrast, cDirection, cBeamF, cTheta1, cTheta2, cResidualNeumann1, cResidualNeumann2
    type(RRMSNormtype) ::	cRRMSNorm
	
    integer(i8b) ::			iBeam, iIter, i, error, indexResidualVectorNumDen
    integer(i8b) ::			iStartT,iStartX,iStartY,iStartZ, LCSM, BiCGSTAB, DFM, CG, ResidualNeumann, NKC, LINEARCASEIMPLEMENTATION
    integer(i8b) ::			iProcN, iProcID, IndexBiCG, indicerefinedapproximation
    integer(i4b) ::		    iProcN4, iProcID4, iErr; !For interaction with MPI commands
    character(LEN=1024)::   acTemp;
    Integer (4) :: err_unit_number
    Character (150) :: error_file, error_file2
    character(len=10) iProcIDString
    
    integer,parameter :: 	REQUIRED = MPI_THREAD_MULTIPLE 	!Added by A.M 13/09/2020
    integer ::			PROVIDED				!Added by A.M 13/09/2020


    ! *****************************************************************************
    !
    !   I/O
    !
    !   Text output to screen and to log files, input/output file communication
    !
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   See the USE, ONLY for a complete list of subroutines and functions called.
    !
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   Initialize: MPI, input parameters, log files, stopwatches, reference beam,
    !     save slices and RRMS norm, inhomogeneity contrast (if applicable)
    !   For each beam in the field space do
    !       Initialize beam
    !       Calculate the primary field solution
    !       If required, output primary field solution to file
    !       Store primary field solution to temp storage on disk
    !       If the beam is not the first one
    !           Get the previous beam from the temp storage
    !           Obtain the contrast source for the previous beam
    !           Calculate beam-to-beam solution of the previous beam to this beam
    !           Store previous beam contribution on temp storage
    !           If the beam is not the first or the second one
    !               Get the first plane of the previous beam CONTRAST from the
    !                 temp storage
    !               Calculate contrast plane-to-beam solution that sums up the
    !                 contrasts of all previous beams except the last one -
    !                 more or less of Huygens' principle, except that we do not
    !                 use a closed surface, and we use it only for the contrast
    !                 field, not the total field.
    !               Add contribution of the previous beam to this and store total
    !                 as the previous beam contribution on temp storage
    !           Add primary field solution to previous beam contributions
    !       If the number of iterations is larger than 1
    !       For each iteration do
    !           Obtain the contrast source from the current beam solution
    !           Obtain the contrast field from the contrast source
    !           Obtain the RRMS error of the contrast field
    !           If the beam is not the first one
    !               Add the previous beam contributions from temp storage
    !                 to the contrast field
    !           If the total number of beams is larger than two
    !               Store the first plane of the contrast field to temp storage
    !           Add the primary field solution to the contrast field to obtain
    !             the new beam solution
    !           If required, output the new beam solution to file
    !       Store the beam solution to temp storage
    !       Export RRMS errors
    !   De-initialize: de-allocate spaces and variables, close log files and end
    !     MPI
    !
    ! =============================================================================
    !
    !------------------------------------------
    ! General Initialization Phase
    !------------------------------------------

    ! call MPI_INIT_THREAD(REQUIRED, PROVIDED, iErr)	!Added by A.M 13/09/2020 [before MPI_INIT(iErr)]
    call MPI_INIT(iErr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, iProcN4, iErr);
    call MPI_COMM_RANK(MPI_COMM_WORLD, iProcID4, iErr);
    
    iProcN = iProcN4
    iProcID = iProcID4

    iMemAllocated=0

    !Welcome message
    write ( *, '("Parnac V",F5.2," - Processor ", I3, " checking in.")') codeversion,iProcID
    write ( *, '(" ")')

    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

    ! Read the command line
    call ReadCommandLineAndEnvironment(iProcID)

    ! Start the stopwatches and the initialization stopwatch
    call StopWatchesInit(iProcID);
    call SWStart(cswIniti)

    ! Read the parameters
    call ReadConfigFile(inputfileversion,error)

    !------------------------------------------
    !------------SWITCH SECTION----------------
    !------------------------------------------

    !nonlinear coeff in kappa on(1) or of(0) [Note: implemented only for C.G. scheme]
    NKC=cSourceParams.NKC

    !LCSM(1) or Standard Neumann(0)
    LCSM=cSourceParams.LCSM             !LCSM L.D. 13-07-2010

    !After LCSM has been set to 1
    !LCSM Bi-CGSTAB(1) or LCSM Neumann(0)
    BiCGSTAB=cSourceParams.BiCGSTAB     !LCSM L.D. 17-09-2010

    !C.G. Scheme [Steepest Descent]
    CG=cSourceParams.CG                 !Steepest Descent L.D. 04-07-2011

    ! If we use a Steepest descent with compleatly linear media this has to be set to 1 because of the alpha parameter
    LINEARCASEIMPLEMENTATION=cSourceParams.LINEARCASEIMPLEMENTATION

    !DFM(1) enables double frequency mode and
    ! generates the stored data Beamsolo2
    !DFM(2) combines the present configuration
    ! with BeamSol2
    !DFM(3) combines the present configuration
    ! with BeamSol2 and save the present in BeamSol2
    DFM=cSourceParams.DFM    !DFM L.D. 05-01-2011

    ResidualNeumann=0
    ! if we want to have the residual as an output set this flag to 1
    !------------------------------------------

    if(iProcID == 0) then
        call ConfigFileInfo(inputfileversion)
    end if
    if(error/=0) then
        stop 'Error: reading from input files unsuccessful.'
    end if
    call SpecialChecksConfigParameters
    call NormalizeConfigParameters

    ! Initialize log and stopwatch, write startup messages to screen and logfile
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !For clock synchronization in OpenLogFile

    call OpenLogFile(iProcID)
    write (acTemp, '("Parnac V",F5.2," - Processor ", I3, " checking in.")') codeversion,iProcID; call PrintToLog(acTemp, 0);
    write (acTemp, '("Input file is ", A )') trim(sInputName); call PrintToLog(acTemp, 0);
    write (acTemp, '("Input directory is ", A )') trim(sInputDir); call PrintToLog(acTemp, 0);
    write (acTemp, '("Output directory is ", A )') trim(sOutputDir); call PrintToLog(acTemp, 0);
    write (acTemp, '("Temporary directory is ", A )') trim(sTempDir); call PrintToLog(acTemp, 0);
    if (cModelParams.UseYSymmetry) then
        call PrintToLog("Symmetry in the Y-dimension observed, grid reduction is employed.",0)
    else
        call PrintToLog("No symmetry in the Y-dimension observed, no grid reduction.",0)
    end if
    
    ! Initialize an empty beam acting as a model for the others
    call InitSpaceFromPhysicalDims(cRefBeam, iSI_FIELD, cModelParams.UseYSymmetry, &
        cModelParams.Lx, cModelParams.Ly, cModelParams.Lt, &
        cModelParams.Sz, cModelParams.Lz/cModelParams.Numbeams, &
        cModelParams.Thetax, cModelParams.Thetay, cModelParams.Thetat, &
        cModelParams.FNyq);        ! Define the part of the space the Near Field refers to
    call InitGrid(cRefBeam, iProcN, iProcID, (/ .true., .true., .true., .true./));                                ! Define the structure to store the data on the near-field

    ! Display some debugging-information
    call PrintToLog("Variables in cBeam",-1)
    call SpaceInfo(cRefBeam)

    ! Initialize save slices
    call InitSaveSlices(cRefBeam) !calculate where in the beams the different save slices are located

    ! Initializes the cRRMSNorm storage array and previous correction slice
    ! RRMS Norm evaluation: only from the last z slice in a beam
    ! If this processor carries some positions in that slice:
    ! Calculate the index start of the indices to be stored in that grid
    ! Initialize an array that acts as storage for the RRMS Norm evaluation
    call InitRRMSNorm(cRRMSNorm, iProcID, cModelParams%Numbeams, maxval(cModelParams%Numiterations))
    call InitRRMSPreviousCorrection(cRRMSNorm,cRefBeam)

    ! Initialize and export the inhomogeneous contrast source space
    if (cModelParams.ContrastSourceType/=iCI_NONLIN .AND. cModelParams.ContrastSourceType/=iCI_SCATTERER ) then

        call InitSpace(cInhomContrast, iSI_INHOMCONTRAST, cRefBeam%bYSymm, &
            1_i8b, cRefBeam.iDimX, cRefBeam.iDimY, cRefBeam.iDimZ,  &
            0_i8b, cRefBeam.iStartX, cRefBeam.iStartY, cRefBeam.iStartZ,  &
            0_i8b, 0_i8b, 0_i8b, &
            cRefBeam%dFnyq, cRefBeam.dTanX, cRefBeam.dTanY, cRefBeam.dTanT )
        call InitGrid(cInhomContrast, cRefBeam%cGrid%iProcN, cRefBeam%cGrid%iProcID, (/ .false., .false., .false., .false./));                                ! Define the structure to store the data on the near-field

        call InitContrastSpace(cInhomContrast,cModelParams.ContrastSourceType)

        call ExportSlice(trim(sOutputDir) // trim(sOutputRoot) // "contrast",&
            "Contrast", &
            cInhomContrast, &
            (/ 0_i8b, 0_i8b, 0_i8b, 0_i8b /), &
            (/ 1_i8b, cInhomContrast.iDimX, cInhomContrast.iDimY, cInhomContrast.iDimZ /), &
            cModelParams.xyzsliceindex, cInhomContrast.iDimT, .true.);

        if (cModelParams.UseAntiAliasing) then
            ! Contrast is stored on disk, don't need the grid anymore..
            call GridDistr0DeAllocate(cInhomContrast%cGrid);
        end if

    end if

    call SWStop(cswIniti)
    !------------------------------------------
    ! End of General Initialization Phase
    !------------------------------------------

    !------------------------------------------
    ! Start of the Pre-loop Phase
    !------------------------------------------
    do iBeam = 0, cModelParams%NumBeams-1

        call PrintToLog("**********************************************************",0)
        write(acTemp,'("****** Calculate Iterative solution for Beam ",I3," ******")') iBeam
        call PrintToLog(acTemp,0)
        call PrintToLog("**********************************************************",0)

        write(acTemp,'("Initialize Beam ",I3)') iBeam; call PrintToLog(acTemp,0)

        ! Obtain the new starting positions for this beam
        iStartT = cRefBeam.iStartT + iBeamOffsetT(cRefBeam, iBeam*cRefBeam%iDimZ)
        iStartX = cRefBeam.iStartX + iBeamOffsetX(cRefBeam, iBeam*cRefBeam%iDimZ)
        iStartY = cRefBeam.iStartY + iBeamOffsetY(cRefBeam, iBeam*cRefBeam%iDimZ)
        iStartZ = cRefBeam.iStartZ + iBeam*cRefBeam%iDimZ

        ! Initialize the field beam
        call InitSpace(cBeam, iSI_FIELD, cRefBeam%bYSymm, &
            cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
            iStartT, iStartX, iStartY, iStartZ,  &
            0_i8b, 0_i8b, 0_i8b, &
            cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
        call InitGrid(cBeam, iProcN, iProcID, (/ .false., .false., .false., .false./));
        ! Define the structure to store the data on the near-field

        !-------------------------------------------------------------------------------------
        !Residual L.D. 13-07-2010 - Initialization for Extra Residuals calculation
        !-------------------------------------------------------------------------------------
        if (ResidualNeumann==1) then

            call InitSpace(cResidualNeumann1, iSI_FIELD, cRefBeam%bYSymm, &
                cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
                iStartT, iStartX, iStartY, iStartZ,  &
                0_i8b, 0_i8b, 0_i8b, &
                cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
            call InitGrid(cResidualNeumann1, iProcN, iProcID, (/ .false., .false., .false., .false./));

            call InitSpace(cResidualNeumann2, iSI_FIELD, cRefBeam%bYSymm, &
                cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
                iStartT, iStartX, iStartY, iStartZ,  &
                0_i8b, 0_i8b, 0_i8b, &
                cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
            call InitGrid(cResidualNeumann2, iProcN, iProcID, (/ .false., .false., .false., .false./));

        end if

        !-------------------------------------------------------------------------------------
        !DFM L.D. 13-07-2010 - Initialization for Dual Frequency Modality
        !-------------------------------------------------------------------------------------
		if (DFM==2 .OR. DFM==3) then
		        call InitSpace(cBeamDFM, iSI_FIELD, cRefBeam%bYSymm, &
		            cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		            iStartT, iStartX, iStartY, iStartZ,  &
		            0_i8b, 0_i8b, 0_i8b, &
		            cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		        call InitGrid(cBeamDFM, iProcN, iProcID, (/ .false., .false., .false., .false./));
		end if

        !-------------------------------------------------------------------------------------
        !LCSM L.D. 13-07-2010 - Initialization for Linearized Contrast Source Method
        !-------------------------------------------------------------------------------------
		if (LCSM==1) then
		    call InitSpace(cBeamST, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cBeamST, iProcN, iProcID, (/ .false., .false., .false., .false./));
		
		    call InitSpace(cBeamTot, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cBeamTot, iProcN, iProcID, (/ .false., .false., .false., .false./));
		end if
    
    

        !--------------------------------------------------------------------------------------
        !Bi-CGSTAB L.D. 17-09-2010 - Initialization for Bi-CGSTAB
        !-------------------------------------------------------------------------------------
        if (BiCGSTAB==1) then
		    call InitSpace(cResidual, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cResidual, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cPCG, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cPCG, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cV, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cV, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cP, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cP, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cXCG, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cXCG, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cT, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cT, iProcN, iProcID, (/ .false., .false., .false., .false./));

    	end if	

        !--------------------------------------------------------------------------------------
        !Steepest Descent L.D. 20-07-2011 - Initialization for Steepest Descent
        !-------------------------------------------------------------------------------------
        if (CG==1) then
		    call PrintToLog("****************************************************************",0)
		    call PrintToLog("****** Initialising Steepest Descent scheme ********************************",0)
		    call PrintToLog("****************************************************************",0)

		    call InitSpace(cResidual, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cResidual, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    if (NKC==1) then
		        call InitSpace(cResidualNKCDX, iSI_FIELD, cRefBeam%bYSymm, &
		            cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		            iStartT, iStartX, iStartY, iStartZ,  &
		            0_i8b, 0_i8b, 0_i8b, &
		            cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		        call InitGrid(cResidualNKCDX, iProcN, iProcID, (/ .false., .false., .false., .false./));

		        call InitSpace(cResidualNKCDY, iSI_FIELD, cRefBeam%bYSymm, &
		            cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		            iStartT, iStartX, iStartY, iStartZ,  &
		            0_i8b, 0_i8b, 0_i8b, &
		            cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		        call InitGrid(cResidualNKCDY, iProcN, iProcID, (/ .false., .false., .false., .false./));

		        call InitSpace(cResidualNKCDZ, iSI_FIELD, cRefBeam%bYSymm, &
		            cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		            iStartT, iStartX, iStartY, iStartZ,  &
		            0_i8b, 0_i8b, 0_i8b, &
		            cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		        call InitGrid(cResidualNKCDZ, iProcN, iProcID, (/ .false., .false., .false., .false./));
		    end if

		    call InitSpace(cDirection, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cDirection, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cBeamF, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cBeamF, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp1, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp1, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp2, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp2, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp3, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp3, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTemp4, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTemp4, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTheta1, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTheta1, iProcN, iProcID, (/ .false., .false., .false., .false./));

		    call InitSpace(cTheta2, iSI_FIELD, cRefBeam%bYSymm, &
		        cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
		        iStartT, iStartX, iStartY, iStartZ,  &
		        0_i8b, 0_i8b, 0_i8b, &
		        cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
		    call InitGrid(cTheta2, iProcN, iProcID, (/ .false., .false., .false., .false./));

		end if
        !--------------------------------------------------------------------------------------
        ! Calculate the primary field solution

        call PrintToLog("*************************************************************",0)
        call PrintToLog("Calculate primary field solution",0)
        call PrintToLog("*************************************************************",0)

        ! Either from source description or as plane wave, this function generates the primary field
        ! given the source type defined in the input file

        if (cModelParams%PrimarySourceType == iEI_LOADFIELD) then
            call GridDistr0CreateEmpty(cBeam.cGrid); 
            cBeam%iSpaceIdentifier = iSI_FIELD; 
            call LoadField(cBeam, trim(cModelParams.FieldFilename))
        elseif (cModelParams%PrimarySourceType /= iEI_PLANEWAVE .AND. cModelParams%PrimarySourceType /= iEI_POINTSOURCECLOUD) then
            call PrimarySourcetoField(cBeam); 
        elseif (cModelParams%PrimarySourceType == iEI_PLANEWAVE) then
            call PlanewaveField(cBeam); 
        elseif (cModelParams%PrimarySourceType == iEI_POINTSOURCECLOUD) then
            call PointSourceCloudField(cBeam); 
        end if

        if (ResidualNeumann==1) then 

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams%PrimarySourceType/=iEI_POINTSOURCECLOUD) then
                call PrimarySourcetoField(cResidualNeumann1);
            elseif (cModelParams%PrimarySourceType==iEI_PLANEWAVE) then
                call PlanewaveField(cResidualNeumann1);
			elseif(cModelParams%PrimarySourceType==iEI_POINTSOURCECLOUD) then
				call PointSourceCloudField(cBeam);
            end if

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams%PrimarySourceType/=iEI_POINTSOURCECLOUD) then
                call PrimarySourcetoField(cResidualNeumann2);
            elseif (cModelParams%PrimarySourceType==iEI_PLANEWAVE) then
                call PlanewaveField(cResidualNeumann2);
			elseif(cModelParams%PrimarySourceType==iEI_POINTSOURCECLOUD) then
				call PointSourceCloudField(cBeam);
            end if

        end if

        if (DFM == 2 .OR. DFM == 3) then

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE .AND. cModelParams%PrimarySourceType/=iEI_POINTSOURCECLOUD) then
                call PrimarySourcetoField(cBeamDFM);
            elseif (cModelParams%PrimarySourceType==iEI_PLANEWAVE) then
                call PlanewaveField(cBeamDFM);
			elseif(cModelParams%PrimarySourceType==iEI_POINTSOURCECLOUD) then
				call PointSourceCloudField(cBeam);
            end if

        end if

        !---------------------------------------------------------------------------------------------------------
        ! Start pre-loop Steepest Descent scheme - Initialization of the starting values for the different spaces
        !---------------------------------------------------------------------------------------------------------
        if (CG==1) then

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cResidual);
				else
				    call PlanewaveField(cResidual);
				end if

				if (NKC==1) then
				    if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				        call PrimarySourcetoField(cResidualNKCDX);
				    else
				        call PlanewaveField(cResidualNKCDX);
				    end if

				    if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				        call PrimarySourcetoField(cResidualNKCDY);
				    else
				        call PlanewaveField(cResidualNKCDY);
				    end if
				    if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				        call PrimarySourcetoField(cResidualNKCDZ);
				    else
				        call PlanewaveField(cResidualNKCDZ);
				    end if
				end if

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cDirection);
				else
				    call PlanewaveField(cDirection);
				end if


				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cBeamF);
				else
				    call PlanewaveField(cBeamF);
				end if


				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTheta1);
				else
				    call PlanewaveField(cTheta1);
				end if


				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTheta2);
				else
				    call PlanewaveField(cTheta2);
				end if


				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTemp);
				else
				    call PlanewaveField(cTemp);
				end if

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTemp1);
				else
				    call PlanewaveField(cTemp1);
				end if

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTemp2);
				else
				    call PlanewaveField(cTemp2);
				end if

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTemp3);
				else
				    call PlanewaveField(cTemp3);
				end if

				if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
				    call PrimarySourcetoField(cTemp4);
				else
				    call PlanewaveField(cTemp4);
				end if

				! test functions in order to make sure that the fields have no NaN value
				call test_isnan(cResidual)
				call test_isnan(cDirection)
				call test_isnan(cBeamF)
				call test_isnan(cTemp)
				call test_isnan(cTemp1)
				call test_isnan(cTemp2)
				call test_isnan(cTemp3)
				call test_isnan(cTemp4)
				call test_isnan(cTheta1)
				call test_isnan(cTheta2)

				call PrintToLog("*************************************************************",0)
				call PrintToLog("****** Calculating Residual0 ********************************",0)
				call PrintToLog("*************************************************************",0)
				! Calculate residual0

				! this fuction calculates the contrast given the field [cResidual in this case] and put
				! the result in the space where the field was. The type of contrast is defined via the
				! variable cModelParams.ContrastSourceType which is an integer
				call FieldtoContrastSource(cResidual, cModelParams.ContrastSourceType, cInhomContrast)

				call PrintToLog("****** Introducing the contrast in kappa ********************************",0)

				! this function calculates the field given the contrast. It basically convolves the contrast with the Green's function.
				call ContrastSourcetoField(cResidual, cResidual, .false.)

				if (NKC==1) then
				    ! cResidualNKCDX/DY/DZ contain P0

				    call FieldtoContrastSource(cResidualNKCDX, 10, cInhomContrast)
				    ! cResidualNKCDX contains dx((dxkappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
				    ! cResidualNKCDX contains G conv (dx((dxkappa)P0^2)))

				    call FieldtoContrastSource(cResidualNKCDY, 11, cInhomContrast)
				    ! cResidualNKCDY contains dy((dykappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
				    ! cResidualNKCDY contains G conv (dy((dykappa)P0^2))

				    call FieldtoContrastSource(cResidualNKCDZ, 12, cInhomContrast)
				    ! cResidualNKCDZ contains dz((dzkappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)
				    ! cResidualNKCDZ contains G conv (dz((dzkappa)P0^2))

				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
				    ! this function does A=B-d*C
				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
				    ! this function does A=B-d*C

				    RandomNumber=-0.5;
				    call SubtractStoredToFieldBis(cResidual,cResidual,cResidualNKCDX,RandomNumber);
				    ! this function does A=B-d*C
				    !The standard Residual schould be summed with it

				end if

				call MPI_BARRIER(MPI_COMM_WORLD, iErr);                                           !L.D. 07-11-2011
				call AbsoluteValueCalculation(TempResidualOutput1,cResidual)                      !L.D. 07-11-2011
				call AbsoluteValueCalculation(TempResidualOutput2,cBeam)                          !L.D. 07-11-2011


				call MPI_BARRIER(MPI_COMM_WORLD, iErr);
				TempResidualOutput3=TempResidualOutput1/TempResidualOutput2

				!---- Export the residuals ------------------------------------------------
				Write (iProcIDString, '(i10.10)') iProcID

				err_unit_number = 10
				error_file = trim(cSourceParams.residualdirectory)//'errorN'//iProcIDString //'.dat'
				Open (Unit=err_unit_number, File=trim(error_file), Status='replace')
				Close (Unit=err_unit_number)
				Call print_error (TempResidualOutput2, err_unit_number, error_file)
				Call print_error (TempResidualOutput1, err_unit_number, error_file)

				call MPI_BARRIER(MPI_COMM_WORLD, iErr);

				RandomNumber=0
				call SubtractStoredToFieldBis(cTemp1,cResidual,cResidual,RandomNumber);
				! this function does A=B-d*C
				cTemp1%iSpaceIdentifier=iSI_FIELDSLICE

				call PrintToLog("**************************************************************",0)
				call PrintToLog("****** Calculating Direction0 ********************************",0)
				call PrintToLog("**************************************************************",0)

				RandomNumber=0
				call SubtractStoredToFieldBis(cDirection,cResidual,cDirection,RandomNumber)             !! L.D. 29-08-2011
				! this function does A=B-d*C

				call ContrastSourcetoFieldConj(cTemp1, cDirection, .false.)                             !! Introduce the conj G

				RandomNumber=0
				call SubtractStoredToFieldBis(cTemp4,cDirection,cDirection,RandomNumber)                !! L.D. 29-08-2011
				! this function does A=B-d*C

				if (NKC==1) then
				    RandomNumber=0

				    call SubtractStoredToFieldBis(cResidualNKCDX,cDirection,cDirection,RandomNumber)        !! L.D. 29-08-2011
				    ! this function does A=B-d*C
				    call SubtractStoredToFieldBis(cResidualNKCDY,cDirection,cDirection,RandomNumber)        !! L.D. 29-08-2011
				    ! this function does A=B-d*C
				    call SubtractStoredToFieldBis(cResidualNKCDZ,cDirection,cDirection,RandomNumber)        !! L.D. 29-08-2011
				    ! this function does A=B-d*C

				end if

				call FieldtoContrastSourceLin(cDirection, cDirection, 7, cInhomContrast)
				call FieldtoContrastSourceLin(cTemp4, cTemp4, 8, cInhomContrast)                !!L.D. 29-08-2011
				cDirection%iSpaceIdentifier=iSI_FIELD
				cTemp4%iSpaceIdentifier=iSI_FIELD                                               !!L.D. 29-08-2011

				call sumandmultiplyCGtris(cDirection,cBeamF,cResidual)
				! this function does A=A*B-C

				RandomNumber=-1                                                                 !!L.D. 29-08-2011
				call SubtractStoredToFieldBis(cDirection,cDirection,cTemp4,RandomNumber);       !!L.D. 29-08-2011
				! this function does A=B-d*C                                                    !!L.D. 29-08-2011

				if (NKC==1) then

				    ! cResidualNKCDX-DY-DZ contain G^H(r(p0)), cBeamF contains p0
				    call FieldtoContrastSourceLin(cResidualNKCDX,cBeamF, 16, cInhomContrast)
				    cResidualNKCDX%iSpaceIdentifier=iSI_FIELD
				    call FieldtoContrastSourceLin(cResidualNKCDY,cBeamF, 17, cInhomContrast)
				    cResidualNKCDY%iSpaceIdentifier=iSI_FIELD
				    call FieldtoContrastSourceLin(cResidualNKCDZ,cBeamF, 18, cInhomContrast)
				    cResidualNKCDZ%iSpaceIdentifier=iSI_FIELD

				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
				    !this function does A=B-d*C
				    RandomNumber=1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
				    !this function does A=B-d*C

				    RandomNumber=1;
				    call SubtractStoredToFieldBis(cDirection,cDirection,cResidualNKCDX,RandomNumber);
				    !this function does A=B-d*C

				end if

				RandomNumber=0
				call SubtractStoredToFieldBis(cTemp2,cDirection,cDirection,RandomNumber);
				! this function does A=B-d*C

				call PrintToLog("**********************************************************",0)
				call PrintToLog("****** Calculating Theta1 ********************************",0)
				call PrintToLog("**********************************************************",0)


				! Calculate Theta1
				call FieldtoContrastSourceLin(cDirection, cTheta1, 6, cInhomContrast)           !!L.D. 29-08-2011

				if (NKC==1) then
				    !cResidualNKCDX needs to be equal to cTheta1=p0 here
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDX,cDirection,cTheta1,RandomNumber)
				    ! this function does A=B-d*C
				    !cResidualNKCDY needs to be equal to cTheta1=p0 here
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDY,cDirection,cTheta1,RandomNumber)
				    ! this function does A=B-d*C
				    !cResidualNKCDZ needs to be equal to cTheta1=p0 here
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDZ,cDirection,cTheta1,RandomNumber)
				    ! this function does A=B-d*C

				    call FieldtoContrastSourceLin(cResidualNKCDX, cTheta1, 13, cInhomContrast)
				    call FieldtoContrastSourceLin(cResidualNKCDY, cTheta1, 14, cInhomContrast)
				    call FieldtoContrastSourceLin(cResidualNKCDZ, cTheta1, 15, cInhomContrast)

				    call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
				    call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
				    call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)

				    RandomNumber=-1
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber)
				    ! this function does A=B-d*C
				    RandomNumber=-1
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber)
				    ! this function does A=B-d*C

				end if

				RandomNumber=0
				call SubtractStoredToFieldBis(cTheta1,cTemp2,cTemp2,RandomNumber);              !!L.D. 29-08-2011
				! this function does A=B-d*C

				call ContrastSourcetoField(cDirection, cDirection, .false.)                     !!L.D. 29-08-2011

				RandomNumber=1
				call SubtractStoredToFieldBis(cTheta1,cDirection,cTheta1,RandomNumber);          !!L.D. 29-08-2011
				! this function does A=B-d*C

				if (NKC==1) then

				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cTheta1,cResidualNKCDX,RandomNumber);
				    ! this function does A=B-d*C

				end if

				RandomNumber=0
				call SubtractStoredToFieldBis(cDirection,cTemp2,cTemp2,RandomNumber);            !!L.D. 29-08-2011
				! this function does A=B-d*C

				!!!!----------------- Note that: Here above theta1 and cDirection have been swaaped

				call PrintToLog("**********************************************************",0)
				call PrintToLog("****** Calculating Theta2 ********************************",0)
				call PrintToLog("**********************************************************",0)

				!cTheta2 needs to be equal to cDirection
				RandomNumber=0
				call SubtractStoredToFieldBis(cTheta2,cDirection,cDirection,RandomNumber);
				! this function does A=B-d*C

				if (NKC==1) then
				    !cResidualNKCDX needs to be equal to cDirection
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDX,cDirection,cDirection,RandomNumber);
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDY,cDirection,cDirection,RandomNumber);
				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNKCDZ,cDirection,cDirection,RandomNumber);
				end if

				! Calculate Theta2
				call FieldtoContrastSourceLin(cTheta2, cTemp2, 9, cInhomContrast)           !!L.D. 29-08-2011
				!call FieldtoContrastSource(cTheta2, cModelParams.ContrastSourceType, cInhomContrast)
				call test_isnan(cTheta2)
				call ContrastSourcetoField(cTheta2, cTheta2, .false.)

				if (NKC==1) then

				    call FieldtoContrastSource(cResidualNKCDX, 10, cInhomContrast)
				    ! cResidualNKCDX contains dx((dxkappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
				    ! cResidualNKCDX contains G conv (dx((dxkappa)P0^2)))

				    call FieldtoContrastSource(cResidualNKCDY, 11, cInhomContrast)
				    ! cResidualNKCDY contains dy((dykappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
				    ! cResidualNKCDY contains G conv (dy((dykappa)P0^2))

				    call FieldtoContrastSource(cResidualNKCDZ, 12, cInhomContrast)
				    ! cResidualNKCDZ contains dz((dzkappa)P0^2)
				    call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)
				    ! cResidualNKCDZ contains G conv (dz((dzkappa)P0^2))

				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
				    ! this function does A=B-d*C
				    RandomNumber=-1;
				    call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
				    ! this function does A=B-d*C

				    RandomNumber=-0.5;
				    call SubtractStoredToFieldBis(cResidual,cResidual,cResidualNKCDX,RandomNumber);
				    ! this function does A=B-d*C
				    !The standard Residual schould be summed with it

				end if

				call PrintToLog("**********************************************************",0)
				call PrintToLog("****** Calculating Alpha0 ********************************",0)
				call PrintToLog("**********************************************************",0)

				! Calculate Phi2
				call TwoRealValueCalculation(Phi2,cResidual,cTheta1)
				call AbsoluteValueCalculation(TempPhi,cTheta2)
				Phi2=Phi2/TempPhi
				Phi2=Phi2/4

				! Calculate Phi3
				call TwoRealValueCalculation(Phi3,cResidual,cTheta2)
				call AbsoluteValueCalculation(TempPhi,cTheta1)
				Phi3=Phi3+TempPhi
				call AbsoluteValueCalculation(TempPhi,cTheta2)
				Phi3=Phi3/TempPhi
				Phi3=Phi3/2

				! Calculate Phi4
				call TwoRealValueCalculation(Phi4,cTheta1,cTheta2)

				Phi4=Phi4/TempPhi

				Phi4=Phi4*3

				Phi4=Phi4/4

				! Calculate alpha0
				z=(((Phi3*Phi4)-(3*Phi2))/6-(Phi4**3)/27)

				q=(Phi3/3)-((Phi4**2)/9)

				acoefficient1=(z)+sqrt((q**3)+(z**2))
				acoefficient2=(z)-sqrt((q**3)+(z**2))

				if (acoefficient1>0) then
				    singacoefficient1=+1;
				else
				    singacoefficient1=-1;
				end if

				if (acoefficient2>0) then
				    singacoefficient2=+1;
				else
				    singacoefficient2=-1;
				end if

				acoefficient1 = (singacoefficient1)*((acoefficient1*singacoefficient1)**(1./3.))
				acoefficient2 = (singacoefficient2)*((acoefficient2*singacoefficient2)**(1./3.))

				AlphaCoeffCGscheme=acoefficient1+acoefficient2-Phi4/3;
				if (isNAN(AlphaCoeffCGscheme)) then

				    call PrintToLog("****** Alpha0 is NAN ********************************",0)
				    AlphaCoeffCGscheme=-1;

				end if

				!---------------------------------------------------------------------------------
				!-------------------------START OPTIMIZZATION PROCESS-------------------------------
				!---------------------------------------------------------------------------------
				deltaepsilonNegative=-0.01;
				deltaepsilon=0.01;

				do indicerefinedapproximation=1, 1000

				    ThirOrderPolinomialValue=(AlphaCoeffCGscheme**3)+(Phi4*(AlphaCoeffCGscheme**2))+(Phi3*AlphaCoeffCGscheme)+(Phi2);
				    ThirOrderPolinomialValueDelta=((AlphaCoeffCGscheme+deltaepsilon)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilon)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilon))+(Phi2);
				    ThirOrderPolinomialValueDeltaNegative=((AlphaCoeffCGscheme+deltaepsilonNegative)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilonNegative)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilonNegative))+(Phi2);

				    do while (((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))).AND.((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))))

				        deltaepsilonNegative=deltaepsilonNegative*10;
				        deltaepsilon=deltaepsilon*10;
				        ThirOrderPolinomialValue=(AlphaCoeffCGscheme**3)+(Phi4*(AlphaCoeffCGscheme**2))+(Phi3*AlphaCoeffCGscheme)+(Phi2);
				        ThirOrderPolinomialValueDelta=((AlphaCoeffCGscheme+deltaepsilon)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilon)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilon))+(Phi2);
				        ThirOrderPolinomialValueDeltaNegative=((AlphaCoeffCGscheme+deltaepsilonNegative)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilonNegative)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilonNegative))+(Phi2);

				    end do

				    if ((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then

				        AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilon/2;
				        deltaepsilonNegative=-deltaepsilon/2;
				        deltaepsilon=deltaepsilon/2;

				    end if

				    if ((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then



				        AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilonNegative/2;
				        deltaepsilon=deltaepsilonNegative/2;
				        deltaepsilonNegative=-deltaepsilonNegative/2;



				    end if


				end do

				!---------------------------------------------------------------------------------
				!-------------------------END OPTIMIZZATION PROCESS-------------------------------
				!---------------------------------------------------------------------------------

				!---------------------------------------------------------------------------------
				!--------------------------------LINEAR CASE--------------------------------------
				!---------------------------------------------------------------------------------
				if (LINEARCASEIMPLEMENTATION==1) then

				    call PrintToLog("**********************************************************",0)
				    call PrintToLog("****** Steepest Descent Linear Case implementation********",0)
				    call PrintToLog("**********************************************************",0)

				    Phi2=0
				    Phi3=0
				    call TwoRealValueCalculation(Phi2,cResidual,cTheta1)
				    Phi2=Phi2/2
				    Phi3=Phi2
				    call AbsoluteValueCalculation(TempPhi,cTheta1)
				    Phi2=Phi2/TempPhi
				    AlphaCoeffCGscheme=-Phi2


				    !-------------------OPTIMIZATION PROCESS------------------------------------------

				    deltaepsilonNegative=-0.01;
				    deltaepsilon=0.01;

				    do indicerefinedapproximation=1, 1000

				        ThirOrderPolinomialValue=AlphaCoeffCGscheme*TempPhi+Phi3
				        ThirOrderPolinomialValueDelta=(AlphaCoeffCGscheme+deltaepsilon)*TempPhi+Phi3
				        ThirOrderPolinomialValueDeltaNegative=(AlphaCoeffCGscheme+deltaepsilonNegative)*TempPhi+Phi3

				        do while (((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))).AND.((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))))

				            deltaepsilonNegative=deltaepsilonNegative*10;
				            deltaepsilon=deltaepsilon*10;
				            ThirOrderPolinomialValue=AlphaCoeffCGscheme*TempPhi+Phi3
				            ThirOrderPolinomialValueDelta=(AlphaCoeffCGscheme+deltaepsilon)*TempPhi+Phi3
				            ThirOrderPolinomialValueDeltaNegative=(AlphaCoeffCGscheme+deltaepsilonNegative)*TempPhi+Phi3

				        end do

				        if ((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then

				            AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilon/2;
				            deltaepsilonNegative=-deltaepsilon/2;
				            deltaepsilon=deltaepsilon/2;

				        end if

				        if ((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then



				            AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilonNegative/2;
				            deltaepsilon=deltaepsilonNegative/2;
				            deltaepsilonNegative=-deltaepsilonNegative/2;



				        end if


				    end do


				end if
				!---------------------------------------------------------------------------------
				!------------------------------END LINEAR CASE------------------------------------
				!---------------------------------------------------------------------------------

				call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

			end if

        !-------------------------------------------------------------------------------------
        ! start pre-loop full nonlinear neumann - Linear incident field generation
        if (LCSM==0) then
            ! Either from source description or as plane wave
            ! Output slices if required
            if((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
                .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
            call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(0_i8b), &
                "p0", cBeam, iBeam, .true.)
            end if

        end if

        ! Test whether there is a NaN in cBeam - indication of an error,
        ! no use to continue the program
        call test_isnan(cBeam)

        !If necessary, store the primary field solution for use in iterative scheme- note that in case
        ! Dual Frequency Modality (DFM) or linearization (LCSM) are on other fields have to be stored
        if(cModelParams%Numiterations(1+iBeam)>1) then
			if (cModelParams%PrimarySourceType /= iEI_LOADFIELD) call StoreField(cBeam, "BeamSol0")

            if (DFM==1) then

                call StoreField(cBeam,"BeamSol2")

            end if

            if (DFM==2 .OR. DFM ==3) then

                call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
                call LoadField(cBeamDFM,"BeamSol2")
                call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
                call AddStoredtoFieldDFM(cBeam,cBeamDFM);
                call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
                call StoreField(cBeam,"BeamSol0")
                
                if (DFM==3) then
                	call StoreField(cBeam,"BeamSol2")
                end if

            end if
            


        end if

		!-------------------------------------------------------------------------------------
		! start Linearization pre-loop code L.D. 20-07-2011
		!-------------------------------------------------------------------------------------

		if (LCSM==1) then
			    ! Either from source description or as plane wave
			    if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
			    call PrimarySourcetoField(cBeamST);
			    else
			    call PlanewaveField(cBeamST);
			    end if

			    if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
			    call PrimarySourcetoField(cBeamTot);
			    else
			    call PlanewaveField(cBeamTot);
			    end if

			    call test_isnan(cBeamTot)
			    call test_isnan(cBeamST)

			    call FieldtoContrastSource(cBeamST, cModelParams.ContrastSourceType, cInhomContrast)
			    call test_isnan(cBeamST)
			    call ContrastSourcetoField(cBeamST, cBeamST, .false.)


			    !If necessary, store the staticterm field solution for use in iterative scheme
			    if(cModelParams%Numiterations(1+iBeam)>1) then
			    call StoreField(cBeamST,"BeamStatic")
			    end if

			    !-------------------------------------------------------------------------------------
			    ! start Bi-CGSTAB pre-loop code L.D. 20-07-2011
			    !-------------------------------------------------------------------------------------
			    if (BiCGSTAB==1) then

					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cResidual);
					else
						call PlanewaveField(cResidual);
					end if
					call test_isnan(cResidual)

					! TEST EXPORT DATA
					!=============================================================================
					!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(2222), &
					!				"p0", cResidual, iBeam, .true.)
					!=============================================================================
					! TEST EXPORT DATA
					!=============================================================================
					!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(5555), &
					!				"p0", cBeam, iBeam, .true.)
					!=============================================================================

					call FieldtoContrastSource(cResidual, cModelParams.ContrastSourceType, cInhomContrast)
					call test_isnan(cResidual)
					call ContrastSourcetoField(cResidual, cResidual, .false.)

					if (ResidualNeumann==1) then

						RandomNumber=-1
						call SubtractStoredToFieldBis(cResidualNeumann1,cResidual,cBeam,RandomNumber);
						! this function does A=B-d*C
						call AbsoluteValueCalculation(TempResidualOutput2,cResidualNeumann1)
						Write (iProcIDString, '(i10.10)') iProcID
						err_unit_number = 10
						error_file = trim(cSourceParams.residualdirectory)//'errorNTest'//iProcIDString //'.dat'
						Open (Unit=err_unit_number, File=trim(error_file), Status='replace')
						Close (Unit=err_unit_number)
						Call print_error (TempResidualOutput2, err_unit_number, error_file)

						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNeumann1,cResidual,cBeam,RandomNumber);
						! this function does A=B-d*C
						call AbsoluteValueCalculation(TempResidualOutput2,cResidualNeumann1)
						err_unit_number = 10
						error_file = trim(cSourceParams.residualdirectory)//'errorN'//iProcIDString //'.dat'
						Open (Unit=err_unit_number, File=trim(error_file), Status='replace')
						Close (Unit=err_unit_number)
						Call print_error (TempResidualOutput2, err_unit_number, error_file)

					end if


					call FieldtoContrastSourceLin(cResidual, cBeam, 6, cInhomContrast)
					call test_isnan(cResidual)
					call ContrastSourcetoField(cResidual, cResidual, .false.)

					if (ResidualNeumann==1) then

						call AbsoluteValueCalculation(TempResidualOutput1,cResidual)


						Write (iProcIDString, '(i10.10)') iProcID
						err_unit_number = 10
						error_file = trim(cSourceParams.residualdirectory)//'errorN'//iProcIDString //'.dat'
						Call print_error (TempResidualOutput1, err_unit_number, error_file)

					end if

					call Copydata(cBeamTot,cResidual)

					! TEST EXPORT DATA
					!=============================================================================
					!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(3333), &
					!				"p0", cResidual, iBeam, .true.)
					!=============================================================================

					! TEST EXPORT DATA
					!=============================================================================
					!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(1111), &
					!				"p0", cBeamST, iBeam, .true.)
					!=============================================================================

					! TEST EXPORT DATA
					!=============================================================================
					!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(4444), &
					!				"p0", cResidual, iBeam, .true.)
					!=============================================================================

					AlphaCG=1
					OmegaCG=1
					RhoOld=1
					RhoNew=0

					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cPCG);
					else
						call PlanewaveField(cPCG);
					end if

					call test_isnan(cPCG)


					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cP);
					else
						call PlanewaveField(cP);
					end if

					call test_isnan(cP)


					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cXCG);
					else
						call PlanewaveField(cXCG);
					end if
					call test_isnan(cXCG)


					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cTemp);
					else
						call PlanewaveField(cTemp);
					end if
					call test_isnan(cTemp)


					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cV);
					else
						call PlanewaveField(cV);
					end if
					call test_isnan(cV)


					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cT);
					else
						call PlanewaveField(cT);
					end if
					call test_isnan(cT)

				end if

			    ! --------------------------------------

			    call LoadField(cBeamST, "BeamSol0") !L.D. 12-10-2010
			    call LoadField(cPCG, "BeamStatic")
			    call LoadField(cXCG, "BeamStatic")
			    call LoadField(cP, "BeamStatic")
			    call LoadField(cTemp, "BeamStatic")
			    call LoadField(cV, "BeamStatic")
			    call LoadField(cT, "BeamStatic")

			    ! Output slices if required
			    if((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
			    .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
			    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(0_i8b), &
			    "p0", cBeamST, iBeam, .true.)
			    end if

			    call LoadField(cBeamST, "BeamStatic") !L.D. 12-10-2010

		end if

		!------------------------------------------
		! end of the Pre-loop Phase
		!------------------------------------------
		


        !--------------------------------------------------------------------------------------
        !Needed only when the domain is splitted in more than a single computational domain

        !!!
        !!!		if (iBeam>0 .and. cModelParams%Numiterations(1+iBeam)>1) then
        !!!			!Apply the previous beam solution to the current one
        !!!			call PrintToLog("************************************************************",0)
        !!!			write(acTemp,'("****** Calculate Beam-to-Beam correction for Beam ",I3," ******")') iBeam
        !!!			call PrintToLog(acTemp,0)
        !!!			call PrintToLog("************************************************************",0)
        !!!
        !!!			!First deallocate the linear field solution to clear memory
        !!!			call GridDistr0Deallocate(cBeam%cGrid)
        !!!
        !!!			! Obtain the new starting positions for the previous beam
        !!!			iStartT = cRefBeam.iStartT + iBeamOffsetT(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!			iStartX = cRefBeam.iStartX + iBeamOffsetX(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!			iStartY = cRefBeam.iStartY + iBeamOffsetY(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!			iStartZ = cRefBeam.iStartZ + (iBeam-1)*cRefBeam%iDimZ
        !!!
        !!!		    ! Initialize the contrast source beam on the previous beam position
        !!!			call InitSpace(cBeamCS, iSI_FIELD, cRefBeam%bYSymm, &
        !!!							cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
        !!!							iStartT, iStartX, iStartY, iStartZ,  &
        !!!							0_i8b, 0_i8b, 0_i8b, &
        !!!							cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
        !!!			call InitGrid(cBeamCS, iProcN, iProcID, (/ .false., .false., .false., .false./));                                ! Define the structure to store the data on the near-field
        !!!			call GridDistr0CreateEmpty(cBeamCS%cGrid)
        !!!
        !!!			call LoadField(cBeamCS, "PrevBeam")
        !!!
        !!!			if (iBeam>1) then
        !!!				! First z-plane in cBeamCS is set to zero; this plane is already
        !!!				! included in the CSPlane calculation
        !!!				call ZeroFirstPlane(cBeamCS)
        !!!			end if
        !!!
        !!!			call FieldtoContrastSource(cBeamCS, cModelParams.ContrastSourceType, cInhomContrast);
        !!!
        !!!            ! Test whether there is a NaN in cBeamCS - indication of an error,
        !!!            ! no use to continue the program
        !!!            call test_isnan(cBeamCS)
        !!!
        !!!			! Calculate the field correction
        !!!			!Memory-efficient method, using out-place/in-place evaluation of the spaces
        !!!			call ContrastSourcetoField(cBeamCS, cBeam, .true.);
        !!!
        !!!            ! Test whether there is a NaN in cBeamCS - indication of an error,
        !!!            ! no use to continue the program
        !!!            call test_isnan(cBeam)
        !!!
        !!!			! Destroy the contrast source beam
        !!!			!! First de-reference pacD2 pointer in cBeamCS - don't deallocate pacD0 grid, pacD0 pointer
        !!!			! in cBeam now points to memory that was originally in cBeamCS
        !!!			nullify(cBeamCS%cGrid%pacD2)
        !!!			call DestructSpace(cBeamCS)
        !!!
        !!!			! Store the previous beam solution in a special temp file for use in the iterative scheme
        !!!			call StoreField(cBeam,"PrevBeamContributions");
        !!!
        !!!			if (iBeam>1) then
        !!!				!Apply the correction plane solution to the current one
        !!!				call PrintToLog("***************************************************************",0)
        !!!				write(acTemp,'("****** Calculate CFPlane-to-Beam correction for Beam ",I3," ******")') iBeam
        !!!				call PrintToLog(acTemp,0)
        !!!				call PrintToLog("***************************************************************",0)
        !!!
        !!!				call GridDistr0Deallocate(cBeam%cGrid)
        !!!
        !!!				! Obtain the new starting positions for the previous beam
        !!!				iStartT = cRefBeam.iStartT + iBeamOffsetT(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!				iStartX = cRefBeam.iStartX + iBeamOffsetX(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!				iStartY = cRefBeam.iStartY + iBeamOffsetY(cRefBeam, (iBeam-1)*cRefBeam%iDimZ)
        !!!				iStartZ = cRefBeam.iStartZ + (iBeam-1)*cRefBeam%iDimZ
        !!!
        !!!				! Initialize the contrast source beam on the 2nd previous beam position
        !!!                ! Although we only have a single plane in z, i.e. iDimZ=1, we need to create a complete
        !!!                ! space since it was obtained and stored from an entire space. In Distr. 0, each
        !!!                ! processor contains a limited number of xyz positions, and to get the right ones in the
        !!!                ! right places we need a space with the same dimensions as from which the plane was
        !!!                ! obtained. Later, when the data is in Distr. 2, we will define a space with iDimZ=1.
        !!!				call InitSpace(cPlaneCF, iSI_FIELDSLICE, cRefBeam%bYSymm, &
        !!!								cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
        !!!								iStartT, iStartX, iStartY, iStartZ,  &
        !!!								0_i8b, 0_i8b, 0_i8b, &
        !!!								cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
        !!!				call InitGrid(cPlaneCF, iProcN, iProcID, (/ .false., .false., .false., .false./));                                ! Define the structure to store the data on the near-field
        !!!
        !!!				call GridDistr0CreateEmpty(cPlaneCF%cGrid)
        !!!				call LoadPlane(cPlaneCF, "PrevPlane")
        !!!
        !!!				! Calculate the field correction
        !!!				call CFPlaneSourcetoField(cPlaneCF, cBeam);
        !!!
        !!!                ! Test whether there is a NaN in cBeam - indication of an error,
        !!!                ! no use to continue the program
        !!!                call test_isnan(cBeam)
        !!!
        !!!				! Destroy the contrast source slice
        !!!				call DestructSpace(cPlaneCF)
        !!!
        !!!				! Add the CF-plane solution to the correction
        !!!				! of the previous beam solution and store the result
        !!!				call AddStoredToField(cBeam,"PrevBeamContributions");
        !!!				call StoreField(cBeam,"PrevBeamContributions");
        !!!
        !!!			end if
        !!!
        !!!			! Add the Primary field solution to the previous beam contributions
        !!!			call AddStoredToField(cBeam,"BeamSol0");
        !!!
        !!!		end if


        !-------------------------------------------------------------------------------------
        ! start Loop Phase for all different Schemes
        !-------------------------------------------------------------------------------------


        !If numiterations>1 then calculate the iterative process for a given number of iterations
        if (cModelParams%numiterations(1+iBeam) > 1) then

            !-----------------------------------------------------------
            ! Steepest Descent case
            !-----------------------------------------------------------
            if (CG==1) then
            	
            	do iIter=1, cModelParams%Numiterations(1+iBeam)-1

					call PrintToLog("***************************************************************",0)
					write(acTemp,'("****** Steepest Descent Scheme calculate Iterative Beam estimate ",I3," for Beam ",I3," *****")') iIter,iBeam
					call PrintToLog(acTemp,0)
					call PrintToLog("***************************************************************",0)


					call sumandmultiplyCG(cBeamF,AlphaCoeffCGscheme,cDirection)
					! this function does A=A+b*(C)

					RandomNumber=0
					call SubtractStoredToFieldBis(cTemp,cBeamF,cBeamF,RandomNumber);
					! this function does A=B-d*C

					!------------------------------------------------------------------------------------
					!----- Output -----------------------------------------------------------------------
					if((    ((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
						.or.(cModelParams%Slicesavespecifier==iAI_LAST))&
						.and.iIter==cModelParams%Numiterations(1+iBeam)-1)&
						.or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
					call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(iIter), &
						"p", cBeamF, iBeam, .true.)
					end if
					!------------------------------------------------------------------------------------
					!------------------------------------------------------------------------------------


					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

					call PrintToLog("*************************************************************",0)
					call PrintToLog("****** Calculating Residual n ********************************",0)
					call PrintToLog("*************************************************************",0)
					! Calculate residual0
					if (NKC==1) then
						RandomNumber=0;
						call SubtractStoredToFieldBis(cResidualNKCDX,cTemp,cResidualNKCDX,RandomNumber);
						! this function does A=B-d*C
						RandomNumber=0;
						call SubtractStoredToFieldBis(cResidualNKCDY,cTemp,cResidualNKCDY,RandomNumber);
						! this function does A=B-d*C
						RandomNumber=0;
						call SubtractStoredToFieldBis(cResidualNKCDZ,cTemp,cResidualNKCDZ,RandomNumber);
						! this function does A=B-d*C
					end if

					call FieldtoContrastSource(cTemp, cModelParams.ContrastSourceType, cInhomContrast)

					call test_isnan(cTemp)

					call ContrastSourcetoField(cTemp, cTemp, .false.)


					call PrintToLog("****** Introducing the contrast in kappa ********************************",0)

					if (NKC==1) then
						call FieldtoContrastSource(cResidualNKCDX, 10, cInhomContrast)
						! cResidualNKCDX contains dx((dxkappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
						! cResidualNKCDX contains G conv (dx((dxkappa)P0^2)))

						call FieldtoContrastSource(cResidualNKCDY, 11, cInhomContrast)
						! cResidualNKCDY contains dy((dykappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
						! cResidualNKCDY contains G conv (dy((dykappa)P0^2))

						call FieldtoContrastSource(cResidualNKCDZ, 12, cInhomContrast)
						! cResidualNKCDZ contains dz((dzkappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)
						! cResidualNKCDZ contains G conv (dz((dzkappa)P0^2))

						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
						! this function does A=B-d*C
						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
						! this function does A=B-d*C

						RandomNumber=-0.5;
						call SubtractStoredToFieldBis(cResidual,cResidual,cResidualNKCDX,RandomNumber);
						! this function does A=B-d*C
						!The standard Residual schould be summed with it

					end if

					call sumandmultiplyCGbis(cResidual,cBeam,cBeamF,cTemp)
					! this function does A=B-C+D

					call MPI_BARRIER(MPI_COMM_WORLD, iErr);                                           !L.D. 07-11-2011
					call AbsoluteValueCalculation(TempResidualOutput1,cResidual)                      !L.D. 07-11-2011
					call AbsoluteValueCalculation(TempResidualOutput2,cBeam)                          !L.D. 07-11-2011

					Call print_error (TempResidualOutput1, err_unit_number, error_file)
					call MPI_BARRIER(MPI_COMM_WORLD, iErr);                                           !L.D. 07-11-2011

					call DestructSpace(cTemp3)

					call InitSpace(cTemp3, iSI_FIELD, cRefBeam%bYSymm, &
						cRefBeam%iDimT, cRefBeam%iDimX, cRefBeam%iDimY, cRefBeam%iDimZ,  &
						iStartT, iStartX, iStartY, iStartZ,  &
						0_i8b, 0_i8b, 0_i8b, &
						cRefBeam%dFnyq, cRefBeam%dTanX, cRefBeam%dTanY, cRefBeam%dTanT )
					call InitGrid(cTemp3, iProcN, iProcID, (/ .false., .false., .false., .false./));

					if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
						call PrimarySourcetoField(cTemp3);
					else
						call PlanewaveField(cTemp3);
					end if

					call test_isnan(cTemp3)

					RandomNumber=0
					call SubtractStoredToFieldBis(cTemp3,cResidual,cResidual,RandomNumber);
					! this function does A=B-d*C
					cTemp3%iSpaceIdentifier=iSI_FIELDSLICE


					call PrintToLog("**************************************************************",0)
					call PrintToLog("****** Calculating Direction n ********************************",0)
					call PrintToLog("**************************************************************",0)
					! Calculate direction0

					call ContrastSourcetoFieldConj(cTemp3, cDirection, .false.) !! put the complex conjugate G
					! convolves the resudual with the green's function

					RandomNumber=0
					call SubtractStoredToFieldBis(cTemp4,cDirection,cDirection,RandomNumber)        !!L.D. 29-08-2011
					! this function does A=B-d*C

					if (NKC==1) then
						!Now CResidualNKCDX-DY-DZ=G^H(r(pn))
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDX,cDirection,cDirection,RandomNumber)        !!L.D. 29-08-2011
						! this function does A=B-d*C
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDY,cDirection,cDirection,RandomNumber)        !!L.D. 29-08-2011
						! this function does A=B-d*C
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDZ,cDirection,cDirection,RandomNumber)        !!L.D. 29-08-2011
						! this function does A=B-d*C

					end if

					call FieldtoContrastSourceLin(cDirection, cDirection, 7, cInhomContrast)
					call FieldtoContrastSourceLin(cTemp4,cTemp4,8,cInhomContrast)                   !!L.D. 29-08-2011

					cDirection%iSpaceIdentifier=iSI_FIELD
					cTemp4%iSpaceIdentifier=iSI_FIELD                                               !!L.D. 29-08-2011

					call sumandmultiplyCGtris(cDirection,cBeamF,cResidual)
					! this function does A=A*B-C

					RandomNumber=-1
					call SubtractStoredToFieldBis(cDirection,cDirection,cTemp4,RandomNumber);       !!L.D. 29-08-2011
					! this function does A=B-d*C

					if (NKC==1) then

						! cResidualNKCDX-DY-DZ contain G^H(r(pn)), cBeamF contains pn
						call FieldtoContrastSourceLin(cResidualNKCDX,cBeamF, 16, cInhomContrast)
						cResidualNKCDX%iSpaceIdentifier=iSI_FIELD
						call FieldtoContrastSourceLin(cResidualNKCDY,cBeamF, 17, cInhomContrast)
						cResidualNKCDY%iSpaceIdentifier=iSI_FIELD
						call FieldtoContrastSourceLin(cResidualNKCDZ,cBeamF, 18, cInhomContrast)
						cResidualNKCDZ%iSpaceIdentifier=iSI_FIELD

						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
						!this function does A=B-d*C
						RandomNumber=1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
						!this function does A=B-d*C

						RandomNumber=1;
						call SubtractStoredToFieldBis(cDirection,cDirection,cResidualNKCDX,RandomNumber);
						!this function does A=B-d*C

					end if

					RandomNumber=0
					call SubtractStoredToFieldBis(cTemp2,cDirection,cDirection,RandomNumber);
					! this function does A=B-d*C

					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Calculating Theta1 n ********************************",0)
					call PrintToLog("**********************************************************",0)

					!cTheta1 needs to be equal to cBeamF
					RandomNumber=0
					call SubtractStoredToFieldBis(cTheta1,cBeamF,cBeamF,RandomNumber);
					! this function does A=B-d*C

					! Calculate Theta1
					call FieldtoContrastSourceLin(cDirection, cTheta1, 6, cInhomContrast)           !!L.D. 29-08-2011

					if (NKC==1) then
						!cResidualNKCDX needs to be equal to cDirection AND cTemp2=cDirection, cTheta1=pn here
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDX,cTemp2,cTheta1,RandomNumber)
						! this function does A=B-d*C
						!cResidualNKCDY needs to be equal to cTheta1=p0 here
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDY,cTemp2,cTheta1,RandomNumber)
						! this function does A=B-d*C
						!cResidualNKCDZ needs to be equal to cTheta1=p0 here
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDZ,cTemp2,cTheta1,RandomNumber)
						! this function does A=B-d*C

						call FieldtoContrastSourceLin(cResidualNKCDX, cTheta1, 13, cInhomContrast)
						call FieldtoContrastSourceLin(cResidualNKCDY, cTheta1, 14, cInhomContrast)
						call FieldtoContrastSourceLin(cResidualNKCDZ, cTheta1, 15, cInhomContrast)

						call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
						call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
						call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)

						RandomNumber=-1
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber)
						! this function does A=B-d*C
						RandomNumber=-1
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber)

					end if



					RandomNumber=0
					call SubtractStoredToFieldBis(cTheta1,cTemp2,cTemp2,RandomNumber);              !!L.D. 29-08-2011
					! this function does A=B-d*C

					call ContrastSourcetoField(cDirection, cDirection, .false.)                     !!L.D. 29-08-2011

					RandomNumber=1
					call SubtractStoredToFieldBis(cTheta1,cDirection,cTheta1,RandomNumber);         !!L.D. 29-08-2011
					! this function does A=B-d*C

					if (NKC==1) then
						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cTheta1,cResidualNKCDX,RandomNumber);
						! this function does A=B-d*C
					end if

					RandomNumber=0
					call SubtractStoredToFieldBis(cDirection,cTemp2,cTemp2,RandomNumber);           !!L.D. 29-08-2011
					! this function does A=B-d*C


					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Calculating Theta2 n *******************************",0)
					call PrintToLog("**********************************************************",0)

					!cTheta2 needs to be equal to cDirection
					RandomNumber=0
					call SubtractStoredToFieldBis(cTheta2,cDirection,cDirection,RandomNumber);
					! this function does A=B-d*C

					if (NKC==1) then
						!cResidualNKCDX needs to be equal to cDirection
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDX,cDirection,cDirection,RandomNumber);
						! this function does A=B-d*C
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDY,cDirection,cDirection,RandomNumber);
						RandomNumber=0
						call SubtractStoredToFieldBis(cResidualNKCDZ,cDirection,cDirection,RandomNumber);
					end if

					! Calculate Theta2
					call FieldtoContrastSourceLin(cTheta2, cTheta2, 9, cInhomContrast)              !!L.D. 29-08-2011
					call test_isnan(cTheta2)
					call ContrastSourcetoField(cTheta2, cTheta2, .false.)

					if (NKC==1) then

						call FieldtoContrastSource(cResidualNKCDX, 10, cInhomContrast)
						! cResidualNKCDX contains dx((dxkappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDX, cResidualNKCDX, .false.)
						! cResidualNKCDX contains G conv (dx((dxkappa)P0^2)))

						call FieldtoContrastSource(cResidualNKCDY, 11, cInhomContrast)
						! cResidualNKCDY contains dy((dykappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDY, cResidualNKCDY, .false.)
						! cResidualNKCDY contains G conv (dy((dykappa)P0^2))

						call FieldtoContrastSource(cResidualNKCDZ, 12, cInhomContrast)
						! cResidualNKCDZ contains dz((dzkappa)P0^2)
						call ContrastSourcetoField(cResidualNKCDZ, cResidualNKCDZ, .false.)
						! cResidualNKCDZ contains G conv (dz((dzkappa)P0^2))

						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDY,RandomNumber);
						! this function does A=B-d*C
						RandomNumber=-1;
						call SubtractStoredToFieldBis(cResidualNKCDX,cResidualNKCDX,cResidualNKCDZ,RandomNumber);
						! this function does A=B-d*C

						RandomNumber=-0.5;
						call SubtractStoredToFieldBis(cResidual,cResidual,cResidualNKCDX,RandomNumber);
						! this function does A=B-d*C
						!The standard Residual schould be summed with it

					end if


					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Calculating Alpha n ********************************",0)
					call PrintToLog("**********************************************************",0)

					! Calculate Phi2
					call TwoRealValueCalculation(Phi2,cResidual,cTheta1)
					call AbsoluteValueCalculation(TempPhi,cTheta2)
					Phi2=Phi2/TempPhi
					Phi2=Phi2/4

					! Calculate Phi3
					call TwoRealValueCalculation(Phi3,cResidual,cTheta2)
					call AbsoluteValueCalculation(TempPhi,cTheta1)
					Phi3=Phi3+TempPhi
					call AbsoluteValueCalculation(TempPhi,cTheta2)
					Phi3=Phi3/TempPhi
					Phi3=Phi3/2

					! Calculate Phi4
					call TwoRealValueCalculation(Phi4,cTheta1,cTheta2)
					Phi4=Phi4/TempPhi
					Phi4=Phi4*3
					Phi4=Phi4/4

					! Calculate alphaN
					z=(((Phi3*Phi4)-(3*Phi2))/6-(Phi4**3)/27)

					q=(Phi3/3)-((Phi4**2)/9)

					acoefficient1=(z)+sqrt((q**3)+(z**2))
					acoefficient2=(z)-sqrt((q**3)+(z**2))

					if (acoefficient1>0) then
						singacoefficient1=+1;
					else
						singacoefficient1=-1;
					end if

					if (acoefficient2>0) then
						singacoefficient2=+1;
					else
						singacoefficient2=-1;
					end if

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

					acoefficient1 = (singacoefficient1)*((acoefficient1*singacoefficient1)**(1./3.))
					acoefficient2 = (singacoefficient2)*((acoefficient2*singacoefficient2)**(1./3.))


					AlphaCoeffCGscheme=acoefficient1+acoefficient2-Phi4/3;

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Phi2 ************************",0)
					call PrintToLog("**********************************************************",0)

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Phi3 ************************",0)
					call PrintToLog("**********************************************************",0)

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Phi4 ************************",0)
					call PrintToLog("**********************************************************",0)

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** acoefficient1 ************************",0)
					call PrintToLog("**********************************************************",0)

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** acoefficient2 ************************",0)
					call PrintToLog("**********************************************************",0)


					if (isNAN(AlphaCoeffCGscheme)) then

						call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
						call PrintToLog("**********************************************************",0)
						call PrintToLog("****** Aplha = NaN ************************",0)
						call PrintToLog("**********************************************************",0)

						AlphaCoeffCGscheme=-1;

					end if

					!AlphaCoeffCGscheme=acoefficient1+acoefficient2
					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Alpha before Optimizzation ************************",0)
					call PrintToLog("**********************************************************",0)

					!---------------------------------------------------------------------------------
					!-------------------------START OPTIMIZZATION PROCESS-------------------------------
					!---------------------------------------------------------------------------------
					deltaepsilonNegative=-0.01;
					deltaepsilon=0.01;

					do indicerefinedapproximation=1, 1000

						ThirOrderPolinomialValue=(AlphaCoeffCGscheme**3)+(Phi4*(AlphaCoeffCGscheme**2))+(Phi3*AlphaCoeffCGscheme)+(Phi2);
						ThirOrderPolinomialValueDelta=((AlphaCoeffCGscheme+deltaepsilon)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilon)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilon))+(Phi2);
						ThirOrderPolinomialValueDeltaNegative=((AlphaCoeffCGscheme+deltaepsilonNegative)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilonNegative)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilonNegative))+(Phi2);

						do while (((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))).AND.((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))))

							deltaepsilonNegative=deltaepsilonNegative*10;
							deltaepsilon=deltaepsilon*10;
							ThirOrderPolinomialValue=(AlphaCoeffCGscheme**3)+(Phi4*(AlphaCoeffCGscheme**2))+(Phi3*AlphaCoeffCGscheme)+(Phi2);
							ThirOrderPolinomialValueDelta=((AlphaCoeffCGscheme+deltaepsilon)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilon)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilon))+(Phi2);
							ThirOrderPolinomialValueDeltaNegative=((AlphaCoeffCGscheme+deltaepsilonNegative)**3)+(Phi4*((AlphaCoeffCGscheme+deltaepsilonNegative)**2))+(Phi3*(AlphaCoeffCGscheme+deltaepsilonNegative))+(Phi2);

						end do

						if ((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then

							AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilon/2;
							deltaepsilonNegative=-deltaepsilon/2;
							deltaepsilon=deltaepsilon/2;


						end if

						if ((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then

							AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilonNegative/2;
							deltaepsilon=deltaepsilonNegative/2;
							deltaepsilonNegative=-deltaepsilonNegative/2;


						end if


					end do

					!---------------------------------------------------------------------------------
					!-------------------------END OPTIMIZZATION PROCESS-------------------------------
					!---------------------------------------------------------------------------------

					!---------------------------------------------------------------------------------
					!--------------------------------LINEAR CASE--------------------------------------
					!---------------------------------------------------------------------------------
					if (LINEARCASEIMPLEMENTATION==1) then

						call PrintToLog("**********************************************************",0)
						call PrintToLog("****** C.G. Linear Case implementation********************",0)
						call PrintToLog("**********************************************************",0)

						Phi2=0
						Phi3=0
						call TwoRealValueCalculation(Phi2,cResidual,cTheta1)
						Phi2=Phi2/2
						Phi3=Phi2
						call AbsoluteValueCalculation(TempPhi,cTheta1)
						Phi2=Phi2/TempPhi
						AlphaCoeffCGscheme=-Phi2


						!-------------------OPTIMIZATION PROCESS------------------------------------------

						deltaepsilonNegative=-0.01;
						deltaepsilon=0.01;

						do indicerefinedapproximation=1, 1000

							ThirOrderPolinomialValue=AlphaCoeffCGscheme*TempPhi+Phi3
							ThirOrderPolinomialValueDelta=(AlphaCoeffCGscheme+deltaepsilon)*TempPhi+Phi3
							ThirOrderPolinomialValueDeltaNegative=(AlphaCoeffCGscheme+deltaepsilonNegative)*TempPhi+Phi3

							do while (((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))).AND.((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))))

								deltaepsilonNegative=deltaepsilonNegative*10;
								deltaepsilon=deltaepsilon*10;
								ThirOrderPolinomialValue=AlphaCoeffCGscheme*TempPhi+Phi3
								ThirOrderPolinomialValueDelta=(AlphaCoeffCGscheme+deltaepsilon)*TempPhi+Phi3
								ThirOrderPolinomialValueDeltaNegative=(AlphaCoeffCGscheme+deltaepsilonNegative)*TempPhi+Phi3

							end do

							if ((ThirOrderPolinomialValueDeltaNegative/abs(ThirOrderPolinomialValueDeltaNegative))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then

								AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilon/2;
								deltaepsilonNegative=-deltaepsilon/2;
								deltaepsilon=deltaepsilon/2;

							end if

							if ((ThirOrderPolinomialValueDelta/abs(ThirOrderPolinomialValueDelta))==(ThirOrderPolinomialValue/abs(ThirOrderPolinomialValue))) then



								AlphaCoeffCGscheme=AlphaCoeffCGscheme+deltaepsilonNegative/2;
								deltaepsilon=deltaepsilonNegative/2;
								deltaepsilonNegative=-deltaepsilonNegative/2;



							end if


						end do


					end if
					!---------------------------------------------------------------------------------
					!------------------------------END LINEAR CASE------------------------------------
					!---------------------------------------------------------------------------------

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

					call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
					call PrintToLog("**********************************************************",0)
					call PrintToLog("****** Aplha after Optimizzation ************************",0)
					call PrintToLog("**********************************************************",0)

				end do

            !-----------------------------------------------------------
            ! Bi-CGSTAB case
            !-----------------------------------------------------------
            else if (BiCGSTAB==1) then
            	do iIter=1, cModelParams%Numiterations(1+iBeam)-1

				call PrintToLog("***************************************************************",0)
				write(acTemp,'("****** Bi-CGSTAB Scheme calculate Iterative Beam estimate ",I3," for Beam ",I3," *****")') iIter,iBeam
				call PrintToLog(acTemp,0)
				call PrintToLog("***************************************************************",0)

				call sumsquare(RhoNew,cBeamTot,cResidual)
				! this function sumsquare does a=sum(sum(B.*C))
				call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in

				BetaCG=(RhoNew*AlphaCG)/(RhoOld*OmegaCG)

				call sumandmultiply(cPCG,cResidual,cV,BetaCG,OmegaCG,iIter)
				! this function does A=B+d*(A-e*C), f is the iteration number

				! TEST EXPORT DATA
				!=============================================================================
				!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(770+iIter), &
				!				"p0", cPCG, iBeam, .true.)
				!=============================================================================

				call Copydata(cTemp,cPCG)
				! this function copy B in A so A=B.

				! TEST EXPORT DATA
				!=============================================================================
				!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(880+iIter), &
				!				"p0", cTemp, iBeam, .true.)
				!=============================================================================

				call LoadField(cBeam, "BeamSol0")

				call FieldtoContrastSourceLin(cPCG, cBeam, 6, cInhomContrast)
				call test_isnan(cPCG)

				call ContrastSourcetoField(cPCG, cPCG, .false.)

				! TEST EXPORT DATA
				!=============================================================================
				!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(880+iIter), &
				!				"p0", cPCG, iBeam, .true.)
				!=============================================================================

				RandomNumber=1
				call SubtractStoredToFieldBis(cV,cTemp,cPCG,RandomNumber);
				! this function does A=B-d*C

				! TEST EXPORT DATA
				!=============================================================================
				!    call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(990+iIter), &
				!				"p0", cV, iBeam, .true.)
				!=============================================================================

				call Copydata(cPCG,cTemp)
				! Now cPCG is copied again in cPCG structure

				call sumsquare(AlphaCG,cBeamTot,cV)
				! this function sumsquare does a=sum(sum(B.*C))

				AlphaCG=RhoNew/AlphaCG

				call SubtractStoredToFieldBis(cResidual,cResidual,cV,AlphaCG);
				! this function sumsquare does A=B-d*C

				call Copydata(cTemp,cResidual)
				! cTemp is now cResidual

				call LoadField(cBeam, "BeamSol0")

				call FieldtoContrastSourceLin(cResidual, cBeam, 6, cInhomContrast)

				call test_isnan(cResidual)

				call ContrastSourcetoField(cResidual, cResidual, .false.)


				RandomNumber=1
				call SubtractStoredToFieldBis(cT,cTemp,cResidual,RandomNumber);
				! this function does A=B-d*C

				call Copydata(cResidual,cTemp)
				! cResidual contains now cResidual again

				call sumsquare(OmegaCG,cT,cResidual)
				! this function sumsquare does a=sum(sum(B.*C))

				call sumsquare(RandomNumber,cT,cT)
				! this function sumsquare does a=sum(sum(B.*C))


				OmegaCG=OmegaCG/RandomNumber

				call LoadField(cBeam, "BeamStatic")

				if (iIter==1)then
					call Copydata(cXCG,cBeam)
				end if

				call SubtractStoredToFieldTris(cXCG,cPCG,cResidual,AlphaCG,OmegaCG);
				! this function does A=A+d*B+e*C

				! TEST EXPORT DATA
				!=============================================================================
				!call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(950+iIter), &
				!			"p0", cXCG, iBeam, .true.)
				!=============================================================================

				call SubtractStoredToFieldBis(cResidual,cResidual,cT,OmegaCG);
				! this function does A=B-d*C

				if (ResidualNeumann==1) then

					call AbsoluteValueCalculation(TempResidualOutput1,cResidual)


					Write (iProcIDString, '(i10.10)') iProcID
					err_unit_number = 10
					error_file = trim(cSourceParams.residualdirectory)//'errorN'//iProcIDString //'.dat'
					Call print_error (TempResidualOutput1, err_unit_number, error_file)

				end if

				RhoOld=RhoNew

				call Copydata(cP,cXCG)

				!-------------------------------------------------------------------------
				if (ResidualNeumann==1) then

					call LoadField(cResidualNeumann1, "BeamSol0")

					RandomNumber=0
					call SubtractStoredToFieldBis(cResidualNeumann2,cXCG,cResidualNeumann2,RandomNumber);
					! this function does A=B-d*C
					call FieldtoContrastSourceLin(cResidualNeumann2, cResidualNeumann1, 6, cInhomContrast)
					call test_isnan(cResidual)
					call ContrastSourcetoField(cResidualNeumann2, cResidualNeumann2, .false.)

					RandomNumber=1
					call SubtractStoredToFieldBis(cResidualNeumann2,cResidualNeumann2,cXCG,RandomNumber);
					! this function does A=B-d*C

					call LoadField(cResidualNeumann1, "BeamStatic")

					RandomNumber=-1
					call SubtractStoredToFieldBis(cResidualNeumann2,cResidualNeumann1,cResidualNeumann2,RandomNumber);
					! this function does A=B-d*C
					call AbsoluteValueCalculation(TempResidualOutput1,cResidualNeumann2)
					Write (iProcIDString, '(i10.10)') iProcID
					err_unit_number = 10
					error_file = trim(cSourceParams.residualdirectory)//'errorNTest'//iProcIDString //'.dat'

					Call print_error (TempResidualOutput1, err_unit_number, error_file)

				end if
				!-------------------------------------------------------------------------

				call AddStoredtoField(cP,"BeamSol0");
				! Output slices if required
				if((    ((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
					.or.(cModelParams%Slicesavespecifier==iAI_LAST))&
					.and.iIter==cModelParams%Numiterations(1+iBeam)-1)&
					.or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
				call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(iIter), &
					"p", cP, iBeam, .true.)
				end if

    		end do

        !-----------------------------------------------------------
        ! start Neumann case
        !-----------------------------------------------------------     
        else      
        	do iIter=1, cModelParams%Numiterations(1+iBeam)-1
			
				call PrintToLog("***************************************************************",0)
				write(acTemp,'("****** Neumann Scheme Calculate Iterative Beam estimate ",I3," for Beam ",I3," *****")') iIter,iBeam
				call PrintToLog(acTemp,0)
				call PrintToLog("***************************************************************",0)
				cModelParams%iIter = iIter
				!-----------------------------------------------------------
				!LCSM - Linearized Neumann
				!-----------------------------------------------------------
				if (LCSM==1) then

				    if (iIter==1) then
				        call LoadField(cBeam, "BeamStatic")
				    else
				        call LoadField(cBeam, "BeamUpgrande")
				    end if

				    call LoadField(cBeamST, "BeamSol0")

				    call FieldtoContrastSourceLin(cBeam, cBeamST, 6, cInhomContrast)

				else
				    call FieldtoContrastSource(cBeam, cModelParams.ContrastSourceType ,cInhomContrast)
				end if
				
				! Test whether there is a NaN in cBeam - indication of an error,
				! no use to continue the program
				call test_isnan(cBeam)
				
				call ContrastSourcetoField(cBeam, cBeam, .false.)
				!! Residual Calculation
				if (ResidualNeumann==1) then

				    call LoadField(cResidualNeumann1, "BeamSol0")

				    if (iIter==1) then

				        call AbsoluteValueCalculation(TempResidualOutput2,cResidualNeumann1)

				        call AbsoluteValueCalculation(TempResidualOutput1,cBeam)

				        Write (iProcIDString, '(i10.10)') iProcID
				        err_unit_number = 10
				        error_file = trim(cSourceParams.residualdirectory)//'errorN'//iProcIDString //'.dat'
				        Open (Unit=err_unit_number, File=trim(error_file), Status='replace')
				        Close (Unit=err_unit_number)

				        Call print_error (TempResidualOutput2, err_unit_number, error_file)
				        Call print_error (TempResidualOutput1, err_unit_number, error_file)

				    end if
				    Write (iProcIDString, '(i10.10)') iProcID

				    RandomNumber=-1
				    call SubtractStoredToFieldBis(cResidualNeumann2,cBeam,cResidualNeumann1,RandomNumber); ! this gives p at the given iteration
				    ! this function does A=B-d*C

				    RandomNumber=0
				    call SubtractStoredToFieldBis(cResidualNeumann1,cResidualNeumann2,cResidualNeumann1,RandomNumber); ! this gives p at the given iteration
				    ! this function does A=B-d*C

		            call FieldtoContrastSource(cResidualNeumann1, cModelParams.ContrastSourceType, cInhomContrast)
		            call ContrastSourcetoField(cResidualNeumann1, cResidualNeumann1, .false.)
					if (cModelParams%Numiterations(1+iBeam)==1) cResidualNeumann1 = cBeam
				    RandomNumber=1
				    call SubtractStoredToFieldBis(cResidualNeumann2,cResidualNeumann2,cResidualNeumann1,RandomNumber);
				    ! this function does A=B-d*C

				    call LoadField(cResidualNeumann1, "BeamSol0")


				    RandomNumber=1
				    call SubtractStoredToFieldBis(cResidualNeumann2,cResidualNeumann1,cResidualNeumann2,RandomNumber);
				    ! this function does A=B-d*C

				    call AbsoluteValueCalculation(TempResidualOutput1,cResidualNeumann2)

				    Call print_error (TempResidualOutput1, err_unit_number, error_file)

				end if
 
				!--------------------------------------------------------------------------------------------------------

				call UpdateRRMSerror(cRRMSNorm, cBeam, iBeam, iIter)
				! If necessary, add the Previous beam contributions to the current solution

				!!----------------------------------------------------------------- In case beam splitting is implemented
				!if(iBeam>0) then
				!	call AddStoredToField(cBeam,"PrevBeamContributions");
				!end if

				! If necessary, store the first plane of the field correction for a next beam
				!if((iIter==cModelParams%Numiterations(1+iBeam)-1).and.(iBeam>0.and.iBeam<cModelParams%NumBeams-1)) then
				!	call StorePlane(cBeam, "PrevPlane")
				!end if

				!-----------------------------------------------------------
				!LCSM - Linearized Neumann
				!-----------------------------------------------------------

				if (LCSM==1) then
				    call AddStoredtoField(cBeam,"BeamStatic");
				    call StoreField(cBeam,"BeamUpgrande")
				    call LoadField(cBeamTot, "BeamSol0")
				    call AddStoredtoField(cBeamTot,"BeamUpgrande");
				    ! Output slices if required
				    if((    ((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
				        .or.(cModelParams%Slicesavespecifier==iAI_LAST))&
				        .and.iIter==cModelParams%Numiterations(1+iBeam)-1)&
				        .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
				    	call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(iIter), &
				        "p", cBeamTot, iBeam, .true.)
				    end if
				else
					if((    ((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
				        .or.(cModelParams%Slicesavespecifier==iAI_LAST))&
				        .and.iIter==cModelParams%Numiterations(1+iBeam)-1)&
				        .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
						call Storeslices(trim(sOutputDir) // trim('ContrastSrc') // int2str(iIter), &
				                    "p", cBeam, iBeam, .true.)
				    endif
				    
				    
					if (cModelParams%PrimarySourceType == iEI_LOADFIELD) then
						call AddStoredtoField(cBeam, trim(cModelParams.FieldFilename)); 
					else
						call AddStoredtoField(cBeam, "BeamSol0"); 
					end if
				    ! Output slices if required
				    if((    ((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) &
				        .or.(cModelParams%Slicesavespecifier==iAI_LAST))&
				        .and. iIter==cModelParams%Numiterations(1+iBeam)-1)&
				        .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
				    	call Storeslices(trim(sOutputDir) // trim(sOutputRoot) // int2str(iIter), &
				        "p", cBeam, iBeam, .true.)
				    end if

				end if
			

				! Test whether there is a NaN in cBeam - indication of an error,
				! no use to continue the program
				call test_isnan(cBeam)
				! if (iIter == cModelParams%Numiterations(1+iBeam)-1) call StoreField(cBeam,"BeamSol0")
				
	    	end do
    	 end if
        !-----------------------------------------------------------
        ! end of the Loop Phase
        !-----------------------------------------------------------

            if(iBeam<cModelParams%NumBeams-1) then
                call StoreField(cBeam,"PrevBeam")
            end if

            call ExportRRMSerror(cRRMSNorm, iProcID)

        end if
        !-----------------------------------------------------------
        ! end of the program
        !-----------------------------------------------------------

        call DestructSpace(cBeam)

        if (LCSM==1) then
            call DestructSpace(cBeamST)
            call DestructSpace(cBeamTot)
        end if

        if (BiCGSTAB==1) then
            call DestructSpace(cResidual)
            call DestructSpace(cPCG)
            call DestructSpace(cV)
            call DestructSpace(cXCG)
            call DestructSpace(cTemp)
            call DestructSpace(cP)
        end if

        if (CG==1) then
            call DestructSpace(cResidual)
            call DestructSpace(cDirection)
            call DestructSpace(cBeamF)
            call DestructSpace(cTemp)
            call DestructSpace(cTemp1)
            call DestructSpace(cTemp2)
            call DestructSpace(cTemp3)
            call DestructSpace(cTemp4)
            call DestructSpace(cTheta1)
            call DestructSpace(cTheta2)
        end if

        if (NKC==1) then
            call DestructSpace(cResidualNKCDX)
            call DestructSpace(cResidualNKCDY)
            call DestructSpace(cResidualNKCDZ)
        end if

        if (ResidualNeumann==1) then
            call DestructSpace(cResidualNeumann1)
            call DestructSpace(cResidualNeumann2)
        end if

        if (DFM==2) then
            call DestructSpace(cBeamDFM)
        end if

    end do

    ! Destruct the spaces used
    call DestructSpace(cRefBeam)
    if (cModelParams.ContrastSourceType /= iCI_NONLIN .AND. cModelParams.ContrastSourceType/=iCI_SCATTERER ) then
        call DestructSpace(cInhomContrast)
    end if

    ! Export and destruct RRMS Norm
    call DestroyRRMSPreviousCorrection(cRRMSNorm)
    call DestroyRRMSNorm(cRRMSNorm)

    !Deallocate arrays in cModelParams
    deallocate(cModelParams%numiterations)
    deallocate(cModelParams%xyzslicedim,cModelParams%xyzslicepos,cModelParams%xyzslicebeam,cModelParams%xyzsliceindex)
    if (cModelParams.ContrastSourceType==iCI_SCATTERER ) then
    	deallocate(ScattererParams%ClusterSliceDim)
    endif
    
    ! Deinitialization
    call StopWatchesStop;
    write (acTemp, '("Parnac - Processor ", I3, " checking out successfully.")') iProcID; call PrintToLog(acTemp, 0);
    write ( *, '("Parnac - Processor ", I3, " checking out successfully.")') iProcID
    write ( *, '(" ")')
    call CloseLogFile;

    call MPI_FINALIZE(ierr)

    END PROGRAM Parnac_Main

