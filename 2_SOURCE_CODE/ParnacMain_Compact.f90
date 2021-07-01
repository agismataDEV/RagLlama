    PROGRAM Parnac_Main


    !*******************************************************************************
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     060612  Original code (KH-LD)
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
    USE ParnacSchemes

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

    call MPI_INIT_THREAD(REQUIRED, PROVIDED, iErr)	!Added by A.M 13/09/2020 [before MPI_INIT(iErr)]
!    call MPI_INIT(iErr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, iProcN4, iErr);
    call MPI_COMM_RANK(MPI_COMM_WORLD, iProcID4, iErr);
!    call OMP_SET_NUM_THREADS(1)				!Added by A.M 13/09/2020
    
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
    if (cModelParams.ContrastSourceType/=iCI_NONLIN .AND. cModelParams.ContrastSourceType/=iCI_BUBBLE) then

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
        call DFMInit()

        !-------------------------------------------------------------------------------------
        !LCSM L.D. 13-07-2010 - Initialization for Linearized Contrast Source Method
        !-------------------------------------------------------------------------------------
		call LCSMINit()

        !--------------------------------------------------------------------------------------
        !Bi-CGSTAB L.D. 17-09-2010 - Initialization for Bi-CGSTAB
        !-------------------------------------------------------------------------------------
        call Bi_CGSTABInit()

        !--------------------------------------------------------------------------------------
        !Steepest Descent L.D. 20-07-2011 - Initialization for Steepest Descent
        !-------------------------------------------------------------------------------------
        call CGInit()
        !--------------------------------------------------------------------------------------
        ! Calculate the primary field solution

        call PrintToLog("*************************************************************",0)
        call PrintToLog("Calculate primary field solution",0)
        call PrintToLog("*************************************************************",0)

        ! Either from source description or as plane wave, this function generates the primary field
        ! given the source type defined in the input file
        if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
            call PrimarySourcetoField(cBeam);
        else
            call PlanewaveField(cBeam);
        end if

        if (ResidualNeumann==1) then

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
                call PrimarySourcetoField(cResidualNeumann1);
            else
                call PlanewaveField(cResidualNeumann1);
            end if

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
                call PrimarySourcetoField(cResidualNeumann2);
            else
                call PlanewaveField(cResidualNeumann2);
            end if

        end if

        if (DFM == 2 .OR. DFM == 3) then

            if(cModelParams%PrimarySourceType/=iEI_PLANEWAVE) then
                call PrimarySourcetoField(cBeamDFM);
            else
                call PlanewaveField(cBeamDFM);
            end if

        end if

        !---------------------------------------------------------------------------------------------------------
        ! Start pre-loop Steepest Descent scheme - Initialization of the starting values for the different spaces
        !---------------------------------------------------------------------------------------------------------
        call CGpre()

        !-------------------------------------------------------------------------------------
        ! start pre-loop full nonlinear neumann - Linear incident field generation
        call NeumannPreLoop()

        ! Test whether there is a NaN in cBeam - indication of an error,
        ! no use to continue the program
        call test_isnan(cBeam)

        !If necessary, store the primary field solution for use in iterative scheme- note that in case
        ! Dual Frequency Modality (DFM) or linearization (LCSM) are on other fields have to be stored
        if(cModelParams%Numiterations(1+iBeam)>1) then
            call StoreField(cBeam,"BeamSol0")

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
		call LCSMPre()


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
            	
            	call CGloop()

            !-----------------------------------------------------------
            ! Bi-CGSTAB case
            !-----------------------------------------------------------
            else if (BiCGSTAB==1) then
            	call Bi_CGSTABloop()

            !-----------------------------------------------------------
            ! start Neumann case
            !-----------------------------------------------------------     
            else      
                call Neumannloop()


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
    if (cModelParams.ContrastSourceType /= iCI_NONLIN .AND. cModelParams.ContrastSourceType/=iCI_BUBBLE) then
        call DestructSpace(cInhomContrast)
    end if

    ! Export and destruct RRMS Norm
    call DestroyRRMSPreviousCorrection(cRRMSNorm)
    call DestroyRRMSNorm(cRRMSNorm)

    !Deallocate arrays in cModelParams
    deallocate(cModelParams%numiterations)
    deallocate(cModelParams%xyzslicedim,cModelParams%xyzslicepos,cModelParams%xyzslicebeam,cModelParams%xyzsliceindex)
    if (cModelParams.ContrastSourceType==iCI_BUBBLE) then
    	deallocate(BubbleParams%ClusterSliceDim)
    endif
    
    ! Deinitialization
    call StopWatchesStop;
    write (acTemp, '("Parnac - Processor ", I3, " checking out successfully.")') iProcID; call PrintToLog(acTemp, 0);
    write ( *, '("Parnac - Processor ", I3, " checking out successfully.")') iProcID
    write ( *, '(" ")')
    call CloseLogFile;

    call MPI_FINALIZE(ierr)

    END PROGRAM Parnac_Main

