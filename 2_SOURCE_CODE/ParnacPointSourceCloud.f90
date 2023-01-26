MODULE ParnacPointSourceCloud

    ! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis
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
    !   The module ParnacContrastFunctions contains subroutines used in the
    !   initialization of the contrast functions, and their evaluation as a
    !   function of a certain field estimate.
    !
    ! *****************************************************************************
    !
    !   MODULES USED
    !

    USE Types
    USE Constants
    USE ParnacGeneral
    USE ParnacDataDef
    USE ParnacParamDef
    USE ParnacIdentifierDef
    USE ParnacDerivative, ONLY: &
        iFDMinOrder, DerivLookup, &
        DerivLookupInit, DerivLookupDestroy, &
        DerivativeReal, DerivativeComplex
    USE ParnacDataRedist, ONLY: &
        ReorderDistr0ToDistr1, ReorderDistr1ToDistr0, &
        ReorderDistr1ToDistr2, ReorderDistr2ToDistr1
    USE ParnacDataSpaceInit, ONLY: &
        iBeamOffsetX, iBeamOffsetY, InitSaveSlices, SpaceInfo, &
        InitSpace, InitGrid, GridDistr0CreateEmpty, DestructSpace
    USE ParnacTransformFilter, ONLY: &
        TransformT, TransformTInv, &
        TransformT_sml, TransformTInv_sml, &
        dTaperingWindow, INTERP1DFREQ, GREENS1DFREQ, &
        FFilterSpatial1D
    USE ParnacParamInit, ONLY: ScattererInit
    USE ParnacOutput, ONLY: ExportSlice
    USE SpecialFun, ONLY: &
        dSinc, INTERP1D, INTERP3D, INTERP3D_SHEPARD, LINSPACE, FIND_CLOSEST_DIVISOR
    USE ParnacDataMap, ONLY: &
        MapVt, MapVtInv
    USE ORDERPACK
    USE RP
    USE ParnacContrastFunctions, ONLY: NonlinContrastOperator_Ali
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
    !   BubbleContrastOperator   sub   Evaluate the contrast source operator
    !                                  for nonlinear microbubbles
    ! =============================================================================
    !
CONTAINS

    SUBROUTINE BubbleContrastOperator(cSpace)

        USE ParnacParamInit, ONLY: &
            ScattererInit
        USE IFPORT
        USE MPI
        ! ===============================================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 200212
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200212  Original code (KH)
        !
        ! ***********************************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine BubbleContrastOperator computes the bubble contrast
        !   source S^B(V) = rho0 * d^2 / dt^2 V_B^2 for the given space. The
        !   data in cSpace should be in Distribution 0. This implementation does not
        !   prevent aliasing in the temporal dimension.
        !
        ! ***********************************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace              io   type(space)  The space in which the point has to be given
        !

        type(Space), target, intent(inout)  ::            cSpace
        type(PointSourceCloudInput)         ::            BubbleAsPointScatterer

        ! ***********************************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   dBubbleLocation             dp    Global coordinate of the location of the bubble
        !   dGlobXYZ                    dp    global coordinates of each bubble scatterer
        !   GlobalPos                   dp    global coordinate of the 8 neighbouring grid points
        !   dLambdaMM                   dp    Normalization factor
        !   dRadius                     dp    Distance of the requested point in space to the location
        !                                     of the bubble as a point scatterer
        !   iIndex                      i8b   loop counter over the positions in the space
        !   i, i_con,i_loc              i8b   loop counter
        !   iDimT/X/Y/Z                 dp    the discrete lengths
        !   iPstart/end                 dp    the starting and ending position of the pulse in the
        !                                     zero-padded array
        !   ProcLoop                    i8b   loop counter over the number of processors used
        !   iBubble                     i8b   loop counter over the number of microbubbles used
        !   XYZ_i1                      i8b   array containing the index from the coordinates of each
        !                                     neighbouring to the bubble point scatterer
        !   XYZIndex                    i8b   Index from the coordinates of one bubble scatterer
        !   IndexTempLoc                i8b   Temporary index containing the number of the location
        !                                     or the contrast source terms, which were read from the
        !                                     file stored
        !   arBufferXYZ                 i8b   Temporary buffer containing a space trace
        !   arBuffer1                   dp    Temporary buffer containing a time trace
        !   dBubbleLocationN            dp    Global coordinate of the location of the bubble,
        !                                     normalized with respect to the wavelength (as employed
        !                                     everywhere in the code)
        !   dBubbleContrast             dp    array containing the contrast source term from all the
        !                                     bubble scatterers
        !   temploc/con                 dp    Temporary array buffers containing the location or the
        !                                     contrast source term of some/all bubbles
        !   n_samples                   i8b   Number used for zero-padding
        !   T_zp_start/end              dp    Values for the cycles before/after the driving pulse
        !   dK_c                        dp    spatial angular cutoff frequency of the filter
        !   RealPressure/Time/R         dp    Array with the pressure/time/radius of the bubble in the
        !                                     intended grid point
        !   RealPressureN/TimeN/R_N     dp    Zero-padded Array with the pressure/time/radius of the
        !                                     bubble in the intended grid point
        !   RealPressure_i1             dp    Array with the pressure of the neighbouring grid points
        !   Interp_PressureGL           dp    Interpolated Pressure from the 8 neigbouring grid points
        !                                     in respect to the global location of the microbubble
        !   Interp_PressureBL           dp    Interpolated Pressure from the 8 neigbouring grid points
        !                                     in respect to the normalized microbbuble location
        !   V_dd_norm                   dp    Array with the calculated contrast source term
        !   exists/exists_loc      logical    variable used to check if files exist in a specific dir
        ! ***********************************************************************************************

        type(Grid), pointer                 ::          pcGrid

        integer                             ::          iErr, iStat(MPI_STATUS_SIZE), ErrCode
        real(dp)                            ::          dRadius(1), Domain_Range(2, 3)
        integer(i8b)                        ::          iDimX, iDimY, iDimZ, iDimT, iDimW, count_neighbs, NeighbGridPoints(16**3, 3)
        integer(i8b)                        ::          iBubble, iIndex, i, j, i_neighbgp, iStart, x, y, z, iSc, dTotalCloudN, CloudIndex
        integer                             ::          ibDimT, BubbleProc, i_loc
        integer                             ::          Bubble_disps(cSpace%cGrid%iProcN), Bubble_condisps(cSpace%cGrid%iProcN), Bubble_count(cSpace%cGrid%iProcN)

        character(len=1024)                 ::          actemp, MPI_FILENAME

        integer(i8b), allocatable           ::          BubbleID(:), tempbuffer(:), temploc(:, :)
        real(dp), allocatable               ::          dBubbleContrast(:), dBubbleLocationN(:, :), tempcon(:)
        integer(i8b), allocatable           ::          CloudID(:)
        real(dp), allocatable               ::          dCloudContrast(:), dCloudLocationN(:, :)

        ! Parameters for Rayleigh - Plesset Solver
        LOGICAL(lgt)                        ::          exists, exists_loc, result

        !Filtering and windowing parameters
        real(dp)                            ::          dTaperSupportWindow(cSpace%iDimT), dTaperMaxFreqWindow(cSpace%iDimT/2 + 1), dLeftBand, dRightBand

        real(dp)                            ::          dTaperX(cSpace%iDimX), dTaperY(cSpace%iDimY), dTaperZ(cSpace%iDimZ), dTaperWidth
        integer(i8b)                        ::          iIndexX, iIndexY, iIndexZ, iOmega

        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries
        !
        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   PrintToLog
        !   RANDOM_SEED
        !   RANDOM_NUMBER
        !   LINSPACE
        !   INTERP1D
        !   INTERP3D
        !   StoreVarContrast
        !   LoadVarContrast
        !   Bubble_Init
        !   RP_SOLVER
        !   dpointScattererContrast
        !   GeneratePSCloud
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   If the iteration is the first one,
        !02       Generate randomly spaced bubbles based on a minimum inbetween distance
        !03   Otherwise,
        !04       Get the bubbles' position and ID from stored file
        !05   For each bubble,
        !06       Find the Beam Index from the normalized position of the bubble,
        !07       Correlate the position of this bubble with the index in the array with the pressure field,
        !08       If the bubble is not on a grid point then
        !09           Find the indexes of the 8 neighbouring points,
        !10           Find in which processor each point is located and save the pressure and location to a file
        !10                          Save only the bubbles that have neighbouring points outside of the processor
        !11          For each bubble,
        !11                     Find in which of the processors' grid is the bubble located and continue with this one,
        !12           Create the time and the pressure arrays for this bubble,
        !13           Load the pressure fields in the 8 neibhouring points and do a trilinear interpolation
        !14           Upsample based on the n_pad variable to improve accuracy of the solver
        !15           Solve the Marmottant model based on the initial radius,
        !16           Restore pressure and time arrays to the initial size,
        !17           Calculate the contrast source term in the wave equation,
        !18           Save the results in arrays respective to the bubble counter ( location , contrast, radius)
        !19   For each of the bubbles(which are located in only one processor)
        !19   save their location, the contrast source term and the radius in a file for every processor
        !20   Load from every processor each file to acquire each bubble location, countrast source term and radius
        !21   Calculate the pressure field in each processor's grid considering the bubbles as point scatterers
        !21   Include on the bubbles that are inside or close to each processor domain
        !22   Upsample with a factor of 2 and downsample with a factor of 2 to get rid of high frequencies
        !22   And correct the value fo the scattered pressure field.
        !
        ! =============================================================================

        pcGrid => cSpace%cGrid

        !********************************************* START OF TAPPERING *********************************************************
        if (cModelParams%UseSupportTapering .EQV. .true.) then
            ! use tapering of the field at start and end to prevent wraparound leakage
            ! build tapering window, the first two periods and the last two periods are tapered
            ! in the contrast source
            dLeftBand = 2.0_dp
            dRightBand = 2.0_dp
            dTaperSupportWindow = dTaperingWindow(cSpace%iDimT, cSpace%dDt, dLeftBand, dRightBand)

            ! multiply the field with it
            iStart = 0
            do iIndex = 0, pcGrid%iD0LocN - 1
                pcGrid%parD0(iStart + 1:iStart + iDimT) = &
                    pcGrid%parD0(iStart + 1:iStart + iDimT)*dTaperSupportWindow
                iStart = iStart + pcGrid%iD0IS
            end do
        end if
        !************************************* START OF 3D INTERPOLATION FOR BUBBLES ***************************************************

        call ReorderDistr0ToDistr1(pcGrid); 
        call MapVt(cSpace)
        call ReorderDistr1ToDistr0(pcGrid); 
        ! This was added for multiple populations ( microbubbles and linear scatterers)
        dTotalCloudN = 0; 
        do iSc = 1, cSourceParams%ScCloudNo
            dTotalCloudN = dTotalCloudN + ScattererParams(iSc)%N
        end do

        ALLOCATE (CloudID(dTotalCloudN), dCloudLocationN(dTotalCloudN, 3))
        iMemAllocated = iMemAllocated  + PRODUCT(SHAPE(dCloudLocationN))*dpS + PRODUCT(SHAPE(CloudID))*i4bs

        ! This was added for multiple populations ( microbubbles and linear scatterers)
        CloudIndex = 0; 
        do iSc = 1, cSourceParams%ScCloudNo
            cSourceParams%iScCloud = iSc
            call ScattererInit()

            INQUIRE (FILE=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)), EXIST=exists_loc)
            if (.NOT. exists_loc) result = MAKEDIRQQ(trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)))

            write (acTemp, "('BubbleContrastOperator')"); call PrintToLog(acTemp, 2)
            !**************************************** INITIALIZE VALUES **********************************************
            iDimT = cSpace%iDimT   !t dimensions
            iDimX = cSpace%iDimX   !x dimensions
            iDimY = cSpace%iDimY   !y dimensions
            iDimZ = cSpace%iDimZ   !z dimensions

            iDimW = iDimT/2 + 1

            !Allocate for memory management
            ALLOCATE (BubbleID(ScattererParams(cSourceParams%iScCloud)%N), dBubbleLocationN(ScattererParams(cSourceParams%iScCloud)%N, 3))

            iMemAllocated = iMemAllocated + PRODUCT(SHAPE(dBubbleLocationN))*dpS + (PRODUCT(SHAPE(BubbleID)))*i4bS

            !Initialize for this variables , close to 0, double precission
            dBubbleLocationN = 1.0D-30
            BubbleID = 0
            i_loc = 0
            !**************************************** CREATE BUBBLE CLUSTER / LOAD BUBBLE PARAMETERS **********************************************

            !!This is only for iBeam =1!! I should change this
            INQUIRE (FILE=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//'Bubble_Location'//int2str(cSourceParams%iScCloud)//int2str(1)), EXIST=exists_loc)
            ! If the file for the first iteration exists , then just load the file from the location
            ! See later for more information about the process of loading
            if (.NOT. exists_loc) then
                write (acTemp, "('Randomly Position the microbubbles in the domain') "); call PrintToLog(acTemp, 2)
                BubbleAsPointScatterer%N = ScattererParams(cSourceParams%iScCloud)%N; BubbleAsPointScatterer%R0 = ScattererParams(cSourceParams%iScCloud)%R0; BubbleAsPointScatterer%ClusterDimsRatio = ScattererParams(cSourceParams%iScCloud)%ClusterDimsRatio; BubbleAsPointScatterer%MinInBetweenDist = ScattererParams(cSourceParams%iScCloud)%MinInBetweenDist
                ! Generate the cluster and save result on slices defined in input file
                call RNDGeneratePSCloud_Eff(cSpace, BubbleAsPointScatterer, dBubbleLocationN, Domain_Range)
                call SaveExtraBubbleSlices(cSpace, dBubbleLocationN, Domain_Range)
                write (acTemp, "('No. of Bubbles consisting the cluster: ', I<INT(log10(real(ScattererParams(cSourceParams%iScCloud)%N,dp)))+1>, ' in a volume of ', F5.3, ' [mL].')") ScattererParams(cSourceParams%iScCloud)%N, PRODUCT(maxval(dBubbleLocationN, 1) - minval(dBubbleLocationN, 1))*(cMediumParams%c0*1E3/cModelParams%freq0)**3*1.0D-3; call PrintToLog(acTemp, 2); 
                ibDimT = iDimT !- FLOOR(minval(dBubbleLocationN(:, 3) + cSpace%iStartZ*cSpace%dDx)/cSpace%dDt) + 20 ! This should not be changed because ibDimT is Integer
            else
                ! Load each Bubble's Position
                open (11, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('Bubble_Location'//int2str(cSourceParams%iScCloud))//int2str(cModelParams%iIter - 1)), status='OLD')
                do iBubble = 1, ScattererParams(cSourceParams%iScCloud)%N
                    read (11, *) dBubbleLocationN(iBubble, 1), dBubbleLocationN(iBubble, 2), dBubbleLocationN(iBubble, 3)
                end do
                close (11)
                ibDimT = iDimT !- FLOOR(minval(dBubbleLocationN(:, 3) + cSpace%iStartZ*cSpace%dDx)/cSpace%dDt) + 20 ! This should not be changed because ibDimT is Integer

                write (acTemp, "('No. of Bubbles consisting the cluster: ', I<INT(log10(real(ScattererParams(cSourceParams%iScCloud)%N,dp)))+1>, ' in a volume of ', F5.3, ' [mL].')") ScattererParams(cSourceParams%iScCloud)%N, PRODUCT(maxval(dBubbleLocationN, 1) - minval(dBubbleLocationN, 1))*(cMediumParams%c0*1E3/cModelParams%freq0)**3*1.0D-3; call PrintToLog(acTemp, 2); 
                ! Load i_loc which is the number of the scatterers in each cpu domain
                open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_iloc'//int2str(cSourceParams%iScCloud))//int2str(pcGrid%iProcID)), status='OLD')
                read (8, *) i_loc
                close (8)
                ALLOCATE (tempcon(i_loc*ibDimT))

                ! Load the radius of the scatterers
                if (ScattererParams(cSourceParams%iScCloud)%Distribution == 'polydisperse') then
                    open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_Radius'//int2str(cSourceParams%iScCloud))//int2str(cModelParams%iIter - 1)), status='OLD')
                    read (8, *) ScattererParams(cSourceParams%iScCloud)%R0
                    close (8)
                end if
            end if

            ! Create the global positions for every point in order to find the proper pressure from the interpolation in respect to Bubble Location
            ! Then use message passing interface in order for the processors to communicate with each other the information for the scatterers
            ! that are positioned in corner cases (close to each cpu domain boundaries)
            ! This is done to reduce the memory usage

            call ProcBubbleCount(cSpace, dBubbleLocationN, tempcon, i_loc, ibDimT)
            if (.NOT. ALLOCATED(tempcon)) ALLOCATE (tempcon(i_loc*ibDimT)); if (.NOT. exists_loc) tempcon = 0.0D0

            if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'around') ALLOCATE (temploc(size(NeighbGridPoints, 1)*ScattererParams(cSourceParams%iScCloud)%N, 2))
            iMemAllocated = iMemAllocated + PRODUCT(SHAPE(tempcon))*dpS + PRODUCT(SHAPE(temploc))*i4bs

            call TrilinInterp_SNDRCV(cSpace, dBubbleLocationN, i_loc, tempcon, temploc, BubbleID, count_neighbs)

            !************************************* COMMUNICATE PRESSURE AND LOCATION OF ALL THE SCATTERERS TO EVERY CPU **************************************************
            ! Create Bubble_count which contains the number of the bubbles inside each processors domain
            ! This is done to enhance speed . Other way is to save and load but it is totally inefficient
            Bubble_count = 0; 
            call MPI_Allgather(i_loc, 1, MPI_INTEGER, Bubble_count, 1, MPI_INTEGER, MPI_COMM_WORLD, iErr)
            if (sum(Bubble_count) .NE. ScattererParams(cSourceParams%iScCloud)%N) then
                write (acTemp, "('Error! Less bubbles that total are considered .')")
                call PrintToLog(acTemp, 1)
                call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
            end if
            call MPI_BARRIER(MPI_COMM_WORLD, iErr)

            write (acTemp, "('Communicate Bubble Location and Pressure')"); call PrintToLog(acTemp, 2)
            Bubble_disps = 0; Bubble_condisps = 0
            do i = 1, pcGrid%iProcN
                if (Bubble_count(i) > 0) then  ! Displacement variable in order to gather the value from every process
                    Bubble_disps(i) = sum(Bubble_count(1:i)) - Bubble_count(i)
                    Bubble_condisps(i) = sum(ibDimT*Bubble_count(1:i)) - ibDimT*Bubble_count(i)        ! The same for contrast source term , multiplied by length of T
                end if
            end do
            ! allocate temporary buffer in order to change the order of scatterers based on the ID of the proc
            ALLOCATE (tempbuffer(i_loc))
            tempbuffer = BubbleID(1:i_loc)
            call MPI_Allgatherv(tempbuffer, i_loc, MPI_INTEGER8, BubbleID, Bubble_count, Bubble_disps, MPI_INTEGER8, MPI_COMM_WORLD, iErr)
            DEALLOCATE (tempbuffer)

            ! Same proceedure but in order to reduce time, do this instead
            dBubbleLocationN = dBubbleLocationN((/BubbleID/), :)
            ScattererParams(cSourceParams%iScCloud)%R0 = ScattererParams(cSourceParams%iScCloud)%R0((/BubbleID/))

            ALLOCATE (dBubbleContrast(ibDimT*ScattererParams(cSourceParams%iScCloud)%N))
            iMemAllocated = iMemAllocated + PRODUCT(SHAPE(dBubbleContrast))*dpS

            ! Same idea with previous. This is done in order each processor to know the contrast and location of all the scatterers
            dBubbleContrast = 1D-30
            call MPI_Allgatherv(tempcon, ibDimT*i_loc, MPI_DOUBLE_PRECISION, dBubbleContrast, ibDimT*Bubble_count, Bubble_condisps, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, iErr)
            call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values

            !Wait until each processor has checked in , in order to finalize the number of bubbles in each processor and the proceed to next step
            iMemAllocated = iMemAllocated - PRODUCT(SHAPE(tempcon))*dpS
            DEALLOCATE (tempcon)

            !****************************************************** CALCULATE TIME SIGNATURE **********************************************
            call CalculateTimeSignature(cSpace, dBubbleContrast, ibDimT)
            !****************************************************** CALCULATE PRESSURE FIELD **********************************************
            write (acTemp, "('Reallocate the variables due to multiple populations ')"); call PrintToLog(acTemp, 2)
             
            if (.NOT. ALLOCATED(dCloudContrast))  iMemAllocated = iMemAllocated + ibDimT*dTotalCloudN*dpS
            if (.NOT. ALLOCATED(dCloudContrast)) ALLOCATE (dCloudContrast(ibDimT*dTotalCloudN))
            dCloudContrast( CloudIndex*ibDimT + 1 : (CloudIndex + ScattererParams(cSourceParams%iScCloud)%N)*ibDimT ) = dBubbleContrast
            CloudID( CloudIndex + 1 : CloudIndex + ScattererParams(cSourceParams%iScCloud)%N ) = BubbleID
            dCloudLocationN( CloudIndex + 1 : CloudIndex + ScattererParams(cSourceParams%iScCloud)%N, : ) = dBubbleLocationN

            CloudIndex = CloudIndex + ScattererParams(cSourceParams%iScCloud)%N

            call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values
            !*************************************************** END OF CALCULATE PRESSURE FIELD *******************************************
            write (acTemp, "('Save variables of population to file ')"); call PrintToLog(acTemp, 2)
            if (pcGrid%iProcID <= 0) then

                ! Save each Bubble's Position
                open (11, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//'Bubble_Location'//int2str(cSourceParams%iScCloud)//int2str(cModelParams%iIter)), status='REPLACE')
                do iBubble = 1, ScattererParams(cSourceParams%iScCloud)%N
                    write (11, *) dBubbleLocationN(iBubble, 1), dBubbleLocationN(iBubble, 2), dBubbleLocationN(iBubble, 3)
                end do
                close (11)

                ! For the polydisperse case , save each bubble's radius
                if (ScattererParams(cSourceParams%iScCloud)%Distribution == 'polydisperse') then
                    open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_Radius'//int2str(cSourceParams%iScCloud))//int2str(cModelParams%iIter)), status='REPLACE')
                    write (8, *) ScattererParams(cSourceParams%iScCloud)%R0
                    close (8)
                end if

            end if
 
            write (acTemp, "('Write the contrast of bubbles to output file')"); call PrintToLog(acTemp, 2)
            ! Multiple writes to split the work load in all the processes. Especially efficient in high number of scatterers
            ! These are unformatted which are a lot more efficient than FORMATTED files but it is more difficult to post-process
            MPI_FILENAME = trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//'Bubble_Contrast'//int2str(cSourceParams%iScCloud)//int2str(cModelParams%iIter))
            call MPI_WRITE_CONTRAST(cSpace, dBubbleContrast, MPI_FILENAME)

            ! This is for reducing memory
            open (8, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//'Bubble_iloc'//int2str(cSourceParams%iScCloud)//int2str(cSpace%cGrid%iProcID)), status='REPLACE')
            write (8, *) i_loc
            close (8)

            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(dBubbleContrast)) + PRODUCT(SHAPE(dBubbleLocationN)))*dpS - PRODUCT(SHAPE(BubbleID))*i4bs

            DEALLOCATE (dBubbleContrast)
            DEALLOCATE (BubbleID)
            DEALLOCATE (dBubbleLocationN)
            call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values

        end do
        ! Iniatilize pressure field with 0
        write (acTemp, "('Calculate the pressure field due to medium nonlinearity ')"); call PrintToLog(acTemp, 2)
        call NonlinContrastOperator_Ali(cSpace); ! Add Nonlinearities ! If beta=0 then it will result on 0
        call MPI_BARRIER(MPI_COMM_WORLD, iErr); 
        !    pcGrid%parD0 = 0.0D0
        !First, transform to W-domain (use small transform, no wraparound regions)
        write (acTemp, "('Calculate pressure field due to the Bubble Cloud')"); call PrintToLog(acTemp, 2)
        if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'around') then
            call CalculateCloudPressure_Around(cSpace, dCloudLocationN, dCloudContrast, CloudID, temploc, NeighbGridPoints, count_neighbs)
        else
            call CalculateCloudPressure(cSpace, dCloudLocationN, dCloudContrast)
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values
        !*************************************************** END OF CALCULATE PRESSURE FIELD *******************************************
        iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(dCloudContrast)) + PRODUCT(SHAPE(dCloudLocationN)))*dpS - PRODUCT(SHAPE(CloudID))*i4bs - PRODUCT(SHAPE(temploc))*i4bs

        DEALLOCATE (dCloudContrast)
        DEALLOCATE (CloudID)
        DEALLOCATE (dCloudLocationN)
        if (ALLOCATED(temploc)) DEALLOCATE (temploc)

        !****************************************************** FILTER AND TAPPER **********************************************

        write (acTemp, "('Tapper and filter contrast source')"); call PrintToLog(acTemp, 2)

        call ReorderDistr0ToDistr1(pcGrid)
        call ReorderDistr1ToDistr2(pcGrid)
        dTaperWidth = 1.0D0

        dTaperX = dTaperingWindow(cSpace%iDimX, cSpace%dDx, dTaperWidth, dTaperWidth)
        dTaperY = dTaperingWindow(cSpace%iDimY, cSpace%dDx, dTaperWidth, dTaperWidth)
        dTaperZ = dTaperingWindow(cSpace%iDimZ, cSpace%dDx, dTaperWidth, dTaperWidth)
        do iOmega = 0, pcGrid%iD2locN - 1
            do iIndexX = 0, cSpace%cGrid%iD2XL - 1
                do iIndexY = 0, cSpace%cGrid%iD2YL - 1
                    do iIndexZ = 0, cSpace%cGrid%iD2ZL - 1

                        iIndex = 1 + iOmega*cSpace%cGrid%iD2IS + iIndexX*cSpace%cGrid%iD2XS + iIndexY*cSpace%cGrid%iD2YS + iIndexZ*cSpace%cGrid%iD2ZS

                        cSpace%cGrid%pacD2(iIndex) = cSpace%cGrid%pacD2(iIndex)*(dTaperX(1 + iIndexX)*dTaperY(1 + iIndexY)*dTaperZ(1 + iIndexZ))

                    end do
                end do
            end do
        END DO
        call ReorderDistr2ToDistr1(pcGrid)
        call MapVtInv(cSpace)
        call ReorderDistr1ToDistr0(pcGrid); 
        call PrintToLog("End of Pressure Calculation", 2)
        call MPI_BARRIER(MPI_COMM_WORLD, iErr); 
    END SUBROUTINE BubbleContrastOperator

    SUBROUTINE GeneratePSCloud(cSpace, PointSourceType, dBubbleLocationN, Domain_Range, cubes)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 200220
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200220  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine GeneratePSCloud creates a cluster of microbubbles based
        !   on the minimum distance between the microbubbles and the domain that encloses
        !   all the bubbles ( maximum distance). This is done by dividing the cube into
        !   smaller cubes and randomly position each bubble inside based on the two
        !        aforementioned dimensions
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace             io   type(space)  Space from which the data array needs
        !                                        to be stored.
        !   MinDistance                i    i8b        Minimum distance between microbubbles
        !   dBubbleLocationN        r    dp                Location of each microbubble (X,Y,Z)
        !   Domain_Range                r    dp         Dimensions of Domain (2,3) that encloses
        !                                                all the microbubbles
        !
        ! *****************************************************************************
        type(Space), intent(in)                     ::              cSpace
        type(PointSourceCloudInput), intent(in)     ::              PointSourceType
        real(dp), intent(out)                       ::              Domain_Range(2, 3)
        real(dp), intent(inout)                     ::              dBubbleLocationN(:, :)

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   GridDist        r      Random number for the positioning of microbubbles
        !   i               i8b    Loop counter for each layer in Y + Z axis
        !   j               i8b    Loop counter for each layer in X axis
        !   Nx                i8b    Number of microbubbles in the X dimension of Domain
        !   Ny                i8b    Number of microbubbles in the Y dimension of Domain
        !   Nz                i8b    Number of microbubbles in the Z dimension of Domain
        !
        ! *****************************************************************************
        ! Cubes
        type(PointSourceCloudInput)                 ::              PointCloud
        integer(i4b)                                ::              iErr, ErrCode
        real(dp), allocatable                       ::              GridDist(:, :)
        integer(i8b), allocatable                    ::              BubbleID(:)
        real(dp)                                    ::              MaxDist(3), dLambdaMM, MinDist
        integer(i8b)                                ::              i, j
        integer(i8b)                                ::              Nx, Ny, Nz
        logical                                     ::              cubes
        character(len=1024)                         ::              acTemp

        ! Lines

        integer(i8b)                                ::              arBubble_i(ScattererParams(cSourceParams%iScCloud)%N), arBubble_j(ScattererParams(cSourceParams%iScCloud)%N)
        integer(i8b)                                ::              iBubble, counter, diff, iBubbleInside
        integer(i8b)                                ::              LocXYZ(ScattererParams(cSourceParams%iScCloud)%N), LocX(ScattererParams(cSourceParams%iScCloud)%N), LocY(ScattererParams(cSourceParams%iScCloud)%N), LocZ(ScattererParams(cSourceParams%iScCloud)%N)
        real(dp)                                    ::              randpos(ScattererParams(cSourceParams%iScCloud)%N)

        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries and storing of the data array to temporary file
        !
        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   ScattererInit
        !   RANDOM_SEED
        !   RANDOM_NUMBER
        !   ALLOCATE/DEALLOCATE
        !        SaveExtraBubbleSlices
        !
        ! =============================================================================
        PointCloud = PointSourceType
        if (cSpace%cGrid%iProcID == 0) then

            dLambdaMM = (cMediumParams.c0*1.0D3)/cModelParams%freq0; ! Normalization factor
            !unnormalized values of the resized Space
            Domain_Range(:, 1) = (PointCloud%ClusterDimsRatio(:, 1)*(cSpace%iDimX - 1)*cSpace%dDx + cSpace%iStartX*cSpace%dDx)*dLambdaMM
            Domain_Range(:, 2) = (PointCloud%ClusterDimsRatio(:, 2)*(cSpace%iDimY - 1)*cSpace%dDx + cSpace%iStartY*cSpace%dDx)*dLambdaMM
            Domain_Range(:, 3) = (PointCloud%ClusterDimsRatio(:, 3)*(cSpace%iDimZ - 1)*cSpace%dDx + cSpace%iStartZ*cSpace%dDx)*dLambdaMM
            MinDist = PointCloud%MinInBetweenDist*1e+3

            if (minval(Domain_Range(1, :) - (/cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ/)*cSpace%dDx*dLambdaMM) < 0.0D0 .OR. &
                maxval(Domain_Range(2, :) - (/cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ/)*cSpace%dDx*dLambdaMM) >= 0.0D0) then
                write (acTemp, "('Error in GeneratePSCloud, Domain out of boundaries')")
                call PrintToLog(acTemp, -1)
                call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
            end if

            if (cubes) then
                MaxDist = ((Domain_Range(2, 1) - Domain_Range(1, 1))*(Domain_Range(2, 2) - Domain_Range(1, 2))*(Domain_Range(2, 3) - Domain_Range(1, 3))/(size(dBubbleLocationN, 1)))**(1.0D0/3)   ! Maximum inbetween distance of bubbles in the specified Domain Range

                Nx = ceiling((Domain_Range(2, 1) - Domain_Range(1, 1))/MaxDist(1))
                Ny = ceiling((Domain_Range(2, 2) - Domain_Range(1, 2))/MaxDist(2))
                Nz = ceiling((Domain_Range(2, 3) - Domain_Range(1, 3))/MaxDist(3)); !ceiling(size(dBubbleLocationN, 1)*1.0/(Nx*Ny))

                if (Nx*Ny*Nz < size(dBubbleLocationN, 1)) then
                    write (acTemp, "('Error in GeneratePSCloud, less discretized points than Bubble number')")
                    call PrintToLog(acTemp, -1)
                    call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                end if

                ALLOCATE (GridDist(Nx*Ny*Nz, 3))

                MaxDist = (Domain_Range(2, :) - Domain_Range(1, :))/(/Nx, Ny, Nz/)
                if (minval(MaxDist) < MinDist) then
                    write (acTemp, "('Error ! Bubble Cluster size too small for ',I5,' Bubbles.')") size(dBubbleLocationN, 1)
                    call PrintToLog(acTemp, -1)
                    write (acTemp, "('Max Bubble number for this Cluster size',I5,' Bubbles.')") size(dBubbleLocationN, 1)/((MinDist/MaxDist(1))**(1.0D0/3))
                    call PrintToLog(acTemp, -1)
                    call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                end if

                call RANDOM_SEED(PUT=[(1, i=1, Nx*Ny*Nz*3)])
                call RANDOM_NUMBER(GridDist)

                do i = 0, Ny*Nz - 1
                    GridDist(1 + i*Nx:Nx + i*Nx, 1) = Domain_Range(1, 1) + 1.0D0*(1.0 - GridDist(1 + i*Nx:Nx + i*Nx, 1))*MinDist + &
                                                      1.0D0*MaxDist(1)*GridDist(1 + i*Nx:Nx + i*Nx, 1) + MaxDist(1)*([(j, j=0, Nx - 1)])
                    GridDist(1 + i*Nx:Nx + i*Nx, 2) = Domain_Range(1, 2) + 1.0D0*(1.0 - GridDist(1 + i*Nx:Nx + i*Nx, 2))*MinDist + &
                                                      1.0D0*MaxDist(2)*GridDist(1 + i*Nx:Nx + i*Nx, 2) + MaxDist(2)*mod(i, Ny)
                    GridDist(1 + i*Nx:Nx + i*Nx, 3) = Domain_Range(1, 3) + 1.0D0*(1.0 - GridDist(1 + i*Nx:Nx + i*Nx, 3))*MinDist + &
                                                      1.0D0*MaxDist(3)*GridDist(1 + i*Nx:Nx + i*Nx, 3) + MaxDist(3)*(i/Ny)
                end do

                ALLOCATE (BubbleID(Nx*Ny*Nz))
                BubbleID = (/(i, i=1, Nx*Ny*Nz)/)
                call FISHER_YATES_SHUFFLE(BubbleID, Nx*Ny*Nz)
                dBubbleLocationN = GridDist(BubbleID(1:size(dBubbleLocationN, 1)), :)
                DEALLOCATE (BubbleID)

                DEALLOCATE (GridDist)

            else  ! Lines random type

                Nx = INT((Domain_Range(2, 1) - Domain_Range(1, 1))/MinDist, 8) + 1
                Ny = INT((Domain_Range(2, 2) - Domain_Range(1, 2))/MinDist, 8) + 1
                Nz = INT((Domain_Range(2, 3) - Domain_Range(1, 3))/MinDist, 8) + 1

                if (Nx*Ny*Nz < PointCloud%N) then
                    write (acTemp, "('Too big minimum in-between distance or too small domain')")
                    call PrintToLog(acTemp, -1)
                    call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                end if

                arBubble_i = 1
                LocXYZ = 0
                call RANDOM_SEED(PUT=[(1, i=1, 10*PointCloud%N*3)]) !! Generate the same random number , when ran
                call RANDOM_SEED()        !GEnerate different random numbers for every run
                call RANDOM_NUMBER(randpos)
                arBubble_i = INT(randpos*Nx*Ny*Nz, 8)
                call I_unista(arBubble_i, iBubble)
                LocXYZ(1:iBubble) = arBubble_i(1:iBubble)

                call RANDOM_SEED(PUT=[(1, i=1, 10*PointCloud%N*3)])
                do while (iBubble < PointCloud%N)
                    counter = 0; arBubble_i = 0; diff = 0; 
                    !                        call RANDOM_SEED()        !GEnerate different random numbers for every run
                    call RANDOM_NUMBER(HARVEST=randpos)
                    diff = PointCloud%N - iBubble
                    arBubble_i(1:diff) = INT(randpos(1:diff)*Nx*Ny*Nz, 8)
                    LocXYZ(iBubble + 1:iBubble + diff) = arBubble_i(1:diff)
                    call I_unista(LocXYZ(1:iBubble + diff), counter)
                    iBubble = counter
                end do

                LocX = mod(LocXYZ, Nx); 
                arBubble_i = (LocXYZ - LocX)/Nx; 
                LocY = mod(arBubble_i, Ny); 
                arBubble_j = (arBubble_i - LocY)/Ny; 
                LocZ = arBubble_j; 
                dBubbleLocationN(:, 1) = LocX*MinDist + Domain_Range(1, 1); 
                dBubbleLocationN(:, 2) = LocY*MinDist + Domain_Range(1, 2); 
                dBubbleLocationN(:, 3) = LocZ*MinDist + Domain_Range(1, 3); 
            end if

            Domain_Range(1, :) = minval(dBubbleLocationN, 1)
            Domain_Range(2, :) = maxval(dBubbleLocationN, 1)

            dBubbleLocationN(:, 1) = dBubbleLocationN(:, 1)/dLambdaMM - cSpace%iStartX*cSpace%dDx
            dBubbleLocationN(:, 2) = dBubbleLocationN(:, 2)/dLambdaMM - cSpace%iStartY*cSpace%dDx
            dBubbleLocationN(:, 3) = dBubbleLocationN(:, 3)/dLambdaMM - cSpace%iStartZ*cSpace%dDx
            dBubbleLocationN = dBubbleLocationN + EPSILON(1.0D0) ! This is added because when the scatterer is exactly positioned at a gridpoint , the code is stucked

            if (minval((/minval(dBubbleLocationN, 1), &
                         (cSpace%iDimX - 1)*cSpace%dDx - maxval(dBubbleLocationN(:, 1), 1), &
                         (cSpace%iDimY - 1)*cSpace%dDx - maxval(dBubbleLocationN(:, 2), 1), &
                         (cSpace%iDimZ - 1)*cSpace%dDx - maxval(dBubbleLocationN(:, 3), 1)/)) <= 0) then
                write (acTemp, "('Error in GeneratePSCloud, Bubble positioned outside domain boundaries')")
                call PrintToLog(acTemp, -1)
                call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
            end if

        end if

        call MPI_BCAST(dBubbleLocationN, PointCloud%N*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)

    END SUBROUTINE GeneratePSCloud

    SUBROUTINE RNDGeneratePSCloud_Eff(cSpace, PointSourceType, dBubbleLocationN, Domain_Range)

        USE MKL_VSL
        USE MKL_VSL_TYPE
        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 200220
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200220  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine RNDGeneratePSCloud creates a cluster of microbubbles based
        !   on the minimum distance between the microbubbles and the domain that encloses
        !   all the bubbles ( maximum distance). This is done by making a grid with the max
        !        distance in each dimension divided by the minimum mutual distance. Then generate
        !        randomly N points in this grid equally to the number of bubbles. Take the unique
        !        values from these numbers and if these are less than N , iterate until all the
        !        unique random numbers equal the number of microbubbles N.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace             io   type(space)  Space from which the data array needs
        !                                        to be stored.
        !   MinDistance                i    i8b        Minimum distance between microbubbles
        !   dBubbleLocationN        r    dp                Location of each microbubble (X,Y,Z)
        !
        ! *****************************************************************************
        real(dp), intent(inout)                     ::      dBubbleLocationN(:, :)
        type(Space), intent(in)                     ::      cSpace
        type(PointSourceCloudInput), intent(in)     ::      PointSourceType
        real(dp), intent(out)                       ::      Domain_Range(2, 3)

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   Domain_Range                r    dp         Dimensions of Domain (2,3) that encloses all the microbubbles
        !   randpos                                r    dp         vector with random values in the generated grid
        !   MinAllowableDist    r    dp                    Value of minimum mutual distance
        !   dLambdaMM                        r    dp                     value of wavelength , normalization factor
        !   Nx,Ny,Nz                        i8b                            Number of gridpoints in new generated domain in x, y and z
        !   iBubble,counter                i8b                            Index counter
        !   diff                                i8b                           difference between N and the duplicate random numbers
        !   i,j                                        r    dp         temporary vectors in order for the proceedure to be clear
        !   LocX,LocY,LocZ                i                                vectors with the location of bubbles in X,Y,Z in the generated domain with mindist
        !   LocXYZ                                i                                vector with one values for a triplet of (X,Y,Z)
        !
        ! *****************************************************************************
        type(PointSourceCloudInput)                 ::      PointCloud
        integer(i4b)                                ::      iErr, ErrCode, errstatus
        real(dp)                                    ::      MultFactor
        integer                                     ::      Bcounter, iBubble, BubbleIdx, i, j, MaxBubbles, temp
        integer                                     ::      Num_counter(cSpace%cGrid%iProcN), Bubble_Outliers(cSpace%cGrid%iProcN)
        real(dp)                                    ::      dLambdaMM, Dist, maxRad(3)

        integer(i8b), allocatable                   ::      Bubble_count(:), BubbleID(:), BubbleCube(:)
        real(dp), allocatable                       ::      randpos(:, :), RadFactor(:)

        character(len=1024)                         ::      acTemp, BCshape
        !Geometric series parameters
        real(dp)                                    ::      a, r
        integer(i8b)                                ::      n, S, iProc, BCshape_int(3)
        integer(i8b), allocatable                   ::      GridDist(:)

        TYPE(VSL_STREAM_STATE)                      ::      STREAM2

        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries and storing of the data array to temporary file
        !
        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   I_Unista
        !   RANDOM_SEED
        !   RANDOM_NUMBER
        !   ALLOCATE/DEALLOCATE
        !        SaveExtraBubbleSlices
        !
        ! =============================================================================

        MultFactor = 2.5D0

        PointCloud = PointSourceType
        S = INT(MultFactor*PointCloud%N, 8)
        n = S; 
        ALLOCATE (BubbleID(S), randpos(S, 3), RadFactor(S))

        do while (S - n < PointCloud%N)

            dLambdaMM = (cMediumParams.c0*1.0D3)/cModelParams%freq0; ! Normalization factor
            !unnormalized values of the resized Space
            Domain_Range(:, 1) = (PointCloud%ClusterDimsRatio(:, 1)*(cSpace%iDimX - 1) + cSpace%iStartX)*cSpace%dDx*dLambdaMM
            Domain_Range(:, 2) = (PointCloud%ClusterDimsRatio(:, 2)*(cSpace%iDimY - 1) + cSpace%iStartY)*cSpace%dDx*dLambdaMM
            Domain_Range(:, 3) = (PointCloud%ClusterDimsRatio(:, 3)*(cSpace%iDimZ - 1) + cSpace%iStartZ)*cSpace%dDx*dLambdaMM

            if (trim(ScattererParams(cSourceParams%iScCloud)%Distribution) == 'monodisperse') then
                RadFactor(1:S - n) = PointCloud%R0(1)
            elseif (trim(ScattererParams(cSourceParams%iScCloud)%Distribution) == 'polydisperse') then
                errstatus = vslnewstream(STREAM2, VSL_BRNG_MCG31, 1)
                errstatus = vdrnggamma(VSL_RNG_METHOD_GAMMA_GNORM_ACCURATE, STREAM2, INT(S), RadFactor, 1.2D0, 0.0D0, 1.0D0)

                RadFactor = (RadFactor - MINVAL(RadFactor))
                RadFactor = RadFactor/MAXVAL(RadFactor)*(ScattererParams(cSourceParams%iScCloud)%PDRange(2) - ScattererParams(cSourceParams%iScCloud)%PDRange(1)) + ScattererParams(cSourceParams%iScCloud)%PDRange(1)
            end if
            if (minval(Domain_Range(1, :) - (/cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ/)*cSpace%dDx*dLambdaMM) < 0.0D0 .OR. &
                maxval(Domain_Range(2, :) - (/cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ/)*cSpace%dDx*dLambdaMM) >= 0.0D0) then
                write (acTemp, "('Error in GeneratePSCloud, Domain out of boundaries')")
                call PrintToLog(acTemp, -1)
                call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
            end if

            !=================================================== RANDOM LOCATION OF SCATTERERS =====================================
            call RANDOM_SEED(PUT=[(1, i=1, 10*S*3)])  ! Random but the same at each run of INCS
            ! Here we generate random points that are located inside the boundaries that we set
            ! call RANDOM_SEED() ! Fully random , even at each run of INCS
            call RANDOM_NUMBER(randpos)
            randpos(1:n, 1) = randpos(1:n, 1)*(Domain_Range(2, 1) - Domain_Range(1, 1))
            randpos(1:n, 2) = randpos(1:n, 2)*(Domain_Range(2, 2) - Domain_Range(1, 2))
            randpos(1:n, 3) = randpos(1:n, 3)*(Domain_Range(2, 3) - Domain_Range(1, 3))
            call MPI_BCAST(randpos, S*3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)

            ! Here we split our domain into n equal parts so each cpu to compare only the point scatterers that are located inside each
            ! own domain, we also add forward the minimum inbetween distance needed so we can compare with the neighbors of the next succesive
            ! domain that can be close to the scatterers in the cpu's domain.
            Dist = (Domain_Range(2, 3) - Domain_Range(1, 3))/cSpace%cGrid%iProcN; 
            Bcounter = COUNT(randpos(1:n, 3) >= cSpace%cGrid%iProcID*Dist .AND. &
                             randpos(1:n, 3) <= (cSpace%cGrid%iProcID + 1)*Dist + PointCloud%MinInBetweenDist*1.0D+3)

            BubbleID(1:n) = (/(i, i=1, n)/)
            ! Here we create our integers ,so each BubbleID
            ALLOCATE (BubbleCube(Bcounter), Bubble_count(INT(Bcounter, 8)))
            BubbleCube = PACK(BubbleID, (randpos(1:n, 3) >= cSpace%cGrid%iProcID*Dist .AND. &
                                         randpos(1:n, 3) <= (cSpace%cGrid%iProcID + 1)*Dist + PointCloud%MinInBetweenDist*1.0D+3))

            ! ============================== FIND HOW MANY DO NOT FULFILL THE REQUIREMENT OF THE HYPOTHESIS OF DISTANCE (REJECTED) ==========================

            Bcounter = 1; 
            Bubble_count = 0; 
            maxRad(2) = MINVAL(Domain_Range(2, 1:3:2) - Domain_Range(1, 1:3:2))
            maxRad(3) = MINVAL(Domain_Range(2, :) - Domain_Range(1, :))
            BCshape = 'cylinder'
            if (cSourceParams%iScCloud==2) BCShape = 'cube'
            BCshape_int = 0; 
            BCshape_int(1) = ABS(BCshape .EQ. 'cube'); BCshape_int(2) = ABS(BCshape .EQ. 'cylinder'); BCshape_int(3) = ABS(BCshape .EQ. 'sphere')
            write (acTemp, "('The scatterers are distributed in a :', A)") BCshape_int(1); call PrintToLog(acTemp, 3); 

            do iBubble = 1, PRODUCT(SHAPE(BubbleCube))

                if (ANY(sqrt((randpos(BubbleCube(iBubble), 1) - randpos(BubbleCube(1:iBubble - 1), 1))**2 + &
                             (randpos(BubbleCube(iBubble), 2) - randpos(BubbleCube(1:iBubble - 1), 2))**2 + &
                             (randpos(BubbleCube(iBubble), 3) - randpos(BubbleCube(1:iBubble - 1), 3))**2) - &
                        PointCloud%MinInBetweenDist*1.0D+3 < 0) .OR. &
                    (SQRT(SUM((randpos(BubbleCube(iBubble), 1:3:2) - (Domain_Range(2, 1:3:2) + Domain_Range(1, 1:3:2))/2.0D0 + Domain_Range(1, 1:3:2))**2)) >= maxRad(2)/2.0D0)*BCshape_int(2) .OR. &
                    (SQRT(SUM((randpos(BubbleCube(iBubble), :) - (Domain_Range(2, :) + Domain_Range(1, :))/2.0D0 + Domain_Range(1, :))**2)) >= maxRad(3)/2.0D0)*BCshape_int(3)) then

                    Bubble_count(Bcounter) = BubbleCube(iBubble)
                    Bcounter = Bcounter + 1
                end if

            end do

            ! ================================== TRANSFER THEM FROM EACH PROCESSOR =======================================================
            call MPI_Allgather(Bcounter - 1, 1, MPI_INTEGER, Num_counter, 1, MPI_INTEGER, MPI_COMM_WORLD, iErr) ! Array with the number of bubble outliers
            Bubble_outliers = 0; 
            do i = 1, cSpace%cGrid%iProcN
                if (Num_counter(i) > 0) then
                    Bubble_outliers(i) = sum(Num_counter(1:i)) - Num_counter(i)                ! Displacement variable in order to gather the value from every process
                end if
            end do

            ALLOCATE (GridDist(Bcounter - 1))
            GridDist = Bubble_count(1:Bcounter - 1)
            DEALLOCATE (Bubble_count)

            ALLOCATE (Bubble_count(sum(Num_counter)))
            Bubble_count = 0
            call MPI_Allgatherv(GridDist, Bcounter - 1, MPI_INTEGER8, Bubble_count, Num_counter, Bubble_outliers, MPI_INTEGER8, MPI_COMM_WORLD, iErr)
            DEALLOCATE (GridDist)
            call I_unista(Bubble_count, n)
            ! n is the number of the rejected unique points in space , here we remove the doubles
            ! (because of the common points between different cpus in order to compare for the neighbours )

            ALLOCATE (GridDist(S))
            GridDist = 1; GridDist(Bubble_count(1:n)) = 0; 
            GridDist(1:S - n) = PACK(BubbleID, GridDist == 1); 
            DEALLOCATE (Bubble_count)
            ! ======================================= GENERATE THE LOCATION MATRIX FROM THE CORRECT VALUES ===================================

            ! if there is no scatterer too close to the others, change the variables of interest and proceed without any more actions
            if (S - n >= PointCloud%N) then
                call FISHER_YATES_SHUFFLE(BubbleID, S - n)
                dBubbleLocationN(1:PointCloud%N, :) = randpos(GridDist(BubbleID(1:PointCloud%N)), :)
                PointCloud%R0(1:PointCloud%N) = RadFactor(GridDist(BubbleID(1:PointCloud%N)))
            else
                write (acTemp, "('More iterations are needed')"); call PrintToLog(acTemp, 3); 
                write (acTemp, "('No. of accepted MBs is :', I<INT(log10(real(PointCloud%N,dp)))+1>)") S - n; call PrintToLog(acTemp, 3); 
            end if

            MultFactor = MultFactor*2.0D0; 
            DEALLOCATE (GridDist, BubbleID, BubbleCube)
            DEALLOCATE (randpos, RadFactor)
        end do   

        dBubbleLocationN(:, 1) = dBubbleLocationN(:, 1) + Domain_Range(1, 1) - 0.8
        dBubbleLocationN(:, 2) = dBubbleLocationN(:, 2) + Domain_Range(1, 2)
        dBubbleLocationN(:, 3) = dBubbleLocationN(:, 3) + Domain_Range(1, 3)

        ! ================================================================================================================================

        ! dBubbleLocationN(1,:)=  (/0.0 ,   0.0 ,    5.187/) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(1,:)=  (/-0.0366470588235294 ,   -0.0366470588235294 ,    3.63602941176471/) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(2,:)=  (/-0.0366470588235294 ,   -0.0366470588235294 ,    5.08194117647059/) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(1,:)=  (/-0.0726470588235294 ,   -0.0726470588235294 ,    3.54033333333333 /) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(1,:)=  (/-0.0823333333333333 ,   -0.0823333333333333 ,    3.458 /) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(2,:)=  (/-0.0726470588235294 ,   -0.0726470588235294 ,    3.63235294117647/) ! At a gridpoint
        ! dBubbleLocationN(1,:)=  (/0.599338235294118 ,   0.108970588235294 ,    5.17610294117647/) ! At a gridpoint
        ! dBubbleLocationN(2,:)=  (/0.572095587323694 ,   0.0726470588235294 ,    3.148860381070307/) ! Βetween gridpoints
        ! dBubbleLocationN(1,:)=  (/0.576333333333333  ,   0.0823333333333333  ,    5.187/)  ! Middle of gridpoints
        ! dBubbleLocationN(2,:)=  (/0.653823529411765 ,  0.0726470588235294 ,   1.45294117647059/) ! At a gridpoint  @ 1E6 freq0 , 10 Fnyq
        ! dBubbleLocationN(1,:)=  (/0.0 ,   0.0 ,   4.11666666666667/) ! At a gridpoint  @ 1E6 freq0 , 6 Fnyq
        ! dBubbleLocationN(1,:)=  (/.576333333333333,0.0823333333333333,5.187  /)
        ! dBubbleLocationN(2,:)=  (/-0.411677777777,0.0823333333333333,6.2573 /)
        ! dBubbleLocationN(1, :) = (/-0.0412222222222222, -0.0412222222222222, 5.7221/) ! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq

        dBubbleLocationN = dBubbleLocationN + EPSILON(1.0D0)
        ! This is added because when the scatterer is exactly positioned at a gridpoint , the code is stucked

        call MPI_BCAST(dBubbleLocationN, INT(PointCloud%N*3), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)
        call MPI_BCAST(PointCloud%R0, INT(PointCloud%N), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, iErr)

        ! if (minval((/ &
        !            Domain_Range(2, :) - maxval(dBubbleLocationN, 1), &
        !            minval(dBubbleLocationN, 1) - Domain_Range(1, :)/)) < 0) then

        !     write (acTemp, "('Error in GeneratePSCloud, Bubble positioned outside domain boundaries')")
        !     call PrintToLog(acTemp, -1)
        !     call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
        ! end if
        Domain_Range(1, :) = minval(dBubbleLocationN, 1)
        Domain_Range(2, :) = maxval(dBubbleLocationN, 1)

        dBubbleLocationN(:, 1) = dBubbleLocationN(:, 1)/dLambdaMM - cSpace%iStartX*cSpace%dDx
        dBubbleLocationN(:, 2) = dBubbleLocationN(:, 2)/dLambdaMM - cSpace%iStartY*cSpace%dDx
        dBubbleLocationN(:, 3) = dBubbleLocationN(:, 3)/dLambdaMM - cSpace%iStartZ*cSpace%dDx

    END SUBROUTINE RNDGeneratePSCloud_Eff

    SUBROUTINE ProcBubbleCount(cSpace, dBubbleLocationN, tempcon, iProc_NumBubbles, ibDimT)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 201128
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200514  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine ProcBubbleCount counts how many bubbles are located in the domain
        !   of each cpu that is used  (Total domain is splitted in N subdomains , each for
        !         each cpu used).Besides this , the values of the time signature of each scatterer
        !        is loaded . Then, only the values of the bubbles inside each subdomain are kept.
        !        The rest are discarded.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace             io   type(space)  Space from which the data array needs
        !                                        to be stored.
        !   dBubbleLocationN   r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   tempcon            r   dp   a vector of the bubbles' time signature
        !   iProc_NumBubbles   i        number of bubbles per process
        !
        ! *****************************************************************************

        real(dp), intent(inout)                 ::      dBubbleLocationN(:, :)
        real(dp), allocatable, intent(inout)    ::      tempcon(:)
        type(Space), intent(in)                 ::      cSpace
        integer, intent(inout)                  ::      iProc_NumBubbles, ibDimT

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   iErr            i8b    Error number
        !   BubbleID        i8b    No. of Bubble
        !   iBubble         i8b    index variable
        !   XYZIndex        i8b    one value to describe (x,y,z) normalized coordinates
        !   iBeamIndex      i8b    Index of a point in the beam domain
        !   iDimT           i      Time Length ( Should be integer)
        !        dBubbleCon                dp           vector that is used for transferring the time signature array
        !
        ! *****************************************************************************

        integer(i4b)                            ::      iErr, FILE_HANDLE, ErrCode
        integer(i8b)                            ::      iBubble, XYZIndex, i, iBeamIndex(3), FILE_IN, FILE_OUT, length_Bubble
        real(dp), allocatable                   ::      dBubbleCon(:)
        character(len=1024)                     ::      acTemp, MPI_FILENAME

        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   FFTW
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   If there are bubbles in the processor ( this is done after iter 1)
        !02          load the saved time signatures
        !03   Initialize the number of bubbles inside each cpu domain
        !04          Count for each scatterer inside the domain
        !04          This is done by generating one value for each triplet of (X,Y,Z)
        !04          and checking if this value is between 0 and max number of gridpoints
        !05          if 2nd iter or higher, transfer the time signatures of the bubbles
        !05          located in each cpu domain.
        !
        ! =============================================================================

        ! ! ! ! !Calculate Number of bubbles below half wavelength of Fnyq
        ! ! ! integer(i8b)                        ::                count5(ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN,275), Bubble_Loc
        ! ! ! integer(i8b)                        ::                  length(4), i_x,i, diff_z(ScattererParams(cSourceParams%iScCloud)%N)
        ! ! ! integer(i8b),allocatable   :: count51(:,:)
        ! ! ! real(dp)                             ::                 diff(ScattererParams(cSourceParams%iScCloud)%N),dLambdaMM

        ! ! ! dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams%freq0;           ! Normalization factor
        ! ! ! count5=0;
        ! ! ! count51=0;
        ! ! ! length_Bubble = ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN
        ! ! ! if (cSpace%cGrid%iProcID == cSpace%cGrid%iProcN-1) length_Bubble = ScattererParams(cSourceParams%iScCloud)%N - ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN * cSpace%cGrid%iProcID
        ! ! ! ALLOCATE(count51(length_Bubble,1500))
        ! ! ! diff_z = 0 ;
        ! ! ! do iBubble =  1, length_Bubble
        ! ! ! diff=0;
        ! ! ! Bubble_Loc = iBubble+ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN*cSpace%cGrid%iProcID
        ! ! ! diff = sqrt((dBubbleLocationN(Bubble_Loc,1) - dBubbleLocationN(:,1))**2 +(dBubbleLocationN(Bubble_Loc,2) - dBubbleLocationN(:,2))**2+ (dBubbleLocationN(Bubble_Loc,3) - dBubbleLocationN(:,3))**2)

        ! ! ! diff_z = INT(ceiling(diff/cSpace%dDx) * SIGN(1.0,dBubbleLocationN(Bubble_Loc,3) - dBubbleLocationN(:,3) ) )

        ! ! ! do i  = 1,maxval(abs(diff_z))
        ! ! ! count51(iBubble,2*i-1) = count(diff_z .EQ. i)
        ! ! ! count51(iBubble,2*i)   = count(diff_z .EQ. -1*i)
        ! ! ! enddo

        ! ! ! if(mod(iBubble,3000)==0) write(*,*) iBubble        , cSpace%cGrid%iProcId
        ! ! ! enddo
        ! ! ! open (7, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('Bubble_Count')//int2str(cModelParams%iIter)//int2str(cSpace%cGrid%iProcID)//'.bin'),status = 'UNKNOWN', action='write',position='append')
        ! ! ! write(7,*) count51
        ! ! ! close(7)

        !***********************************CALCULATE HOW MANY MICROBUBBLES ARE LOCATED IN EACH PROCESSOR***********************************************

        write (acTemp, "('Count number of microbubbles that exist in every task') "); call PrintToLog(acTemp, 2)
        MPI_FILENAME = trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('Bubble_Contrast'//int2str(cSourceParams%iScCloud))//int2str(cModelParams%iIter - 1))
        call MPI_FILE_OPEN(MPI_COMM_WORLD, MPI_FILENAME, MPI_MODE_RDWR, MPI_INFO_NULL, FILE_HANDLE, iErr)

        if (iProc_NumBubbles > 0) then
            ALLOCATE (dBubbleCon(ibDimT*ScattererParams(cSourceParams%iScCloud)%N))
            call MPI_FILE_READ(FILE_HANDLE, dBubbleCon, INT(ibDimT*ScattererParams(cSourceParams%iScCloud)%N, KIND=8), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, iErr)
        end if
        call MPI_FILE_CLOSE(FILE_HANDLE, iErr)

        iProc_NumBubbles = 0
        do iBubble = 1, ScattererParams(cSourceParams%iScCloud)%N
            ! Determine Beam Index of the point source
            ! In ParnacDataDef it is specified that cSpace%iStartT/X/Y/Z should be included
            ! The difference here is that the dBubbleLocationN parameter is in the normalized space
            ! So if we change this in real space then
            ! dBubbleLocation = (dBubbleLocationN + (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaNN
            ! So there is no difference if the starting position is not added in all the equations
            iBeamIndex(3) = nint(dBubbleLocationN(iBubble, 3)/cSpace%dDx)
            iBeamIndex(1) = nint(dBubbleLocationN(iBubble, 1)/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex(3))
            iBeamIndex(2) = nint(dBubbleLocationN(iBubble, 2)/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex(3))

            !Find the Index in the cSpace%cGrid%aiD0Loc of the BubbleLocation
            XYZIndex = iBeamIndex(3)*cSpace%iDimY*cSpace%iDimX + iBeamIndex(2)*cSpace%iDimX + iBeamIndex(1) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN
            if (XYZIndex >= 0 .AND. XYZINdex < cSpace%cGrid%iD0LocN) then
                iProc_NumBubbles = iProc_NumBubbles + 1
                if (ALLOCATED(tempcon)) tempcon(1 + (iProc_NumBubbles - 1)*ibDimT:iProc_NumBubbles*ibDimT) = dBubbleCon(1 + (iBubble - 1)*ibDimT:iBubble*ibDimT)
            end if

        end do

        if (ALLOCATED(dBubbleCon)) DEALLOCATE (dBubbleCon)
    END SUBROUTINE ProcBubbleCount

    SUBROUTINE TrilinInterp_SNDRCV(cSpace, dBubbleLocationN, i_loc, tempcon, temploc, BubbleID, count_neighbs)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 201128
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function CalculateCloudPressure calculates the pressure of the total bubble
        !        cloud.  This is done in two phases. The first phase is for each iteration ,except
        !   the last one. For this phase, all the neighbouring point of microbubbles are used
        !        and then the pressure is calculated in these points by all the microbubbles.
        !        In order to keep memory low , this is done by each cpu using the neighbouring points
        !        that are located inside its subdomain. The second phase is the last iteration.
        !        At this iteration ,the pressure is calculated at each point of the domain.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                       io   type(space)  The space in which the point has to be given
        !   dBubbleLocationN         r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   dBubbleContrast                 r   dp   a vector with the time signatures of the microbubbles
        !   NeighbGPIndex                 i   i8b  a vector with the location of the neighbouring points of microbubbles in each cpu domain
        !   i_neighbgp                          i   i8b  number of neighbouring points ( Do not get confused, not all elements of NeighbGPIndex
        !                                                                          have the values , only 1:i_neighbgp)
        !   last_iter                         l        logical value which states if it is the last iteration or not
        !
        ! *****************************************************************************
        type(Space), intent(inout)      ::      cSpace
        real(dp), intent(inout)         ::      dBubbleLocationN(:, :), tempcon(:)
        integer, intent(inout)          ::      i_loc
        integer(i8b), intent(inout)     ::      BubbleID(:), temploc(:, :), count_neighbs

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   P_Scattered                       dpc   Scattered Pressure Matrix for MKL
        !   dSincMult                             dpc   a matrix that contains spatially filtered dirac
        !   NeighbLocInd                        i8b   a vector which contains the indexes of the neighbouring points
        !   dGlobBXYZ                        r     a matrix of the gridpoints location of each cpu domain (same dim as dBubbleLocationN)
        !   M, ibDimT                               i     integer value for MKL , Time length
        !   N, GridPointsNum            i          integer value for MKL , Space length
        !   K, BLayer                                i           integer value for MKL , Bubbles per layer
        !   iIndex, i, iBubble            i8b   index variables
        !   B_Div                                           i8b   divisor of total numbers of variables
        !        dLambdaMM                                dp          value of wavelength , normalization factor
        !        dK_c                                        dp          Kutoff frequency
        !
        ! ************************************************************************
        integer                         ::      iErr, iStat(MPI_STATUS_SIZE), ErrCode, log_isnan
        integer(i8b)                    ::      iBeamIndex(3), iBeamIndexN(3), XYZindex, iTargetIndex(3), iBubble
        integer(i8b)                    ::      trilin(2), ones(size(trilin, 1)**3, 3), XYZ_i1(size(ones, 1))  ! This is for the Trilinear Interpolation
        integer(i8b)                    ::      NeighbGridPoints(size(temploc, 1)/ScattererParams(cSourceParams%iScCloud)%N, 3), p(size(NeighbGridPoints, 1)**(1.0D0/3.0D0)), XYZ_i2(size(p, 1)**3)
        integer                         ::      BubbleProc, ibDimT
        integer(i8b)                    ::      iDimT, iDimX, iDimY, iDimZ, iDimW, i, j, l
        integer(i8b), parameter         ::      OS_Factor = 1

        real(dp)                        ::      dGlobXYZ(3), GlobalPos(size(ones, 1), 3), xyz_factor(3), diff(size(ones, 1), 3), dRadius(1)
        real(dp)                        ::      RealPressure(OS_Factor*cSpace%iDimT)
        complex(dpc)                    ::      RealPressureC(OS_Factor*cSpace%iDimT)
        real(dp)                        ::      RealPressure_i1(size(ones, 1), cSpace%iDimT), RealPressure_Own(size(ones, 1), cSpace%iDimT), Interp_PressureBL(cSpace%iDimT)
        real(dp), allocatable           ::      temp_globloc(:), LocPres(:, :)
        real(dp)                        ::      dLambdaMM, dK_c
        real(dp)                        ::      dDOmega, dFFTFactor, dMultFactor(cSpace%iDimT)

        character(len=1024)             ::      actemp
        complex(dpc), allocatable       ::      InputValWT(:), OutputValWT(:)

        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   FFTW
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   First translate the NeighbLocInd into normalized coordinates in (x,y,z) if not in last iteration
        !01          and set the points for calculation only to neighbouring, otherwise calculate for every point of domain
        !02          Set a good choice for the number of Bubbles per Layer ( should check memory usage)
        !03          Iterate for the number of Layers (Total Bubbles / No. Bubbles per Layer)
        !04          Calculate dirac function matrix
        !05          Do the matrix multiplication based on the space points for calculation
        !06          if not in last iteration , transfer the calculated values to a matrix for easy replacement to cSpace var
        !07          else just replace the calculated values
        !
        ! =============================================================================

        iDimT = cSpace%iDimT   !t dimensions
        iDimX = cSpace%iDimX   !x dimensions
        iDimY = cSpace%iDimY   !y dimensions
        iDimZ = cSpace%iDimZ   !z dimensions

        iDimW = iDimT/2 + 1
        if (i_loc == 0) then
            ibDimT = 0; 
        else
            ibDimT = PRODUCT(SHAPE(tempcon))/i_loc; 
        end if

        dK_c = pi/cSpace%dDx
        dDOmega = two_pi*2.0_dp*cSpace%dFnyq/real(iDimT, dp)
        dMultFactor = (/(i, i=0, iDimW - 1), (i - iDimT, i=iDimW, iDimT - 1)/)*dDOmega

        ALLOCATE (LocPres(i_loc, size(ones, 1)*(3 + iDimT)), temp_globloc(3 + cSpace%iDimT))
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(temp_globloc)) + PRODUCT(SHAPE(LocPres)))*dpS

        temp_globloc = 1.0D-30
        LocPres = 1D-30
        GlobalPos = 1D-30
        i_loc = 0
        count_neighbs = 0

        ! This is to generate all the 8 neighbouring points for the trilinear interpolation
        ! ones(:,1) = (/-1, 0,-1, 0,-1, 0,-1,0/);
        ! ones(:,2) = (/-1,-1, 0, 0,-1,-1, 0,0/);
        ! ones(:,3) = (/-1,-1,-1,-1, 0, 0, 0,0/)
        trilin = (/(i, i=-INT(CEILING(size(trilin, 1)*1.0D0/2)), INT(FLOOR(size(trilin, 1)*1.0D0/2 - 1)))/); 
        ones = RESHAPE((/(((trilin(i), trilin(j), trilin(l), i=1, size(trilin, 1)), j=1, size(trilin, 1)), l=1, size(trilin, 1))/), (/size(trilin, 1)**3, 3/), order=(/2, 1/))
        ! This is to produce results on sufficient gridpoints to have accurate results. By increasing the size of p (even numbers) , the accuracy increases

        p = (/(i, i=-INT(CEILING(size(p, 1)*1.0D0/2)), INT(FLOOR(size(p, 1)*1.0D0/2 - 1)))/); 
        NeighbGridPoints = RESHAPE((/(((p(i), p(j), p(l), i=1, size(p, 1)), j=1, size(p, 1)), l=1, size(p, 1))/), (/size(p, 1)**3, 3/), order=(/2, 1/))

        ALLOCATE (InputValWT(ibDimT + 1), OutputValWT(ibDimT + 1))
        call dfftw_plan_dft_r2c_1d(cSpace%cGrid%cTransforms%iPlanTransform1D, iDimT, RealPressure, RealPressureC, fftw_estimate + fftw_unaligned)
        call dfftw_plan_dft_c2r_1d(cSpace%cGrid%cTransforms%iPlanTransform1D_inv, iDimT, RealPressureC, RealPressure, fftw_estimate + fftw_unaligned)
        call dfftw_plan_dft_r2c_1d(cSpace%cGrid%cTransforms%iPlanTransform1D_Own, 2*ibDimT, InputValWT, InputValWT, FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_1d(cSpace%cGrid%cTransforms%iPlanTransform1D_Own_inv, 2*ibDimT, OutputValWT, OutputValWT, FFTW_ESTIMATE)

        write (acTemp, "('Calculating Trilinear Interpolation') "); call PrintToLog(acTemp, 2)
        do iBubble = 1, ScattererParams(cSourceParams%iScCloud)%N
            ! Determine Beam Index of the point source
            ! In ParnacDataDef it is specified that cSpace%iStartT/X/Y/Z should be included
            ! The difference here is that the dBubbleLocationN parameter is in the normalized space
            ! So if we change this in real space then we should also change
            ! dBubbleLocation = (dBubbleLocationN + (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaNN
            ! So there is no difference if the starting position is not added in all the equations
            iBeamIndex(3) = nint(dBubbleLocationN(iBubble, 3)/cSpace%dDx)
            iBeamIndex(1) = nint(dBubbleLocationN(iBubble, 1)/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex(3))
            iBeamIndex(2) = nint(dBubbleLocationN(iBubble, 2)/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex(3))

            !Find the Index in the arBufferXYZ of the BubbleLocation
            XYZIndex = iBeamIndex(3)*iDimY*iDimX + iBeamIndex(2)*iDimX + iBeamIndex(1) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN
            if (XYZIndex >= 0 .AND. XYZIndex < cSpace%cGrid%iD0LocN) i_loc = i_loc + 1  ! Update i_loc if Bubble inside processor

            ! Determine current location in global grid
            ! Compare the global coordinates with the scatterer location
            dGlobXYZ(1) = (iBeamIndex(1) + iBeamOffsetX(cSpace, iBeamIndex(3)))*cSpace%dDx
            dGlobXYZ(2) = (iBeamIndex(2) + iBeamOffsetY(cSpace, iBeamIndex(3)))*cSpace%dDx
            dGlobXYZ(3) = (iBeamIndex(3))*cSpace%dDx

            dRadius(1) = sqrt(sum((dGlobXYZ - dBubbleLocationN(iBubble, :))**2))

            ! If the scatterer is not in the grid point (dRadius(1)>1e-10),
            ! then calculate pressure based on the 8 neighbouring points, trilinear interpolation
            if (dRadius(1) > 1.0D-10) then

                !each if the point is found in a processor different that XYZIndex location
                ! Xyz_factor is to find if the bubble real location is below or above the location of the grid point closer to the bubble
                xyz_factor = floor(dBubbleLocationN(iBubble, :) - dGlobXYZ(:)) + 1.0D0/2.0D0

                do i = 1, size(ones, 1)
                    iBeamIndexN = iBeamIndex + ones(i, :) + int(xyz_factor + 1.0D0/2.0D0)
                    XYZ_i1(i) = iBeamIndexN(3)*iDimY*iDimX + iBeamIndexN(2)*iDimX + iBeamIndexN(1) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN

                    if (XYZ_i1(i) >= 0 .AND. XYZ_i1(i) < cSpace%cGrid%iD0LocN) then

                        RealPressure_i1(i, :) = cSpace%cGrid%parD0(1 + XYZ_i1(i)*iDimT:(XYZ_i1(i) + 1)*iDimT)
                        GlobalPos(i, 1) = (cSpace%cGrid%aiD0Loc(XYZ_i1(i) + 1, 1) + iBeamOffsetX(cSpace, cSpace%cGrid%aiD0Loc(XYZ_i1(i) + 1, 3)))*cSpace%dDx
                        GlobalPos(i, 2) = (cSpace%cGrid%aiD0Loc(XYZ_i1(i) + 1, 2) + iBeamOffsetY(cSpace, cSpace%cGrid%aiD0Loc(XYZ_i1(i) + 1, 3)))*cSpace%dDx
                        GlobalPos(i, 3) = (cSpace%cGrid%aiD0Loc(XYZ_i1(i) + 1, 3))*cSpace%dDx
                        temp_globloc(1:3) = GlobalPos(i, :)
                        temp_globloc(4:iDimT + 3) = RealPressure_i1(i, :)

                        ! Bubbles' neighbours that are positioned in a different processor, are transferred via MPI_Send to the correct cpu
                        if (XYZIndex < 0 .OR. XYZIndex >= cSpace%cGrid%iD0LocN) then
                            BubbleProc = (iBeamIndex(1) + iBeamIndex(2)*iDimX + iBeamIndex(3)*iDimX*iDimY)/cSpace%cGrid%iD0LocN
                            call MPI_Send(temp_globloc, iDimT + 3, MPI_DOUBLE_PRECISION, BubbleProc, (iBubble - 1)*size(ones, 1) + i, MPI_COMM_WORLD, iErr)
                        else
                            LocPres(i_loc, 1 + (3 + iDimT)*(i - 1):(3 + iDimT)*i) = temp_globloc
                        end if
                    end if

                end do

                if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'around') then
                    do i = 1, size(NeighbGridPoints, 1)
                        iBeamIndexN = iBeamIndex + NeighbGridPoints(i, :) + int(xyz_factor + 1.0D0/2.0D0)
                        XYZ_i2(i) = iBeamIndexN(3)*iDimY*iDimX + iBeamIndexN(2)*iDimX + iBeamIndexN(1) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN

                        if (XYZ_i2(i) >= 0 .AND. XYZ_i2(i) < cSpace%cGrid%iD0LocN) then
                            count_neighbs = count_neighbs + 1
                            temploc(count_neighbs, 1) = iBubble
                            temploc(count_neighbs, 2) = i
                        end if

                    end do
                end if

                if (XYZIndex >= 0 .AND. XYZIndex < cSpace%cGrid%iD0LocN) then
                    ! Find the scatterer position and do inside if all the necessary actions because of spatial decomposition
                    ! Place the Bubble in respect to i_loc for later use for storing its values to a file
                    RealPressure = cSpace%cGrid%parD0(1 + XYZIndex*iDimT:(XYZIndex + 1)*iDimT); ! Array buffer with length size of iDimT for every XYZ position

                    do i = 1, size(ones, 1)
                        ! Save to temp_globloc which is a temporary buffer the position and pressure of the bubbles that are not positioned close to
                        ! each processor's grid.
                        RealPressure_i1(i, :) = 1D-30; GlobalPos(i, :) = 1D-30; temp_globloc = 1D-30; 
                        temp_globloc = LocPres(i_loc, 1 + (3 + iDimT)*(i - 1):(3 + iDimT)*i)

                        ! If some bubble is positioned close to the boundaries, receive the sent message
                        if (maxval(ABS(temp_globloc(4:iDimT + 3))) < 1.0D-20 .OR. maxval(abs(dBubbleLocationN(iBubble, :) - temp_globloc(1:3))) > size(trilin, 1)/2*cSpace%dDx) then
                            xyz_factor = floor(dBubbleLocationN(iBubble, :) - dGlobXYZ(:)) + 1.0D0/2.0D0; 
                            iBeamIndexN = iBeamIndex + ones(i, :) + int(xyz_factor + 1.0D0/2.0D0)
                            BubbleProc = (iBeamIndexN(3)*iDimY*iDimX + iBeamIndexN(2)*iDimX + iBeamIndexN(1))/cSpace%cGrid%iD0LocN
                            call MPI_Recv(temp_globloc, iDimT + 3, MPI_DOUBLE_PRECISION, BubbleProc, (iBubble - 1)*size(ones, 1) + i, MPI_COMM_WORLD, iStat, iErr) ! IMprove with MPI_Irecv
                        end if

                        ! Assign the values of the temp_globloc to the parameters that will be used for the interpolation
                        GlobalPos(i, :) = temp_globloc(1:3)
                        RealPressure_i1(i, :) = temp_globloc(4:iDimT + 3)

                        ! if Location is transferred in a correct way. If the difference with the bubble location is bigger than a space step,
                        ! then it means that the transfer did not succeed.
                        log_isnan = sum(isnan(GlobalPos(i, :))*-1) + sum(isnan(RealPressure_i1(i, :))*-1)
                        if (maxval(abs(dBubbleLocationN(iBubble, :) - GlobalPos(i, :))) > size(trilin, 1)/2*cSpace%dDx .OR. maxval(ABS(RealPressure_i1(i, :))) < 1E-20 .OR. log_isnan > 0) then
                            write (acTemp, "('Error in Trilinear Interpolation in Location , Bubble No. ', I7, ' Neighbouring ',I5, ' .')") iBubble, i
                            call PrintToLog(acTemp, 1)
                            write (*, *) dBubbleLocationN(iBubble, :)
                            write (*, *) GlobalPos(i, :)
                            write (*, *) "Exists : ", XYZIndex, cSpace%cGrid%iProcID, maxval(abs(dBubbleLocationN(iBubble, :) - GlobalPos(i, :))) > size(trilin, 1)/2*cSpace%dDx, maxval(ABS(RealPressure_i1(i, :))) < 1D-20, log_isnan
                            call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                        end if

                        RealPressure_Own(i, :) = 0.0D0; 
                        diff(i, :) = dBubbleLocationN(iBubble, :) - GlobalPos(i, :)

                        iTargetIndex(3) = nint(GlobalPos(i, 3)/cSpace%dDx)
                        iTargetIndex(1) = nint(GlobalPos(i, 1)/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex(3))
                        iTargetIndex(2) = nint(GlobalPos(i, 2)/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex(3))

                        if (MAXVAL(ABS(tempcon(1 + (i_loc - 1)*ibDimT:i_loc*ibDimT))) > 0.0D0) &
                            call GREENS1DFREQ(cSpace, tempcon(1 + (i_loc - 1)*ibDimT:i_loc*ibDimT), RealPressure_Own(i, iDimT - ibDimT + 1:iDimT), iTargetIndex, diff(i, :))
                        RealPressure_Own(i, :) = RealPressure_Own(i, :)*(cModelParams%freq0/cMediumParams%c0)*(dK_c/pi)**3

                        ! if ( MAXVAL(ABS(dBubbleLocationN(iBubble,:) - &
                        ! ( (/-0.0366470588235294 ,   -0.0366470588235294 ,    5.12194117647059/)/dLambdaMM  - (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ /)* cSpace%dDx))) < 1E-7) then! At a gridpoint  @ 1.7E6 freq0 , 6 Fnyq
                        ! write(*,*) "==================================="
                        ! write(*,*) iTargetIndex - iBeamIndex
                        ! write(*,*) RealPressure_Own(i,:)
                        ! write(*,*) "==================================="
                        ! write(*,*) iTargetIndex - iBeamIndex
                        ! write(*,*) RealPressure_i1(i,:)
                        ! endif
                        RealPressure_i1(i, :) = RealPressure_i1(i, :) - RealPressure_Own(i, :)
                        ! This is for the non-comoving time window. There should be a time shift
                        RealPressure_i1(i, 1:iDimT - (ones(i, 3) + size(trilin, 1)/2)) = RealPressure_i1(i, 1 + ones(i, 3) + size(trilin, 1)/2:iDimT)

                    end do
                    ! Do the 3D interpolation in order to acquire the pressure , RealPressureFinal
                    call INTERP3D(GlobalPos, RealPressure_i1, dBubbleLocationN(iBubble, :), Interp_PressureBL)
                    RealPressure = Interp_PressureBL

                    call dfftw_execute(cSpace%cGrid%cTransforms%iPlanTransform1D, RealPressure, RealPressureC)
                    RealPressureC = RealPressureC*exp(-im*dMultFactor*diff(1, 3))/real(iDimT, dp)*dTaperingWindow(iDimT, 1.0_dp/(2*iDimT*cSpace%dDt), 0.0D0, 0.1D0)
                    call dfftw_execute(cSpace%cGrid%cTransforms%iPlanTransform1D_inv, RealPressureC, RealPressure)

                    tempcon(1 + (i_loc - 1)*ibDimT:i_loc*ibDimT) = RealPressure(iDimT - ibDimT + 1:iDimT); 
                    BubbleID(i_loc) = iBubble
                end if
                !   else

                !         if (XYZIndex>=0 .AND. XYZIndex<cSpace%cGrid%iD0LocN) then
                !                 i=1;
                !                 RealPressure_Own(i,:) = 0.0D0;
                !                 diff(i,:) = dBubbleLocationN(iBubble,:) - dGlobXYZ

                !                 if (MAXVAL(ABS(tempcon(1+(i_loc-1)*ibDimT:i_loc*ibDimT)))>0.0D0) call GREENS1DFREQ(cSpace,tempcon(1+(i_loc-1)*ibDimT:i_loc*ibDimT), RealPressure_Own(i,iDimT - ibDimT + 1: iDimT), iTargetIndex, diff(i,:))
                !                 RealPressure_Own(i,:) = RealPressure_Own(i,:) * (cModelParams%freq0/cMediumParams%c0) * (dK_c/pi)**3

                !                 RealPressure = cSpace%cGrid%parD0(1+XYZIndex*iDimT: (XYZIndex+1)*iDimT);   ! Array buffer with length size of iDimT for every XYZ position
                !                 RealPressure  = RealPressure - RealPressure_Own(i,:)

                !                 tempcon(1+(i_loc-1)*ibDimT:i_loc*ibDimT) = RealPressure(iDimT - ibDimT + 1: iDimT);
                !                 BubbleID(i_loc) = iBubble
                !          endif

            end if
        end do
        call MPI_BARRIER(MPI_COMM_WORLD, iErr); 
        !************************************* 3D INTEPR  **************************************************
        call dfftw_destroy_plan(cSpace%cGrid%cTransforms%iPlanTransform1D)
        call dfftw_destroy_plan(cSpace%cGrid%cTransforms%iPlanTransform1D_inv)
        call dfftw_destroy_plan(cSpace%cGrid%cTransforms%iPlanTransform1D_Own)
        call dfftw_destroy_plan(cSpace%cGrid%cTransforms%iPlanTransform1D_Own_inv)

        ! Deallocate to save memory
        if (ALLOCATED(LocPres)) then
            iMemAllocated = iMemAllocated - PRODUCT(SHAPE(LocPres))*dpS
            DEALLOCATE (LocPres)
        end if
        iMemAllocated = iMemAllocated - PRODUCT(SHAPE(temp_globloc))*dpS
        DEALLOCATE (temp_globloc)

    END SUBROUTINE TrilinInterp_SNDRCV

    SUBROUTINE CalculateTimeSignature(cSpace, dBubbleContrast, ibDimT)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 201128
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function CalculateCloudPressure calculates the pressure of the total bubble
        !        cloud.  This is done in two phases. The first phase is for each iteration ,except
        !   the last one. For this phase, all the neighbouring point of microbubbles are used
        !        and then the pressure is calculated in these points by all the microbubbles.
        !        In order to keep memory low , this is done by each cpu using the neighbouring points
        !        that are located inside its subdomain. The second phase is the last iteration.
        !        At this iteration ,the pressure is calculated at each point of the domain.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                       io   type(space)  The space in which the point has to be given
        !   dBubbleLocationN         r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   dBubbleContrast                 r   dp   a vector with the time signatures of the microbubbles
        !   NeighbGPIndex                 i   i8b  a vector with the location of the neighbouring points of microbubbles in each cpu domain
        !   i_neighbgp                          i   i8b  number of neighbouring points ( Do not get confused, not all elements of NeighbGPIndex
        !                                                                          have the values , only 1:i_neighbgp)
        !   last_iter                         l        logical value which states if it is the last iteration or not
        !
        ! *****************************************************************************
        type(Space), intent(inout), target      ::      cSpace
        real(dp), intent(inout)                 ::      dBubbleContrast(:)
        integer, intent(inout)                  ::      ibDimT

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   P_Scattered                       dpc   Scattered Pressure Matrix for MKL
        !   dSincMult                             dpc   a matrix that contains spatially filtered dirac
        !   NeighbLocInd                        i8b   a vector which contains the indexes of the neighbouring points
        !   dGlobBXYZ                        r     a matrix of the gridpoints location of each cpu domain (same dim as dBubbleLocationN)
        !   M, ibDimT                               i     integer value for MKL , Time length
        !   N, GridPointsNum            i          integer value for MKL , Space length
        !   K, BLayer                                i           integer value for MKL , Bubbles per layer
        !   iIndex, i, iBubble            i8b   index variables
        !   B_Div                                           i8b   divisor of total numbers of variables
        !        dLambdaMM                                dp          value of wavelength , normalization factor
        !        dK_c                                        dp          Kutoff frequency
        !
        ! *****************************************************************************

        type(Grid), pointer                     ::      pcGrid

        integer                                 ::      iErr, iStat(MPI_STATUS_SIZE), ErrCode

        integer(i8b)                            ::      iBubble
        integer                                 ::      BubbleProc

        integer(i8b), parameter                 ::      n_pad = 2, OS_Factor = 1

        real(dp)                                ::      RealPressure(OS_Factor*ibDimT), RealTime(ibDimT), LinearAmplitude(ibDimT)
        real(dp)                                ::      RealPressurePad(ibDimT*n_pad), RealTimePad(ibDimT*n_pad), RealTimePadOut(ibDimT*n_pad)
        real(dp)                                ::      V_dd_norm(ibDimT), V_dd_Pad(ibDimT*n_pad), Time_Signature(ibDimT)

        integer                                 ::      Bubble_count(cSpace%cGrid%iProcN), Bubble_condisps(cSpace%cGrid%iProcN)

        real(dp), allocatable                   ::      tempcon(:)

        integer(i8b)                            ::      iDimT, i, Bubble_per_CPU
        character(len=1024)                     ::      actemp
        type(Space)                             ::      cSpaceScatterer

        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   FFTW
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   First translate the NeighbLocInd into normalized coordinates in (x,y,z) if not in last iteration
        !01          and set the points for calculation only to neighbouring, otherwise calculate for every point of domain
        !02          Set a good choice for the number of Bubbles per Layer ( should check memory usage)
        !03          Iterate for the number of Layers (Total Bubbles / No. Bubbles per Layer)
        !04          Calculate dirac function matrix
        !05          Do the matrix multiplication based on the space points for calculation
        !06          if not in last iteration , transfer the calculated values to a matrix for easy replacement to cSpace var
        !07          else just replace the calculated values
        !
        ! =============================================================================

        pcGrid => cSpace%cGrid

        iDimT = cSpace%iDimT   !t dimensions

        V_dd_norm = 0; 
        V_dd_Pad = 0; 
        LinearAmplitude = 0; 
        Bubble_per_CPU = ScattererParams(cSourceParams%iScCloud)%N/pcGrid%iProcN + (pcGrid%iProcID + 1)/pcGrid%iProcN*mod(ScattererParams(cSourceParams%iScCloud)%N, pcGrid%iProcN)

        ALLOCATE (tempcon(Bubble_per_CPU*ibDimT))

        RealTime = [(i + iDimT - ibDimT, i=cSpace%iStartT, ibDimT - 1 + cSpace%iStartT)]*cSpace%dDt/cModelParams%freq0     ! Array with the real time unnormalized
        RealTimePad = 1D-30
        call linspace(RealTimePad, RealTime(1), RealTime(size(RealTime, 1)), ibDimT*n_pad)

        if (trim(ScattererParams(cSourceParams%iScCloud)%ScType) == 'microbubble') then
            write (acTemp, "('Calculate the time signatures of the microbubbles')"); call PrintToLog(acTemp, 2)
            do i = 1, Bubble_per_CPU
                iBubble = i + ScattererParams(cSourceParams%iScCloud)%N/pcGrid%iProcN*pcGrid%iProcID
                RealPressure = dBubbleContrast(1 + (iBubble - 1)*ibDimT:iBubble*ibDimT)

                ! You can use INTERP1DFREQ (frequency interpolation) instead, but this is more accurate and efficient
                ! call INTERP1DFREQ(RealPressure,RealPressurePad, cSpace,0)
                call INTERP1D(RealTime, RealPressure, RealTimePad, RealPressurePad)

                ! Marmottant solver ( ODE Solver) , Result is V_dd_Pad (volume acceleration) and RealTimePadOut which is the updated time values
                call RP_SOLVER(RealPressurePad, RealTimePad, ibDimT*n_pad, iBubble, V_dd_Pad, RealTimePadOut)

                ! Return to the initial dimensions!
                ! call INTERP1DFREQ(V_dd_Pad, V_dd_norm, cSpace,0)
                call INTERP1D(RealTimePadOut, V_dd_Pad, RealTime, V_dd_norm)
                Time_Signature = cMediumParams%rho0*V_dd_norm; 
                ! Check to see if results are realistic
                !*************************** CHANGED THIS with addition of V_dd_pad
                ! ! if ( .NOT. ALL(ABS(V_dd_norm)<1E100_dp) .OR. .NOT. ALL(ABS(V_dd_pad)<1E100_dp)) then ! Check Nan , Large values or Infinity
                ! ! write(acTemp,"('Error in ODE SOLVER in V_dd_pad , Bubble No. ', I7, ' .')") iBubble
                ! call PrintToLog(acTemp,1)
                ! open (7, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('output_file_pressure')//int2str(cModelParams%iIter)//int2str(pcGrid%iProcID)//'.txt'),status = 'UNKNOWN', action='write',position='append')
                ! write(7,*) iBubble
                ! write(7,*) RealTime, RealPressure,RealTimePad,RealPressurePad
                ! close(7)
                ! open (8, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('v_dd_norm')//int2str(cModelParams%iIter)//int2str(pcGrid%iProcID)//'.txt'),status = 'UNKNOWN', action='write',position='append')
                ! write(8,*) RealTime,V_dd_norm,RealTimePadOut,V_dd_pad
                ! close(8)
                ! call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                ! endif
                if (.NOT. ALL(ABS(Time_Signature) < 1E100_dp)) call PrintToLog("Error: dBubbleContrast in CalculateTimeSignature!!!", 1)

                tempcon(1 + (i - 1)*ibDimT:i*ibDimT) = Time_Signature; 
            end do
        elseif (trim(ScattererParams(cSourceParams%iScCloud)%ScType) == 'linear') then
            write (acTemp, "('Calculate the time signatures of the linear scatterers')"); call PrintToLog(acTemp, 2)

            call InitSpace(cSpaceScatterer, iSI_FIELD, cModelParams%UseYSymmetry, &
                           int(ibDimT, i8b), 1, 1, 1, &
                           0_i8b, cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ, 0_i8b, 0_i8b, 0_i8b, &
                           cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
            !iDimT is such that each processor gets one complex, i.e. two real slices..
            call InitGrid(cSpaceScatterer, 1, 1, (/.false., .false., .false., .false./)); 
            call GridDistr0CreateEmpty(cSpaceScatterer%cGrid); 
            do i = 1, Bubble_per_CPU
                iBubble = i + ScattererParams(cSourceParams%iScCloud)%N/pcGrid%iProcN*pcGrid%iProcID
                RealPressure = dBubbleContrast(1 + (iBubble - 1)*ibDimT:iBubble*ibDimT)

                LinearAmplitude = -4.0/3.0*pi*(ScattererParams(cSourceParams%iScCloud)%R0(iBubble))**3*cMediumParams%rho0* &
                                  (1.0D0/(ScattererParams(cSourceParams%iScCloud)%rho1*ScattererParams(cSourceParams%iScCloud)%c1**2) - 1.0D0/(cMediumParams%rho0*cMediumParams%c0**2))
                call LinScatterer(cSpace, RealPressure, Time_Signature, LinearAmplitude*cModelParams%freq0**2, 2)
                cSpaceScatterer%cGrid%parD0 = Time_Signature; 
                ! cSpaceScatterer%cGrid%parD0  = RealPressure;
                ! call NonLinScatterer(cSpaceScatterer, LinearAmplitude* cModelParams%freq0**2, 2, 1 )

                ! Check to see if results are realistic
                !*************************** CHANGED THIS with addition of V_dd_pad
                ! if ( .NOT. ALL(ABS(V_dd_norm)<1E100_dp) .OR. .NOT. ALL(ABS(V_dd_pad)<1E100_dp)) then ! Check Nan , Large values or Infinity
                ! write(acTemp,"('Error in ODE SOLVER in V_dd_pad , Bubble No. ', I7, ' .')") iBubble
                ! call PrintToLog(acTemp,1)
                ! open (7, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('output_file_pressure')//int2str(cModelParams%iIter)//int2str(pcGrid%iProcID)//'.txt'),status = 'UNKNOWN', action='write',position='append')
                ! write(7,*) iBubble
                ! write(7,*) RealTime, RealPressure,RealTimePad,RealPressurePad
                ! close(7)
                ! open (8, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//trim('v_dd_norm')//int2str(cModelParams%iIter)//int2str(pcGrid%iProcID)//'.txt'),status = 'UNKNOWN', action='write',position='append')
                ! write(8,*) RealTime,Time_Signature,RealTimePadOut,V_dd_pad
                ! close(8)
                ! call MPI_Abort(MPI_COMM_WORLD, ErrCode, iErr)
                ! endif
                if (.NOT. ALL(ABS(Time_Signature) < 1E100_dp)) call PrintToLog("Error: dBubbleContrast in CalculateTimeSignature!!!", 1)

                tempcon(1 + (i - 1)*ibDimT:i*ibDimT) = cSpaceScatterer%cGrid%parD0; 
            end do
            call DestructSpace(cSpaceScatterer)
        end if
        Bubble_count = 0; Bubble_condisps = 0; 
        do i = 1, pcGrid%iProcN
            Bubble_count(i) = ScattererParams(cSourceParams%iScCloud)%N/pcGrid%iProcN + i/pcGrid%iProcN*mod(ScattererParams(cSourceParams%iScCloud)%N, pcGrid%iProcN)
            Bubble_condisps(i) = sum(ibDimT*Bubble_count(1:i)) - ibDimT*Bubble_count(i)        ! The same for contrast source term , multiplied by length of T
        end do

        dBubbleContrast = 1D-30
        call MPI_Allgatherv(tempcon, ibDimT*Bubble_per_CPU, MPI_DOUBLE_PRECISION, dBubbleContrast, ibDimT*Bubble_count, Bubble_condisps, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, iErr)
        call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values
        DEALLOCATE (tempcon)

    END SUBROUTINE CalculateTimeSignature

    SUBROUTINE CalculateCloudPressure_Around(cSpace, dBubbleLocationN, dBubbleContrast, BubbleID, temploc, NeighbGridPoints, count_neighbs)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 201128
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function CalculateCloudPressure calculates the pressure of the total bubble
        !        cloud.  This is done in two phases. The first phase is for each iteration ,except
        !   the last one. For this phase, all the neighbouring point of microbubbles are used
        !        and then the pressure is calculated in these points by all the microbubbles.
        !        In order to keep memory low , this is done by each cpu using the neighbouring points
        !        that are located inside its subdomain. The second phase is the last iteration.
        !        At this iteration ,the pressure is calculated at each point of the domain.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                       io   type(space)  The space in which the point has to be given
        !   dBubbleLocationN         r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   dBubbleContrast                 r   dp   a vector with the time signatures of the microbubbles
        !   NeighbGPIndex                 i   i8b  a vector with the location of the neighbouring points of microbubbles in each cpu domain
        !   i_neighbgp                          i   i8b  number of neighbouring points ( Do not get confused, not all elements of NeighbGPIndex
        !                                                                          have the values , only 1:i_neighbgp)
        !   last_iter                         l        logical value which states if it is the last iteration or not
        !
        ! *****************************************************************************
        real(dp), intent(inout)             ::          dBubbleLocationN(:, :)
        real(dp), intent(in)                ::          dBubbleContrast(:)
        integer(i8b), intent(in)            ::          BubbleID(:), NeighbGridPoints(:, :), count_neighbs, temploc(:, :)
        type(Space), intent(inout)          ::          cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   P_Scattered                       dpc   Scattered Pressure Matrix for MKL
        !   dSincMult                             dpc   a matrix that contains spatially filtered dirac
        !   NeighbLocInd                        i8b   a vector which contains the indexes of the neighbouring points
        !   dGlobBXYZ                        r     a matrix of the gridpoints location of each cpu domain (same dim as dBubbleLocationN)
        !   M, ibDimT                               i     integer value for MKL , Time length
        !   N, GridPointsNum            i          integer value for MKL , Space length
        !   K, BLayer                                i           integer value for MKL , Bubbles per layer
        !   iIndex, i, iBubble            i8b   index variables
        !   B_Div                                           i8b   divisor of total numbers of variables
        !        dLambdaMM                                dp          value of wavelength , normalization factor
        !        dK_c                                        dp          Kutoff frequency
        !
        ! *****************************************************************************

        integer                             ::          i, j, iBubble, iDimX, iDimY, iDimZ, iDimT
        integer(i8b)                        ::          iBeamIndex(3), iBeamIndexN(3), XYZ_i2
        real(dp)                            ::          dLambdaMM, dK_c, dGlobXYZ(3), xyz_factor(3), diff(8, 3), GlobalPos(1, 3), Dirac_factor
        character(len=1024)                 ::          actemp

        write (acTemp, "('Calculate the contrast source term around each scatterer')"); call PrintToLog(acTemp, 3)

        iDimX = cSpace%iDimX
        iDimY = cSpace%iDimY
        iDimZ = cSpace%iDimZ
        iDimT = cSpace%iDimT

        dK_c = pi/cSpace%dDx

        do j = 1, count_neighbs
            ! Determine Beam Index of the point source
            iBubble = MINLOC(BubbleID, 1, BubbleID == temploc(j, 1))
            i = temploc(j, 2)

            iBeamIndex(3) = nint(dBubbleLocationN(iBubble, 3)/cSpace%dDx)
            iBeamIndex(1) = nint(dBubbleLocationN(iBubble, 1)/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex(3))
            iBeamIndex(2) = nint(dBubbleLocationN(iBubble, 2)/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex(3))

            dGlobXYZ(1) = (iBeamIndex(1) + iBeamOffsetX(cSpace, iBeamIndex(3)))*cSpace%dDx
            dGlobXYZ(2) = (iBeamIndex(2) + iBeamOffsetY(cSpace, iBeamIndex(3)))*cSpace%dDx
            dGlobXYZ(3) = (iBeamIndex(3))*cSpace%dDx

            xyz_factor = floor(dBubbleLocationN(iBubble, :) - dGlobXYZ) + 1.0D0/2.0D0; 
            iBeamIndexN = iBeamIndex + NeighbGridPoints(i, :) + int(xyz_factor + 1.0D0/2.0D0)
            XYZ_i2 = iBeamIndexN(3)*iDimY*iDimX + iBeamIndexN(2)*iDimX + iBeamIndexN(1) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN
            ! if (XYZ_i2(i)<0 .OR. XYZ_i2(i)>=cSpace%cGrid%iD0LocN  ) write(*,*) "ERROR", iBubble,i

            GlobalPos(1, 1) = (cSpace%cGrid%aiD0Loc(XYZ_i2 + 1, 1) + iBeamOffsetX(cSpace, cSpace%cGrid%aiD0Loc(XYZ_i2 + 1, 3)))*cSpace%dDx
            GlobalPos(1, 2) = (cSpace%cGrid%aiD0Loc(XYZ_i2 + 1, 2) + iBeamOffsetY(cSpace, cSpace%cGrid%aiD0Loc(XYZ_i2 + 1, 3)))*cSpace%dDx
            GlobalPos(1, 3) = (cSpace%cGrid%aiD0Loc(XYZ_i2 + 1, 3))*cSpace%dDx

            diff(1, :) = dBubbleLocationN(iBubble, :) - GlobalPos(1, :)
            Dirac_factor = PRODUCT(dSinc(dK_c*diff(1, :)))

            cSpace%cGrid%parD0(1 + XYZ_i2*iDimT:XYZ_i2*iDimT + iDimT) = &
                cSpace%cGrid%parD0(1 + XYZ_i2*iDimT:XYZ_i2*iDimT + iDimT) + dBubbleContrast(1 + (iBubble - 1)*iDimT:iBubble*iDimT)*Dirac_factor
        end do

        cSpace%cGrid%parD0 = cSpace%cGrid%parD0*(cModelParams%freq0/cMediumParams%c0)*(dK_c/pi)**3
    END SUBROUTINE CalculateCloudPressure_Around

    SUBROUTINE CalculateCloudPressure(cSpace, dBubbleLocationN, dBubbleContrast)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 201128
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function CalculateCloudPressure calculates the pressure of the total bubble
        !        cloud.  This is done in two phases. The first phase is for each iteration ,except
        !   the last one. For this phase, all the neighbouring point of microbubbles are used
        !        and then the pressure is calculated in these points by all the microbubbles.
        !        In order to keep memory low , this is done by each cpu using the neighbouring points
        !        that are located inside its subdomain. The second phase is the last iteration.
        !        At this iteration ,the pressure is calculated at each point of the domain.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                       io   type(space)  The space in which the point has to be given
        !   dBubbleLocationN         r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   dBubbleContrast                 r   dp   a vector with the time signatures of the microbubbles
        !   NeighbGPIndex                 i   i8b  a vector with the location of the neighbouring points of microbubbles in each cpu domain
        !   i_neighbgp                          i   i8b  number of neighbouring points ( Do not get confused, not all elements of NeighbGPIndex
        !                                                                          have the values , only 1:i_neighbgp)
        !   last_iter                         l        logical value which states if it is the last iteration or not
        !
        ! *****************************************************************************
        real(dp), intent(inout)             ::          dBubbleLocationN(:, :)
        real(dp), intent(in)                ::          dBubbleContrast(:)
        type(Space), intent(inout)          ::          cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   P_Scattered                       dpc   Scattered Pressure Matrix for MKL
        !   dSincMult                             dpc   a matrix that contains spatially filtered dirac
        !   NeighbLocInd                        i8b   a vector which contains the indexes of the neighbouring points
        !   dGlobBXYZ                        r     a matrix of the gridpoints location of each cpu domain (same dim as dBubbleLocationN)
        !   M, ibDimT                               i     integer value for MKL , Time length
        !   N, GridPointsNum            i          integer value for MKL , Space length
        !   K, BLayer                                i           integer value for MKL , Bubbles per layer
        !   iIndex, i, iBubble            i8b   index variables
        !   B_Div                                           i8b   divisor of total numbers of variables
        !        dLambdaMM                                dp          value of wavelength , normalization factor
        !        dK_c                                        dp          Kutoff frequency
        !
        ! *****************************************************************************

        real(dp), allocatable               ::          P_scattered(:, :), P_scattered_Total(:, :), dSincMult(:, :)
        real(dp), allocatable               ::          dGlobBXYZ(:, :)
        integer(i8b), allocatable           ::          NeighbLocInd(:), XYZIndex(:), NeighbGPIndex(:)
        integer, allocatable                ::          GridPoint_OnProc(:, :), TimePoints(:)

        integer                             ::          M, N, K, iDimT, BLayer, iDimX, iDimY, iDimZ, iDimXYZ, i, j, BubbleProc, byte_precision
        integer                             ::          GridPoint_Count(cSpace%cGrid%iProcN), GridPoint_Disps(cSpace%cGrid%iProcN)
        integer                             ::          iErr, iStat(MPI_STATUS_SIZE), iREQUEST, MPI_NEWTYPE, ErrCode
        integer(i8b)                        ::          iIndex, iLI, x, y, z, iBubble, B_Div, GridPointsNum, i_count, iBeamIndex_min(3), iBeamIndex_max(3), i_neighbgp
        integer(i8b)                        ::          TotalCloudN
        real(dp)                            ::          dLambdaMM, dK_c

        character(len=1024)                 ::          actemp

        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   FFTW
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   First translate the NeighbLocInd into normalized coordinates in (x,y,z) if not in last iteration
        !01          and set the points for calculation only to neighbouring, otherwise calculate for every point of domain
        !02          Set a good choice for the number of Bubbles per Layer ( should check memory usage)
        !03          Iterate for the number of Layers (Total Bubbles / No. Bubbles per Layer)
        !04          Calculate dirac function matrix
        !05          Do the matrix multiplication based on the space points for calculation
        !06          if not in last iteration , transfer the calculated values to a matrix for easy replacement to cSpace var
        !07          else just replace the calculated values
        !
        ! =============================================================================

        iDimT = cSpace%iDimT
        iDimX = cSpace%iDimX
        iDimY = cSpace%iDimY
        iDimZ = cSpace%iDimZ

        dLambdaMM = (cMediumParams.c0*1.0D3)/cModelParams%freq0; ! Normalization factor λ in [mm]
        dK_c = pi/cSpace%dDx
        TotalCloudN = size(dBubbleLocationN, 1)
        ! At this point, the normalized location of the gridpoints of interest are calculated
        if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'all' .OR. cModelParams%iIter == cModelParams%Numiterations(1) - 1) then
            ! In this case, every gridpoint's location is computed
            write (acTemp, "('Calculate the contrast source term at each spatial gridpoint')"); call PrintToLog(acTemp, 3)
            GridPointsNum = cSpace%cGrid%iD0LocN
            ALLOCATE (dGlobBXYZ(GridPointsNum, 3))
            ! This is the way to compute them and in this order. This is important because later they will be used
            ! for the computation of the total scattered pressure and they should be placed in this order
            do i = 0, GridPointsNum - 1
                dGlobBXYZ(1 + i, 1) = (cSpace%cGrid%aiD0Loc(i + 1, 1) + iBeamOffsetX(cSpace, cSpace%cGrid%aiD0Loc(i + 1, 3)))*cSpace%dDx
                dGlobBXYZ(1 + i, 2) = (cSpace%cGrid%aiD0Loc(i + 1, 2) + iBeamOffsetY(cSpace, cSpace%cGrid%aiD0Loc(i + 1, 3)))*cSpace%dDx
                dGlobBXYZ(1 + i, 3) = (cSpace%cGrid%aiD0Loc(i + 1, 3))*cSpace%dDx
            end do

            BLayer = 1000; 
        else

            j = 6
            iBeamIndex_min(3) = nint(minval(dBubbleLocationN(:, 3))/cSpace%dDx)
            iBeamIndex_min(1) = nint(minval(dBubbleLocationN(:, 1))/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex_min(3))
            iBeamIndex_min(2) = nint(minval(dBubbleLocationN(:, 2))/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex_min(3))
            iBeamIndex_min = iBeamIndex_min - j

            iBeamIndex_max(3) = nint(maxval(dBubbleLocationN(:, 3))/cSpace%dDx)
            iBeamIndex_max(1) = nint(maxval(dBubbleLocationN(:, 1))/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex_max(3))
            iBeamIndex_max(2) = nint(maxval(dBubbleLocationN(:, 2))/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex_max(3))
            iBeamIndex_max = iBeamIndex_max + j

            i_neighbgp = (iBeamIndex_max(1) - iBeamIndex_min(1) + 1)*(iBeamIndex_max(2) - iBeamIndex_min(2) + 1)*(iBeamIndex_max(3) - iBeamIndex_min(3) + 1)

            iMemAllocated = iMemAllocated + i_neighbgp*i4bs
            ALLOCATE (NeighbGPIndex(i_neighbgp))
            NeighbGPIndex(1:i_neighbgp) = (/(((z*iDimY*iDimX + y*iDimX + x, x=iBeamIndex_min(1), iBeamIndex_max(1)), y=iBeamIndex_min(2), iBeamIndex_max(2)), z=iBeamIndex_min(3), iBeamIndex_max(3))/)

            ! In this case, only the gridpoints that are added in the main module above, are computed
            ! The number of the neighbouring points for every scatterer equals to size(p,1)**3
            write (acTemp, "('Calclulate the contrast source term in the gridpoints inside the boundaries of the cloud')"); call PrintToLog(acTemp, 3)
            GridPointsNum = i_neighbgp/cSpace%cGrid%iProcN + (cSpace%cGrid%iProcID + 1)/cSpace%cGrid%iProcN*mod(i_neighbgp, cSpace%cGrid%iProcN)! Do this for all processors
            ALLOCATE (dGlobBXYZ(GridPointsNum, 3))
            ! Because all cpus see all the gridpoints, we have to use a way to separate those
            ! that are located in a specific processor. Generate a matrix, where first column is the global value of each gridpoint
            ! and the 2nd value is the processor, where it is located.
            ALLOCATE (GridPoint_OnProc(2, GridPointsNum))

            do i = 0, GridPointsNum - 1
                ! Divide the number of gridpoints based on the number of cpus
                j = NeighbGPIndex((i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID + i + 1); 
                GridPoint_OnProc(1, i + 1) = i
                GridPoint_OnProc(2, i + 1) = j/cSpace%cGrid%iD0LocN

                ! This is the way to compute the normalized global values of each gridpoint  in space.
                dGlobBXYZ(1 + i, 1) = mod(j, cSpace%iDimX)*cSpace%dDx; ; 
                j = (j - mod(j, cSpace%iDimX))/cSpace%iDimX; 
                dGlobBXYZ(1 + i, 2) = mod(j, cSpace%iDimY)*cSpace%dDx; 
                j = (j - mod(j, cSpace%iDimY))/cSpace%iDimY; 
                dGlobBXYZ(1 + i, 3) = j*cSpace%dDx; 
            end do

            ! Increase the number of scatterers per layer, because there are less points than final iteration
            BLayer = 5E3; 
        end if
        M = PRODUCT(SHAPE(dBubbleContrast))/TotalCloudN; N = GridPointsNum; K = BLayer

        ! Because of big size arrays for clusters higher than 10^5 bubbles
        ! due to memory , a big array can not be created ,so I implement the idea
        ! of layering in the matrix, So I do a loop for every K bubbles and calculate
        ! the scattered pressure. I do this until all the bubbles hve been included
        if (TotalCloudN <= K .OR. MOD(TotalCloudN, K) .NE. 0) K = TotalCloudN
        B_Div = TotalCloudN/K; 
        ALLOCATE (P_scattered(M, N))
        ALLOCATE (P_scattered_Total(M + 1, N))
        ALLOCATE (dSincMult(K, N))

        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(dSincMult)))*sizeof(dSincMult(1, 1)) + PRODUCT(SHAPE(dGlobBXYZ))*dpS
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(P_scattered)) + PRODUCT(SHAPE(P_scattered_Total)))*sizeof(P_scattered(1, 1))

        write (acTemp, "('Create ', I3,' Layer(s) of ', I<INT(log10(K*1.0D0))+1>, ' Bubble(s) for Dirac Array')") B_Div, K; call PrintToLog(acTemp, 3)

        ! Here we generate only the contrast source term of interest! This array will include the nonlinear contrast source term for a
        ! number of time steps equal to (cSpace%iDimT+1)/iProcN in every processor. This means that instead of dividing the domain based on
        ! the number of cores, the time is divided. In this way we have the same number of computations per processor and every processor
        ! takes part in the computations
        P_scattered_Total = 1.0D-30; 
        do iBubble = 1, B_Div ! This is the loop for the number of layers of the bubbles
            ! Calculate the spatial signature in the point source definition for every scatterer
            do i = 0, N - 1
                dSincMult(:, 1 + i) = dSinc(dK_c*(dBubbleLocationN(K*(iBubble - 1) + 1:K*iBubble, 1) - dGlobBXYZ(1 + i, 1))) &
                                      *dSinc(dK_c*(dBubbleLocationN(K*(iBubble - 1) + 1:K*iBubble, 2) - dGlobBXYZ(1 + i, 2))) &
                                      *dSinc(dK_c*(dBubbleLocationN(K*(iBubble - 1) + 1:K*iBubble, 3) - dGlobBXYZ(1 + i, 3)))
            end do
            ! DO the multiplication with the temporal factor , in order to calculte the scattered pressure
            ! Each processor computes the scattered pressure in the gridpoints of interest only for the iterated layer of microbubbles
            ! Uncoment for double precision! You should change the type of P_scattered and P_scattered_Total too to dp
            ! CALL DGEMM('N','N',M,N,K,1.0D0,dBubbleContrast(:,1+K*(iBubble-1):K*iBubble) * (dK_c/pi)**3, M, dSincMult, K, 0.0D0, P_scattered, M)
            CALL DGEMM('N', 'N', M, N, K, 1.0D0, RESHAPE(dBubbleContrast(1 + K*(iBubble - 1)*M:K*iBubble*M)*(dK_c/pi)**3, (/M, K/)), M, dSincMult, K, 0.0D0, P_scattered, M)
            ! CALL SGEMM('N','N', M, N, K, 1.0, REAL(dBubbleContrast(:,1+K*(iBubble-1):K*iBubble) * (dK_c/pi)**3,sp), M, dSincMult, K, 0.0, P_scattered, M)
            P_scattered_Total(2:M + 1, :) = P_scattered_Total(2:M + 1, :) + P_Scattered
        end do
        if (.NOT. ALL(ABS(P_scattered_Total) < 1E100_dp)) call PrintToLog("Error: P_scattered !!!", 1)

        iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(dSincMult)))*sizeof(dSincMult(1, 1)) - PRODUCT(SHAPE(dGlobBXYZ))*dpS
        iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(P_scattered)))*sizeof(P_scattered(1, 1))
        DEALLOCATE (P_Scattered)
        DEALLOCATE (dSincMult)
        DEALLOCATE (dGlobBXYZ)

        write (acTemp, "('Matrix Multiplication has finished.')"); call PrintToLog(acTemp, 3)
        write (acTemp, "('GridPoints Distribution to each processor starts.')"); call PrintToLog(acTemp, 3)
        if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'all' .OR. cModelParams%iIter == cModelParams%Numiterations(1) - 1) then
            ALLOCATE (P_scattered(iDimT, N))
            P_scattered = 0.0D0; 
            P_scattered(iDimT - M + 1:iDimT, :) = P_scattered_Total(2:M + 1, :)
            cSpace%cGrid%parD0 = cSpace%cGrid%parD0 + RESHAPE(P_scattered*(cModelParams%freq0/cMediumParams%c0), (/iDimT*N/))
            DEALLOCATE (P_scattered)
        else
            ! We need this row in order to know for which value of the gridpoint, the pressure of each row(gridpoint) is computed so we can send it to the right processor
            P_scattered_Total(1, :) = NeighbGPIndex(1 + (i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID:(i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID + N)*1.0D0

            ! Here, the array for the number of elements GridPoint_Count and
            ! the array for the starting index of the location of the elements per processor GridPoint_Disps is generated.
            GridPoint_Count = 0; 
            GridPoint_Disps = 0; 
            do BubbleProc = 0, cSpace%cGrid%iProcN - 1
                j = i_neighbgp/cSpace%cGrid%iProcN
                i = i_neighbgp/cSpace%cGrid%iProcN + (BubbleProc + 1)/cSpace%cGrid%iProcN*mod(i_neighbgp, cSpace%cGrid%iProcN); 
                ! Compute the index location for the total number of gridpoints that are computed by BubbleProc but are located in iProcID
                GridPoint_Count(BubbleProc + 1) = COUNT(NeighbGPIndex(j*BubbleProc + 1:j*BubbleProc + i)/cSpace%cGrid%iD0LocN == cSpace%cGrid%iProcID)

                ! Compute the index location for the total number of gridpoints for each processor
                if (GridPoint_Count(BubbleProc + 1) > 0) then
                    GridPoint_Disps(BubbleProc + 1) = sum((M + 1)*GridPoint_Count(1:BubbleProc + 1)) - (M + 1)*GridPoint_Count(BubbleProc + 1)
                end if
            end do

            N = sum(GridPoint_Count); 
            ALLOCATE (P_Scattered((M + 1)*N, 1))
            iMemAllocated = iMemAllocated + PRODUCT(SHAPE(P_Scattered))*sizeof(P_scattered(1, 1))
            P_scattered = 0.0D0

            ! This is to automate the proceedure between the selection of single or double precision
            byte_precision = sizeof(P_scattered(1, 1))
            MPI_NEWTYPE = MPI_REAL8; if (byte_precision == spS) MPI_NEWTYPE = MPI_REAL

            do BubbleProc = 0, cSpace%cGrid%iProcN - 1
                ! Count how many points on the running processor are located at the iterated value of processor
                j = count(GridPoint_OnProc(2, :) == BubbleProc)
                ALLOCATE (XYZIndex(j))
                XYZIndex = (/pack(GridPoint_OnProc(1, :) + 1, GridPoint_OnProc(2, :) == BubbleProc)/)
                ! This has the value of the gridpoints relevant to the iterated value of processor
                ! Add +1 because it will be the index for P_scattered_Total
                ! Send then the relevant values from each processor to the iterated value of processor.
                ! In this way, each processor will get the computed pressure field  only on the gridpoints that are located in the respective processor
                call MPI_Gatherv(P_scattered_Total(:, XYZIndex), (M + 1)*j, MPI_NEWTYPE, P_scattered, (M + 1)*GridPoint_Count, GridPoint_Disps, MPI_NEWTYPE, BubbleProc, MPI_COMM_WORLD, iErr)
                DEALLOCATE (XYZIndex)
            end do

            P_scattered_Total = RESHAPE(P_scattered, (/M + 1, N/))
            iMemAllocated = iMemAllocated - PRODUCT(SHAPE(P_Scattered))*sizeof(P_scattered(1, 1))
            DEALLOCATE (P_Scattered)
            DEALLOCATE (GridPoint_OnProc)

            if (N > 0) then
                ALLOCATE (NeighbLocInd(M*N)); ALLOCATE (TimePoints(M))
                TimePoints = (/(i + iDimT - M, i=1, M)/); NeighbLocInd = (/(INT(P_scattered_Total(1, j) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN)*iDimT + TimePoints, j=1, N)/)
                cSpace%cGrid%parD0(NeighbLocInd) = cSpace%cGrid%parD0(NeighbLocInd) + RESHAPE(P_scattered_Total(2:M + 1, :)*(cModelParams%freq0/cMediumParams%c0), (/M*N/))
                DEALLOCATE (NeighbLocInd); DEALLOCATE (TimePoints)
            end if
        end if

        iMemAllocated = iMemAllocated - PRODUCT(SHAPE(P_scattered_Total))*sizeof(P_scattered_Total(1, 1))
        DEALLOCATE (P_scattered_Total)
        iMemAllocated = iMemAllocated - PRODUCT(shape(NeighbGPIndex))*i4bs
        iF (ALLOCATED(NeighbGPIndex)) DEALLOCATE (NeighbGPIndex)
        write (acTemp, "('Proceedure finished.')"); call PrintToLog(acTemp, 3)

    END SUBROUTINE CalculateCloudPressure

    SUBROUTINE LinScatterer(cSpace, InputVal, OutputVal, Magnitude, OrderFFT)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function INTERP1DFREQ returns a signal upsampled or downsampled
        !   based on the input signal and the output length.
        !   If (InputLen > FinalLen) then decimate ( downsample)
        !   If (InputLen < FinalLen) then upsample
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   InputVal         r   dp   a vector of the (time) values of the initial signal
        !   OutputLen        i   i8b    Length of output signal
        !   OutputVal        r   dp   a vector of OutputLen length of the  output signal
        !
        real(dp), intent(in)                ::              InputVal(:)
        real(dp), intent(out)               ::              OutputVal(:)
        integer(i8b), intent(in)            ::              OrderFFT
        real(dp), intent(in)                ::              Magnitude(:)
        type(Space), intent(in)             ::              cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   InputValW                       dpc   Frequency components of input signal
        !   OutputValW                             dpc   Frequency components of output signal
        !   InputLen                                   i8b   Length of Input Signal
        !   EvenOrOdd                       i8b   Variable, 0 if even or 1 if odd
        !   iNumPlanTransform             i8b   plan for FFT
        !   iNumPlanTransform_inv        i8b   plan for IFFT
        !   dTaperMaxFreqWindow     dp    Tapering of Max Frequency to remove artifact
        !   dLeftBand                             dp    Left limit of banded tapering
        !   dRightBand                                   dp    Right limit of banded tapering
        !

        complex(dpc)                        ::              InputValW(size(InputVal) + 1), dMultFactor(size(InputVal) + 1)
        complex(dpc), allocatable           ::              OutputValW(:)
        integer(i8b)                        ::              InputLen, OutputLen, EvenOrOdd, iNumPlanTransform, iNumPlanTransform_inv, iDimW

        !Filtering and windowing parameters
        real(dp)                            ::              dTaperMaxFreqWindow(size(InputVal)/2 + 1), dLeftBand, dRightBand, dDomega, dDt
        integer                             ::              i, iErr, ErrCode

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

        ! Same proceedure from ReorderDistr1to0 or the other way , but it is implemented for 1d arrays
        OutputLen = size(OutputVal)
        InputLen = size(InputVal)
        EvenOrOdd = mod(InputLen, 2_i8b)
        iDimW = InputLen/2 + 1
        dDt = cSpace%dDt !/ (1.0D0*InputLen/cSpace%iDimT)

        InputValW = 0.0
        InputValW(1:InputLen/2) = InputVal(1:InputLen - EvenOrOdd:2) + im*InputVal(2:InputLen - EvenOrOdd:2)
        if (EvenOrOdd == 1) then
            InputValW(InputLen/2 + EvenOrOdd) = InputVal(InputLen)
        end if

        ! Do a FFT of InputLen-point (Initial Variable Length) with complex numbers
        call dfftw_plan_dft_r2c_1d(iNumPlanTransform, InputLen, InputValW, InputValW, FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(iNumPlanTransform, InputValW, InputValW)
        call dfftw_destroy_plan(iNumPlanTransform)

        !Allocate to upsample or decimate . OutputValW has to have more or equal points to the initial signal
        ALLOCATE (OutputValW(max(InputLen + 1, OutputLen + 1)))

        ! Tapering of the max frequency because of artifact created due to higher frequencies
        dLeftBand = 0.0
        dRightBand = 0.1_dp
        dTaperMaxFreqWindow = dTaperingWindow(iDimW, 1.0_dp/(InputLen*dDt), dLeftBand, dRightBand)

        dMultFactor = 1.0D0
        if (OrderFFT /= 0) then
            dDOmega = two_pi/real(InputLen*dDt, dp)
            dMultFactor(1:iDimW) = ((/(i, i=0, iDimW - 1)/)*dDOmega*im)**real(OrderFFT, dp)
        end if

        ! Zero-padding after the values which will be used for interpolation
        ! This is done by initializing OutputValW with 0 ( Really small number due to precision)
        ! Replace the first half part with the values of initial variable.
        OutputValW = 0.0D0
        OutputValW(1:iDimW) = InputValW(1:iDimW)*dMultFactor(1:iDimW)/real(InputLen, dp)*dTaperMaxFreqWindow

        ! Inverse FFT of OutputLen-point
        call dfftw_plan_dft_c2r_1d(iNumPlanTransform_inv, OutputLen, OutputValW, OutputValW, FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(iNumPlanTransform_inv, OutputValW, OutputValW)
        call dfftw_destroy_plan(iNumPlanTransform_inv)

        OutputVal(1:OutputLen:2) = real(OutputValW(1:OutputLen/2 + EvenOrOdd), dp)
        OutputVal(2:OutputLen:2) = dimag(OutputValW(1:OutputLen/2))
        OutputVal = OutputVal*Magnitude; 
        ! OutputVal = OutputValD(1:InputLen) * Magnitude !* dTaperingWindow(InputLen,dDt,2.0_dp, 2.0_dp) ;

        DEALLOCATE (OutputValW)
    END SUBROUTINE LinScatterer

    SUBROUTINE NonLinScatterer(cSpace, Magnitude, OrderFFT, OrderP)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function INTERP1DFREQ returns a signal upsampled or downsampled
        !   based on the input signal and the output length.
        !   If (InputLen > FinalLen) then decimate ( downsample)
        !   If (InputLen < FinalLen) then upsample
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   InputVal         r   dp     a vector of the (time) values of the initial signal
        !   OutputLen        i   i8b    Length of output signal
        !   OutputVal        r   dp     a vector of OutputLen length of the  output signal
        !
        integer(i8b), intent(in)                    ::              OrderFFT, OrderP
        real(dp), intent(in)                        ::              Magnitude(:)
        type(Space), intent(inout)                  ::              cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   InputValW               dpc   Frequency components of input signal
        !   OutputValW              dpc   Frequency components of output signal
        !   InputLen                i8b   Length of Input Signal
        !   EvenOrOdd               i8b   Variable, 0 if even or 1 if odd
        !   iNumPlanTransform       i8b   plan for FFT
        !   iNumPlanTransform_inv   i8b   plan for IFFT
        !   dTaperMaxFreqWindow     dp    Tapering of Max Frequency to remove artifact
        !   dLeftBand               dp    Left limit of banded tapering
        !   dRightBand              dp    Right limit of banded tapering
        !

        complex(dpc), allocatable                   ::              dMultFactor(:)
        integer(i8b)                                ::              InputLen, OutputLen, EvenOrOdd, iNumPlanTransform, iNumPlanTransform_inv, iDimW

        !Filtering and windowing parameters
        real(dp)                                    ::              dTaperMaxFreqWindow(size(cSpace%cGrid%parD0)/2 + 1), dLeftBand, dRightBand, dDomega, dFFTFactor
        real(dp)                                    ::              arBuffer1(size(cSpace%cGrid%parD0) + 1), arBuffer2(size(cSpace%cGrid%parD0) + 1)
        integer                                     ::              i, iErr, iStart, iIndex

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

        ! Same proceedure from ReorderDistr1to0 or the other way , but it is implemented for 1d arrays
        InputLen = size(cSpace%cGrid%parD0)
        EvenOrOdd = mod(InputLen, 2_i8b)
        iDimW = InputLen/2 + 1

        call ReorderDistr0ToDistr1(cSpace%cGrid)
        call TransformT_sml(cSpace%cGrid)

        if (cModelParams.UseFreqTapering .EQV. .true.) then
            !Tapering of the highest frequency part; otherwise the chopoff noise around the
            ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
            ! and is blown up by the double derivative...

            !build tapering window, the highest .1*f0 are tapered in the contrast source
            dLeftBand = 0
            dRightBand = 0.1_dp
            dTaperMaxFreqWindow = dTaperingWindow(iDimW, 1.0_dp/(cSpace%iDimT*cSpace%dDt), dLeftBand, dRightBand)

            !multiply the field with it
            iStart = 0
            do iIndex = 0, cSpace%cGrid%iD1LocN - 1
                cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW) = &
                    cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW)*dTaperMaxFreqWindow

                iStart = iStart + cSpace%cGrid%iD1IS
            end do

        end if
        call TransformTInv(cSpace%cGrid)

        iStart = 0
        do iIndex = 0, cSpace.cGrid.iD1LocN - 1
            arBuffer1 = real(cSpace%cGrid%pacD1(iStart + 1:iStart + cSpace%iDimT), dp)
            arBuffer2 = dimag(cSpace%cGrid%pacD1(iStart + 1:iStart + cSpace%iDimT))

            cSpace%cGrid%pacD1(iStart + 1:iStart + cSpace%iDimT) = arBuffer1**OrderP + im*arBuffer2**OrderP; 
            iStart = iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex

        end do

        ALLOCATE (dMultFactor(iDimW))

        dMultFactor = 1.0D0
        dFFTFactor = 1.0_dp/real(2.0_dp*real(cSpace%cGrid%iD0TL, dp)**(2 + OrderP - 1), dp)

        if (OrderFFT /= 0) then
            dDOmega = two_pi/real(InputLen*cSpace%dDt, dp)
            dMultFactor(1:iDimW) = real(((/(i, i=0, iDimW - 1)/)*dDOmega)**OrderFFT*dFFTFactor, dp)*im**OrderFFT
        end if

        !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
        !Only real multiplication..
        if (cModelParams.UseFreqTapering .EQV. .true.) then
            !Here, the max freq tapering is applied as well
            dMultFactor = dMultFactor*dTaperMaxFreqWindow
        end if

        !Transform to W-domain and multiply with -w**2 and with
        ! FFT normalization factor (1/N for each FORWARD transform,
        ! and an extra square because of the square factor above;
        ! this gives (1/N)**2*(1/(2*N)) )
        call TransformT(cSpace%cGrid)

        iStart = 0
        do iIndex = 0, cSpace.cGrid.iD1LocN - 1
            arBuffer1 = real(cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW), dp)
            arBuffer2 = dimag(cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW))

            cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW) = cSpace%cGrid%pacD1(iStart + 1:iStart + iDimW)*dMultFactor; 
            iStart = iStart + cSpace.cGrid.iD1TL; 
        end do
        call TransformTInv_sml(cSpace%cGrid)
        call ReorderDistr1ToDistr0(cSpace%cGrid)
        cSpace%cGrid%parD0 = cSpace%cGrid%parD0*Magnitude; 
        DEALLOCATE (dMultFactor)

    END SUBROUTINE NonLinScatterer

    SUBROUTINE SaveExtraBubbleSlices(cSpace, dBubbleLocationN, Domain_Range)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 200514
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200514  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine SaveExtraBubbleSlices generates extra slices for output
        !         based on the values inserted in input. The slices are generated based on the
        !        number of slices in each dimension (x,y or z). Then the domain of the bubble
        !        cloud is divided in slices of a number equal to the aforementioned number.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace             i    type(space)     Space from which the data array needs to be stored.
        !   dBubbleLocationN   io   dp              a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   Domain_Range       i    dp              Range of the domain that includes the microbubble cloud
        !
        ! *****************************************************************************
        real(dp), intent(inout)                 ::                dBubbleLocationN(:, :)
        real(dp), intent(in)                    ::                Domain_Range(2, 3)
        type(Space), intent(in)                 ::                cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   iErr            i8b    Error number
        !   iLi             i8b    Loop counter over the blocks to save
        !   iLast           i8b    Index of the last full block
        !   iBlockL         i8b    Size of the blocks
        !   acFilename      char   Total filename including path and suffix
        !   acTemp          char   Temporary char array for output log messages
        !
        ! *****************************************************************************
        character(len=1024)                     ::              filename, actemp
        real(dp)                                ::              DimStart(3), CountPos(3), dLambdaMM
        integer(i8b)                            ::              i, j(3), slicedim, ExtraSliceNum
        character(LEN=1)                        ::              xyzslicedim(cModelParams%Numslices)                !Which dimension: x, y, z or t
        real(dp)                                ::              xyzslicepos(cModelParams%Numslices)                !Position in [mm]
        integer(i8b)                            ::              xyzslicebeam(cModelParams%Numslices)        !Which beam the slice is located; -1 if it may be all
        integer(i8b)                            ::              xyzsliceindex(cModelParams%Numslices)        !Index within the beam
        character(len=1)                        ::              allslicexyzdim(ScattererParams(cSourceParams%iScCloud)%ClusterSlicesN)

        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries and storing of the data array to temporary file
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   FFT in time domain
        !02          Spectral difference by multiplying (iw)^n, ( n = 2 for second derivative)
        !03   IFFT and return to time domain signal
        !04          Multiply with the constant factor A (magnitude)
        !
        ! =============================================================================

        ExtraSliceNum = ScattererParams(cSourceParams%iScCloud)%ClusterSlicesN
        allslicexyzdim = ScattererParams(cSourceParams%iScCloud)%ClusterSliceDim
        dLambdaMM = (cMediumParams.c0*1.0D3)/cModelParams%freq0; ! Normalization factor

        if (ExtraSliceNum == 0) then
            write (acTemp, "('0 Slices were chosen in the input file')")
            call PrintToLog(acTemp, 2)
            return
        end if

        xyzslicedim = (cModelParams%xyzslicedim)
        xyzslicepos = (cModelParams%xyzslicepos)
        xyzslicebeam = (cModelParams%xyzslicebeam)
        xyzsliceindex = (cModelParams%xyzsliceindex)

        ! Deallocate to change dimension
        DEALLOCATE (cModelParams%xyzslicedim)
        DEALLOCATE (cModelParams%xyzslicepos)
        DEALLOCATE (cModelParams%xyzslicebeam)
        DEALLOCATE (cModelParams%xyzsliceindex)

        ! Change dimensions and reallocate
        cModelParams%Numslices = cModelParams%Numslices + ExtraSliceNum

        ALLOCATE (cModelParams%xyzslicedim(cModelParams%Numslices))
        ALLOCATE (cModelParams%xyzslicepos(cModelParams%Numslices))
        ALLOCATE (cModelParams%xyzslicebeam(cModelParams%Numslices))
        ALLOCATE (cModelParams%xyzsliceindex(cModelParams%Numslices))

        cModelParams%xyzslicedim(1:cModelParams%Numslices - ExtraSliceNum) = xyzslicedim
        cModelParams%xyzslicepos(1:cModelParams%Numslices - ExtraSliceNum) = xyzslicepos
        cModelParams%xyzslicebeam(1:cModelParams%Numslices - ExtraSliceNum) = xyzslicebeam
        cModelParams%xyzsliceindex(1:cModelParams%Numslices - ExtraSliceNum) = xyzsliceindex
        ! Count how many slices are for each dimension to divide the domain
        cModelParams%xyzslicedim(cModelParams%Numslices - ExtraSliceNum + 1:cModelParams%Numslices) = allslicexyzdim

        CountPos = (/count(allslicexyzdim == 'x'), &
                     count(allslicexyzdim == 'y'), &
                     count(allslicexyzdim == 'z')/)

        j = 1
        DimStart = (/cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ/)*cSpace%dDx
        ! Separate the domain of the bubble cloud in a number of slices equal to the number of slices given by the user in the input file
        do i = ExtraSliceNum, 1, -1

            slicedim = sum(abs((/'x', 'y', 'z'/) == cModelParams%xyzslicedim(cModelParams%Numslices - i + 1))*(/1, 2, 3/))
            ! If the counter is 1 , save the average position of the bubble cluster in this dimension
            if (ExtraSliceNum == 1) then
                cModelParams%xyzslicepos(cModelParams%Numslices - i + 1) = sum((sum(dBubbleLocationN, 1)/size(dBubbleLocationN, 1) &
                                                                                + DimStart)*(abs((/'x', 'y', 'z'/) == cModelParams%xyzslicedim(cModelParams%Numslices - i + 1))))
            else
                cModelParams%xyzslicepos(cModelParams%Numslices - i + 1) = Domain_Range(1, slicedim)/dLambdaMM &
                                                                           + (Domain_Range(2, slicedim) - Domain_Range(1, slicedim))/(CountPos(slicedim) - 1)*(j(slicedim) - 1)/dLambdaMM
            end if

            j(slicedim) = j(slicedim) + 1

        end do

        if ((cModelParams%Slicesavespecifier == iAI_FIRSTANDLAST) .or. (cModelParams%Slicesavespecifier == iAI_ALL)) then
            ! Initialize slices, because of the addition of the extra slices
            call InitSaveSlices(cSpace)

            do i = ExtraSliceNum, 1, -1
                filename = trim(trim(sOutputDir)//trim(sOutputRoot)//int2str(0_i8b))//'_'//cModelParams%xyzslicedim(cModelParams%Numslices - i + 1)// &
                           int2str(cModelParams%Numslices - i + 1)//int2str(0)
                ! Export slice

                if (cModelParams%xyzslicedim(cModelParams%Numslices - i + 1) == 'x') then
                    call ExportSlice(trim(filename), "p0", cSpace, &
                                     (/0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices - i + 1), 0_i8b, 0_i8b/), &
                                     (/cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ/), &
                                     cModelParams%xyzsliceindex, cSpace%iDimX, .true.); 
                elseif (cModelParams%xyzslicedim(cModelParams%Numslices - i + 1) == 'y') then
                    call ExportSlice(trim(filename), "p0", cSpace, &
                                     (/0_i8b, 0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices - i + 1), 0_i8b/), &
                                     (/cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ/), &
                                     cModelParams%xyzsliceindex, cSpace%iDimY, .true.); 
                elseif (cModelParams%xyzslicedim(cModelParams%Numslices - i + 1) == 'z') then
                    call ExportSlice(trim(filename), "p0", cSpace, &
                                     (/0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices - i + 1)/), &
                                     (/cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b/), &
                                     cModelParams%xyzsliceindex, cSpace%iDimZ, .true.); 
                end if
            end do
        end if

    END SUBROUTINE SaveExtraBubbleSlices

    SUBROUTINE MPI_WRITE_CONTRAST(cSpace, dBubbleContrast, MPI_FILENAME)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 200514
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200514  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine SaveExtraBubbleSlices generates extra slices for output
        !         based on the values inserted in input. The slices are generated based on the
        !        number of slices in each dimension (x,y or z). Then the domain of the bubble
        !        cloud is divided in slices of a number equal to the aforementioned number.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                  io      type(space)     Space from which the data array needs to be stored.
        !   dBubbleLocationN        r       dp              a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   Domain_Range            r       dp              Range of the domain that includes the microbubble cloud
        !
        ! *****************************************************************************
        real(dp), intent(in)                            ::          dBubbleContrast(:)
        type(Space), intent(in)                         ::          cSpace
        character(len=1024), intent(in)                 ::          MPI_FILENAME

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   iErr                    i8b    Error number
        !   iLi                     i8b    Loop counter over the blocks to save
        !   iLast                   i8b    Index of the last full block
        !   iBlockL                 i8b    Size of the blocks
        !   acFilename              char   Total filename including path and suffix
        !   acTemp                  char   Temporary char array for output log messages
        !
        ! *****************************************************************************
        integer                                         ::          ibDimT
        integer                                         ::          FILE_HANDLE, FILE_IN, FILE_OUT, FILE_STATUS(MPI_STATUS_SIZE), MPI_ERROR, MPI_BUF_COUNT
        integer(KIND=MPI_OFFSET_KIND)                   ::          FILE_OFFSET
        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries and storing of the data array to temporary file
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   FFT in time domain
        !02          Spectral difference by multiplying (iw)^n, ( n = 2 for second derivative)
        !03   IFFT and return to time domain signal
        !04          Multiply with the constant factor A (magnitude)
        !
        ! =============================================================================
        ibDimT = size(dBubbleContrast, 1)/ScattererParams(cSourceParams%iScCloud)%N

        FILE_IN = ibDimT*ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN*cSpace%cGrid%iProcID
        FILE_OUT = ibDimT*ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN*(cSpace%cGrid%iProcID + 1)
        IF (cSpace%cGrid%iProcID == cSpace%cGrid%iProcN - 1) FILE_OUT = ibDimT*ScattererParams(cSourceParams%iScCloud)%N
        MPI_BUF_COUNT = FILE_OUT - FILE_IN; 
        FILE_OFFSET = INT(ibDimT*ScattererParams(cSourceParams%iScCloud)%N/cSpace%cGrid%iProcN*cSpace%cGrid%iProcID*dpS, kind=MPI_OFFSET_KIND)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, MPI_FILENAME, MPI_MODE_CREATE + MPI_MODE_RDWR, MPI_INFO_NULL, FILE_HANDLE, MPI_ERROR)
        call MPI_FILE_WRITE_AT_ALL(FILE_HANDLE, FILE_OFFSET, dBubbleContrast(FILE_IN + 1:FILE_OUT), MPI_BUF_COUNT, MPI_DOUBLE_PRECISION, FILE_STATUS, MPI_ERROR)
        call MPI_FILE_CLOSE(FILE_HANDLE, MPI_ERROR)

    END SUBROUTINE MPI_WRITE_CONTRAST

    SUBROUTINE PointSourceCloudOperator(cSpace)

        USE IFPORT
        USE MPI
        ! ===============================================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 211108
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   3.0     200212  Original code (KH)
        !
        ! ***********************************************************************************************
        !
        !   DESCRIPTION
        !
        !   The subroutine PointSourceCloudOperator computes the pressure field generated from a point
        !   source cloud. The position of the point sources is generated here, too.
        !   The data in cSpace should be in Distribution 0.
        !
        ! ***********************************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace              io   type(space)  The space in which the point has to be given
        !

        type(Space), target, intent(inout)      ::          cSpace

        ! ***********************************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   dProtonLocationN            dp    Global coordinate of the location of the point source
        !   dPointSourceContrast        dp    array containing the contrast source term from all the
        !   GlobalPos                   dp    global coordinate of the 8 neighbouring grid points
        !   dLambdaMM                   dp    Normalization factor
        !   DomaiRange                  dp    Output of the point source cloud generator, Domain boundaries
        !   i                           i8b   loop counter
        !   iDimT                       dp    the length of the discrete time vector
        !   arBufferXYZ                 i8b   Temporary buffer containing a space trace
        !
        ! ***********************************************************************************************

        integer                                 ::          iErr
        character(len=1024)                     ::          actemp
        real(dp)                                ::          dLambdaMM, Domain_Range(2, 3), dFFTFactor
        integer(i8b)                            ::          iDimT, i, iPS

        real(dp), allocatable                   ::          dPointSourceContrast(:), dProtonLocationN(:, :)

        !Filtering and windowing parameters
        integer(i8b)                            ::          iDimW, ibDimW, iStart, iIndex
        real(dp)                                ::          dTaperSupportWindow(cSpace%iDimT), dTaperMaxFreqWindow(cSpace%iDimT), dLeftBand, dRightBand
        LOGICAL(lgt)                            ::          exists, exists_loc, result

        ! *****************************************************************************
        !
        !   I/O
        !
        !   log file entries
        !
        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   PrintToLog
        !   CalculateCloudMonoEnergeticPressure
        !   GeneratePSCloud
        !   DeltaPulse
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   Generate the point source cloud using RNDGeneratePSCloud or GenPSCloud
        !02   Generate the time signature for one point source using DeltaPulse
        !03   Calculate the generated pressure field from all the point sources
        !
        ! =============================================================================
        write (acTemp, "('CloudPointSourceOperator')"); call PrintToLog(acTemp, 2)
        call ScattererInit()
        INQUIRE (FILE=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)), EXIST=exists_loc)
        if (.NOT. exists_loc) result = MAKEDIRQQ(trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)))
        !**************************************** INITIALIZE VALUES **********************************************
        iDimT = cSpace%iDimT   !t dimensions

        iDimW = iDimT/2 + 1
        dLambdaMM = (cMediumParams.c0*1.0D3)/cModelParams%freq0; ! Normalization factor λ in [mm]

        !Allocate for memory management
        ALLOCATE (dProtonLocationN(PointSourceCloudParams%N, 3))
        ALLOCATE (dPointSourceContrast(iDimT))

        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(dProtonLocationN)))*dpS

        !Initialize for this variables , close to 0, double precission
        dProtonLocationN = 1.0D-30
        dPointSourceContrast = 1.0D-30
        !**************************************** CREATE POINT SOURCE CLUSTER / LOAD  POINT SOURCE PARAMETERS **********************************************

        write (acTemp, "('Randomly Position the protons in the domain') "); call PrintToLog(acTemp, 2)

        ! Generate the cluster and save result on slices defined in input file
        ! call RNDGeneratePSCloud_Eff(cSpace,PointSourceCloudParams,dProtonLocationN, Domain_Range)
        dProtonLocationN(1, :) = (/0.0, 0.0, 4.11666666666667/)

        dProtonLocationN(:, 1) = dProtonLocationN(:, 1)/dLambdaMM - cSpace%iStartX*cSpace%dDx
        dProtonLocationN(:, 2) = dProtonLocationN(:, 2)/dLambdaMM - cSpace%iStartY*cSpace%dDx
        dProtonLocationN(:, 3) = dProtonLocationN(:, 3)/dLambdaMM - cSpace%iStartZ*cSpace%dDx

        write (acTemp, "('No. of Protons consisting the cluster: ', I<INT(log10(real(ScattererParams(cSourceParams%iScCloud)%N,dp)))+1>, ' in a volume of ', F5.3, ' [mL].')") PointSourceCloudParams%N, PRODUCT(maxval(dProtonLocationN, 1) - minval(dProtonLocationN, 1))*dLambdaMM**3*1.0D-3; call PrintToLog(acTemp, 2); 
        if (cSpace%cGrid%iProcID <= 0) then

            ! Save each Bubble's Position
            open (11, file=trim(trim(sOutputDir)//trim(ScattererParams(cSourceParams%iScCloud)%sBubbleDir)//'PS_Location'//int2str(cModelParams%iIter)//int2str(cSpace%cGrid%iProcID)), status='NEW')
            do iPS = 1, PointSourceCloudParams%N
                write (11, *) dProtonLocationN(iPS, 1), dProtonLocationN(iPS, 2), dProtonLocationN(iPS, 3)
            end do
            close (11)
        end if

        write (acTemp, "('Calculate the time signature of the point source ')"); call PrintToLog(acTemp, 2)
        call DeltaPulse(cSpace, PointSourceCloudParams%PointSourceAmplitude, INT(2*cModelParams%Fnyq*cSourceParams%Tdelay, 8), dPointSourceContrast, 'gaussian')

        ! Tappering
        dPointSourceContrast = dPointSourceContrast*dTaperingWindow(cSpace%iDimT, cSpace%dDt, 2.0_dp, 2.0_dp)

        !****************************************************** CALCULATE PRESSURE FIELD **********************************************
        ! Iniatilize pressure field with 0
        ! Iterate through each bubble and find the pressure field in the whole domain
        ! Add each field to the previous to find the total pressure distribution in the grid
        cSpace%cGrid%parD0 = 1.0D-20
        write (acTemp, "('Calculate pressure field due to the Bubble Cloud')"); call PrintToLog(acTemp, 2)
        call CalculateCloudMonoEnergeticPressure(cSpace, dProtonLocationN, dPointSourceContrast)
        iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(dProtonLocationN)))*dpS

        ! ! ! call ReorderDistr0ToDistr1( cSpace%cGrid)
        ! ! ! call TransformT_sml( cSpace%cGrid)

        ! ! ! if (cModelParams.UseFreqTapering .EQV. .true.) then
        ! ! ! ! Tapering of the highest frequency part; otherwise the chopoff noise around the
        ! ! ! ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
        ! ! ! ! and is blown up by the double derivative...

        ! ! ! ! build tapering window, the highest .1*f0 are tapered in the contrast source
        ! ! ! dLeftBand = 0
        ! ! ! dRightBand = 0.1_dp
        ! ! ! dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)

        ! ! ! dFFTFactor = 1.0_dp/real(cSpace%cGrid%iD0TL,dp)
        ! ! ! ! multiply the field with it
        ! ! ! iStart = 0
        ! ! ! do iIndex=0, cSpace%cGrid%iD1LocN-1
        ! ! ! cSpace%cGrid%pacD1(iStart+1:iStart+iDimW) = &
        ! ! ! cSpace%cGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow        * dFFTFactor

        ! ! ! iStart=iStart+ cSpace%cGrid%iD1IS
        ! ! ! end do

        ! ! ! end if

        ! ! ! L.D.
        ! ! ! ! Now, transform back to T-domain as though it was a grid with wraparound regions
        ! ! ! ! however, the wraparound regions are now used as anti-aliasing regions...
        ! ! ! call TransformTInv_sml( cSpace%cGrid)
        ! ! ! call ReorderDistr1ToDistr0( cSpace%cGrid)

        DEALLOCATE (dProtonLocationN)
        DEALLOCATE (dPointSourceContrast)

        call PrintToLog("End of Pressure Calculation", 2)

    END SUBROUTINE PointSourceCloudOperator

    SUBROUTINE DeltaPulse(cSpace, PointSourceAmplitude, Tindex, dPointSourceContrast, signal_type)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 211108
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function DeltaPulse returns the time signature of a signal generated from
        !   a dirac function centered at Tindex
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !
        !   cSpace                      io      type(space)  The space in which the point has to be given
        !   PointSourceAmplitude        r       dp           The amplitude of one point source. Assuming MonoEnergetic point sources
        !   dPointSourceContrast        r       dp            a vector with the time signatures of the point sources
        !
        real(dp), intent(out)                   ::          dPointSourceContrast(:)
        real(dp), intent(in)                    ::          PointSourceAmplitude
        integer(i8b), intent(in)                ::          Tindex
        type(Space), intent(in)                 ::          cSpace
        character(len=*), intent(in)            ::          signal_type

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   RealTime                    r       dp      The discrete time vector
        !   dK_c_t                      r       dp      cutoff frequency for time variable
        !   i                           i       i8b     Loop variable
        !   iDimT                       i       i8b     The length of the discrete time vector
        !   DistFromPS                  r       dp      Distance from PS
        !

        integer(i8b)                            ::          iDimT, i
        real(dp)                                ::          dK_c_t, RealTime(cSpace%iDimT)
        real(dp)                                ::          DistFromPS
        complex(dpc)                            ::          fdPointSourceContrast(1:cSpace%iDimT/2 + 1)
        complex(dpc)                            ::          fshapeinPS(cSpace%iDimT, 1), fshapeoutPS(cSpace%iDimT, 1)

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
        !   dSinc
        !
        ! =============================================================================

        iDimT = cSpace%iDimT

        RealTime = [(i, i=cSpace%iStartT, iDimT - 1 + cSpace%iStartT)]*cSpace%dDt   ! Array with the real time unnormalized
        dK_c_t = pi/cSpace%dDt

        ! Because the source amplitude is defined at 1 [um] distance, it should be adjusted with the description of the point source
        ! The peak value is A @ the location of the scatterer, so in order to find A, we need to compute A/(4*pi*1E-6) = PointSourceAmplitude
        DistFromPS = PointSourceCloudParams%Dist_Amplitude; 
        if (signal_type == 'delta') then
            dPointSourceContrast = PointSourceAmplitude*dSinc(dK_c_t*(RealTime - Tindex*cSpace%dDt))
        elseif (signal_type == 'bipolar') then
            dPointSourceContrast = PointSourceAmplitude*exp(-((RealTime - iDimT/2*cSpace%dDt + 2*cSourceParams%Tdelay)/ &
                                                              (0.5D0/2.0_dp))**cSourceParams%power) &
                                   *sin(-two_pi*(RealTime - cSourceParams%Tdelay))*(1.0_dp + dsign(1.0_dp, RealTime))/2.0_dp
        elseif (signal_type == 'gaussian') then
            dPointSourceContrast = PointSourceAmplitude*exp(-((RealTime - cSourceParams%Tdelay)/ &
                                                              (cSourceParams%Twidth/2.0_dp))**cSourceParams%power) &
                                   *sin(two_pi*(RealTime - cSourceParams%Tdelay))*(1.0_dp + dsign(1.0_dp, RealTime))/2.0_dp
        end if
        dPointSourceContrast = dPointSourceContrast*(4*pi*DistFromPS)

    END SUBROUTINE DeltaPulse

    SUBROUTINE CalculateCloudMonoEnergeticPressure(cSpace, dPointSourceLocationN, dPointSourceContrast)

        ! =============================================================================
        !
        !   Programmer: Agisilaos Matalliotakis // Date : 211108
        !
        !   Language: Fortran 90
        !
        !   Version Date    Comment
        !   ------- -----   -------
        !   1.0     200608  Original code (KH)
        !
        ! *****************************************************************************
        !
        !   DESCRIPTION
        !
        !   The function CalculateCloudMonoEnergeticPressure calculates the pressure of the total point
        !        source cloud.  This is done by subdividing the total points of interest for the computation
        !   in a number equal to the number of the cpus used.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace                       io   type(space)  The space in which the point has to be given
        !   dPointSourceLocationN    r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
        !   dPointSourceContrast         r   dp   a vector with the time signatures of the microbubbles
        !
        ! *****************************************************************************
        real(dp), intent(inout)                 ::              dPointSourceLocationN(:, :)
        real(dp), intent(inout)                 ::              dPointSourceContrast(:)
        type(Space), intent(inout)              ::              cSpace

        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   P_Scattered                 dpc     Scattered Pressure Matrix for MKL
        !   dSincMult                   dpc     a matrix that contains spatially filtered dirac
        !   dSincMult                   dpc     a matrix that contains spatially filtered dirac from all point sources
        !   dGlobBXYZ                   r       a matrix of the gridpoints location of each cpu domain (same dim as dPointSourceLocationN)
        !   M, ibDimT                   i       integer value for MKL , Time length
        !   N, GridPointsNum            i       integer value for MKL , Space length
        !   K, BLayer                   i       integer value for MKL , Bubbles per layer
        !   iIndex, i, iBubble          i8b     index variables
        !   B_Div                       i8b     divisor of total numbers of variables
        !   dLambdaMM                   dp      value of wavelength , normalization factor
        !   dK_c                        dp      Kutoff frequency
        !
        ! *****************************************************************************

        real(dp), allocatable                   ::              dSincMult(:, :), dSincMultTotal(:), P_scattered(:, :), P_scattered_Total(:, :)
        real(dp), allocatable                   ::              dGlobBXYZ(:, :)
        integer(i8b), allocatable               ::              NeighbLocInd(:), XYZIndex(:), NeighbGPIndex(:)
        integer, allocatable                    ::              GridPoint_OnProc(:, :), TimePoints(:)

        integer                                 ::              M, N, K, iDimT, BLayer, iDimX, iDimY, iDimZ, i, j, BubbleProc, byte_precision
        integer                                 ::              GridPoint_Count(cSpace%cGrid%iProcN), GridPoint_Disps(cSpace%cGrid%iProcN)
        integer                                 ::              iErr, iStat(MPI_STATUS_SIZE), iREQUEST, MPI_NEWTYPE, ErrCode
        integer(i8b)                            ::              iIndex, iLI, x, y, z, i_count, iBeamIndex_min(3), iBeamIndex_max(3), i_neighbgp
        integer(i8b)                            ::              iBubble, B_Div, GridPointsNum
        real(dp)                                ::              dK_c

        character(len=1024)                     ::              actemp

        ! *****************************************************************************
        !
        !   SUBROUTINES/FUNCTIONS CALLED
        !
        !   FFTW
        !
        ! *****************************************************************************
        !
        !   PSEUDOCODE
        !
        !01   First translate the location indexes into normalized coordinates in (x,y,z)
        !01          and calculate for every point of domain
        !02          Set a good choice for the number of Bubbles per Layer ( should check memory usage)
        !03          Iterate for the number of Layers (Total Bubbles / No. Bubbles per Layer)
        !04          Calculate dirac function matrix
        !05          Do the matrix multiplication based on the space points for calculation
        !
        ! =============================================================================

        iDimT = cSpace%iDimT
        iDimX = cSpace%iDimX
        iDimY = cSpace%iDimY
        iDimZ = cSpace%iDimZ

        dK_c = pi/cSpace%dDx

        if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'all' .OR. cModelParams%iIter == cModelParams%Numiterations(1) - 1) then
            ! In this case, every gridpoint's location is computed
            write (acTemp, "('Calculate the contrast source term at each spatial gridpoint')"); call PrintToLog(acTemp, 3)
            GridPointsNum = cSpace%cGrid%iD0LocN
            ALLOCATE (dGlobBXYZ(GridPointsNum, 3))
            ! This is the way to compute them and in this order. This is important because later they will be used
            ! for the computation of the total scattered pressure and they should be placed in this order
            do i = 0, GridPointsNum - 1
                dGlobBXYZ(1 + i, 1) = (cSpace%cGrid%aiD0Loc(i + 1, 1) + iBeamOffsetX(cSpace, cSpace%cGrid%aiD0Loc(i + 1, 3)))*cSpace%dDx
                dGlobBXYZ(1 + i, 2) = (cSpace%cGrid%aiD0Loc(i + 1, 2) + iBeamOffsetY(cSpace, cSpace%cGrid%aiD0Loc(i + 1, 3)))*cSpace%dDx
                dGlobBXYZ(1 + i, 3) = (cSpace%cGrid%aiD0Loc(i + 1, 3))*cSpace%dDx
            end do

            BLayer = 1000; 
        else

            j = 6
            iBeamIndex_min(3) = nint(minval(dPointSourceLocationN(:, 3))/cSpace%dDx)
            iBeamIndex_min(1) = nint(minval(dPointSourceLocationN(:, 1))/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex_min(3))
            iBeamIndex_min(2) = nint(minval(dPointSourceLocationN(:, 2))/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex_min(3))
            iBeamIndex_min = iBeamIndex_min - j

            iBeamIndex_max(3) = nint(maxval(dPointSourceLocationN(:, 3))/cSpace%dDx)
            iBeamIndex_max(1) = nint(maxval(dPointSourceLocationN(:, 1))/cSpace%dDx) - iBeamOffsetX(cSpace, iBeamIndex_max(3))
            iBeamIndex_max(2) = nint(maxval(dPointSourceLocationN(:, 2))/cSpace%dDx) - iBeamOffsetY(cSpace, iBeamIndex_max(3))
            iBeamIndex_max = iBeamIndex_max + j

            i_neighbgp = (iBeamIndex_max(1) - iBeamIndex_min(1) + 1)*(iBeamIndex_max(2) - iBeamIndex_min(2) + 1)*(iBeamIndex_max(3) - iBeamIndex_min(3) + 1)

            iMemAllocated = iMemAllocated + i_neighbgp*i4bs
            ALLOCATE (NeighbGPIndex(i_neighbgp))
            NeighbGPIndex(1:i_neighbgp) = (/(((z*iDimY*iDimX + y*iDimX + x, x=iBeamIndex_min(1), iBeamIndex_max(1)), y=iBeamIndex_min(2), iBeamIndex_max(2)), z=iBeamIndex_min(3), iBeamIndex_max(3))/)

            ! In this case, only the gridpoints that are added in the main module above, are computed
            ! The number of the neighbouring points for every scatterer equals to size(p,1)**3
            write (acTemp, "('Calclulate the contrast source term in the gridpoints inside the boundaries of the cloud')"); call PrintToLog(acTemp, 3)
            GridPointsNum = i_neighbgp/cSpace%cGrid%iProcN + (cSpace%cGrid%iProcID + 1)/cSpace%cGrid%iProcN*mod(i_neighbgp, cSpace%cGrid%iProcN)! Do this for all processors
            ALLOCATE (dGlobBXYZ(GridPointsNum, 3))
            ! Because all cpus see all the gridpoints, we have to use a way to separate those
            ! that are located in a specific processor. Generate a matrix, where first column is the global value of each gridpoint
            ! and the 2nd value is the processor, where it is located.
            ALLOCATE (GridPoint_OnProc(2, GridPointsNum))

            do i = 0, GridPointsNum - 1
                ! Divide the number of gridpoints based on the number of cpus
                j = NeighbGPIndex((i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID + i + 1); 
                GridPoint_OnProc(1, i + 1) = i
                GridPoint_OnProc(2, i + 1) = j/cSpace%cGrid%iD0LocN

                ! This is the way to compute the normalized global values of each gridpoint  in space.
                dGlobBXYZ(1 + i, 1) = mod(j, cSpace%iDimX)*cSpace%dDx; ; 
                j = (j - mod(j, cSpace%iDimX))/cSpace%iDimX; 
                dGlobBXYZ(1 + i, 2) = mod(j, cSpace%iDimY)*cSpace%dDx; 
                j = (j - mod(j, cSpace%iDimY))/cSpace%iDimY; 
                dGlobBXYZ(1 + i, 3) = j*cSpace%dDx; 
            end do

            ! Increase the number of scatterers per layer, because there are less points than final iteration
            BLayer = 5E3; 
        end if
        M = iDimT; N = GridPointsNum; K = BLayer

        ! Because of big size arrays for clusters higher than 10^5 bubbles
        ! due to memory , a big array can not be created ,so I implement the idea
        ! of layering in the matrix, So I do a loop for every K bubbles and calculate
        ! the scattered pressure. I do this until all the bubbles hve been included
        if (PointSourceCloudParams%N <= K) K = PointSourceCloudParams%N
        B_Div = PointSourceCloudParams%N/K; 
        ALLOCATE (dSincMult(K, N))
        ALLOCATE (dSincMultTotal(N))

        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(dSincMult)) + PRODUCT(SHAPE(dSincMultTotal)))*sizeof(dSincMult(1, 1)) + PRODUCT(SHAPE(dGlobBXYZ))*dpS

        write (acTemp, "('Create ', I3,' Layer(s) of ', I<INT(log10(K*1.0D0))+1>, ' Point Source(s) for Dirac Array')") B_Div, K; call PrintToLog(acTemp, 3)

        if (.NOT. ALL(ABS(dPointSourceContrast) < 1E100_dp)) call PrintToLog("Error: P_scattered !!!", 1)
        dSincMultTotal = 0.0D0
        do iBubble = 1, B_Div ! This is the loop for the number of layers of the bubbles
            ! Calculate the spatial signature in the point source definition for every scatterer
            do i = 0, N - 1
                dSincMult(:, 1 + i) = dSinc(dK_c*(dPointSourceLocationN(K*(iBubble - 1) + 1:K*iBubble, 1) - dGlobBXYZ(1 + i, 1))) &
                                      *dSinc(dK_c*(dPointSourceLocationN(K*(iBubble - 1) + 1:K*iBubble, 2) - dGlobBXYZ(1 + i, 2))) &
                                      *dSinc(dK_c*(dPointSourceLocationN(K*(iBubble - 1) + 1:K*iBubble, 3) - dGlobBXYZ(1 + i, 3)))
            end do
            dSincMultTotal = dSincMultTotal + SUM(dSincMult, 1)
        end do

        write (acTemp, "('Spatial Computation of the pressure field has finished.')"); call PrintToLog(acTemp, 3)

        DEALLOCATE (dSincMult)
        iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(dSincMult)))*sizeof(dSincMult(1, 1))
        ALLOCATE (P_scattered_Total(M + 1, N))
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(P_scattered_Total)))*sizeof(P_scattered_Total(1, 1))
        ! P_scattered  = SPREAD(dSincMultTotal,1, M) * (dK_c/pi)**3
        ! do i = 1,M
        ! P_scattered(i,:) = P_scattered(i,:) * dPointSourceContrast(i)
        ! enddo
        CALL DGEMM('N', 'N', M, N, 1, 1.0D0, RESHAPE(dPointSourceContrast, (/M, 1/)), M, RESHAPE(dSincMultTotal*(dK_c/pi)**3, (/1, N/)), 1, 0.0D0, P_scattered_Total(2:M + 1, :), M)
        ! CALL SGEMM('N','N', M, N, K, 1.0, REAL(RESHAPE(dBubbleContrast(1+K*(iBubble-1)*M:K*iBubble*M)* (dK_c/pi)**3,(/M,K/)),sp), M, dSincMult, K, 0.0, P_scattered, M)
        if (.NOT. ALL(ABS(P_scattered_Total) < 1E100_dp)) call PrintToLog("Error: P_scattered !!!", 1)

        iMemAllocated = iMemAllocated - PRODUCT(SHAPE(dGlobBXYZ))*dpS
        DEALLOCATE (dGlobBXYZ)

        write (acTemp, "('Matrix Multiplication has finished.')"); call PrintToLog(acTemp, 3)
        write (acTemp, "('GridPoints Distribution to each processor starts.')"); call PrintToLog(acTemp, 3)
        if (ScattererParams(cSourceParams%iScCloud)%GridPointsPressure == 'all' .OR. cModelParams%iIter == cModelParams%Numiterations(1) - 1) then
            ALLOCATE (P_scattered(iDimT, N))
            P_scattered = 0.0D0; 
            P_scattered(iDimT - M + 1:iDimT, :) = P_scattered_Total(2:M + 1, :)
            cSpace%cGrid%parD0 = cSpace%cGrid%parD0 + RESHAPE(P_scattered*(cModelParams%freq0/cMediumParams%c0), (/iDimT*N/))
            DEALLOCATE (P_scattered)
        else
            ! We need this row in order to know for which value of the gridpoint, the pressure of each row(gridpoint) is computed so we can send it to the right processor
            P_scattered_Total(1, :) = NeighbGPIndex(1 + (i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID:(i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID + N)*1.0D0

            ! Here, the array for the number of elements GridPoint_Count and
            ! the array for the starting index of the location of the elements per processor GridPoint_Disps is generated.
            GridPoint_Count = 0; 
            GridPoint_Disps = 0; 
            do BubbleProc = 0, cSpace%cGrid%iProcN - 1
                j = i_neighbgp/cSpace%cGrid%iProcN
                i = i_neighbgp/cSpace%cGrid%iProcN + (BubbleProc + 1)/cSpace%cGrid%iProcN*mod(i_neighbgp, cSpace%cGrid%iProcN); 
                ! Compute the index location for the total number of gridpoints that are computed by BubbleProc but are located in iProcID
                GridPoint_Count(BubbleProc + 1) = COUNT(NeighbGPIndex(j*BubbleProc + 1:j*BubbleProc + i)/cSpace%cGrid%iD0LocN == cSpace%cGrid%iProcID)

                ! Compute the index location for the total number of gridpoints for each processor
                if (GridPoint_Count(BubbleProc + 1) > 0) then
                    GridPoint_Disps(BubbleProc + 1) = sum((M + 1)*GridPoint_Count(1:BubbleProc + 1)) - (M + 1)*GridPoint_Count(BubbleProc + 1)
                end if
            end do

            N = sum(GridPoint_Count); 
            ALLOCATE (P_Scattered((M + 1)*N, 1))
            iMemAllocated = iMemAllocated + PRODUCT(SHAPE(P_Scattered))*sizeof(P_scattered(1, 1))
            P_scattered = 0.0D0

            ! This is to automate the proceedure between the selection of single or double precision
            byte_precision = sizeof(P_scattered(1, 1))
            MPI_NEWTYPE = MPI_REAL8; if (byte_precision == spS) MPI_NEWTYPE = MPI_REAL

            do BubbleProc = 0, cSpace%cGrid%iProcN - 1
                ! Count how many points on the running processor are located at the iterated value of processor
                j = count(GridPoint_OnProc(2, :) == BubbleProc)
                ALLOCATE (XYZIndex(j))
                XYZIndex = (/pack(GridPoint_OnProc(1, :) + 1, GridPoint_OnProc(2, :) == BubbleProc)/)
                ! This has the value of the gridpoints relevant to the iterated value of processor
                ! Add +1 because it will be the index for P_scattered_Total
                ! Send then the relevant values from each processor to the iterated value of processor.
                ! In this way, each processor will get the computed pressure field  only on the gridpoints that are located in the respective processor
                call MPI_Gatherv(P_scattered_Total(:, XYZIndex), (M + 1)*j, MPI_NEWTYPE, P_scattered, (M + 1)*GridPoint_Count, GridPoint_Disps, MPI_NEWTYPE, BubbleProc, MPI_COMM_WORLD, iErr)
                DEALLOCATE (XYZIndex)
            end do

            P_scattered_Total = RESHAPE(P_scattered, (/M + 1, N/))
            iMemAllocated = iMemAllocated - PRODUCT(SHAPE(P_Scattered))*sizeof(P_scattered(1, 1))
            DEALLOCATE (P_Scattered)
            DEALLOCATE (GridPoint_OnProc)

            if (N > 0) then
                ALLOCATE (NeighbLocInd(M*N)); ALLOCATE (TimePoints(M))
                TimePoints = (/(i + iDimT - M, i=1, M)/); NeighbLocInd = (/(INT(P_scattered_Total(1, j) - cSpace%cGrid%iProcID*cSpace%cGrid%iD0LocN)*iDimT + TimePoints, j=1, N)/)
                cSpace%cGrid%parD0(NeighbLocInd) = cSpace%cGrid%parD0(NeighbLocInd) + RESHAPE(P_scattered_Total(2:M + 1, :)*(cModelParams%freq0/cMediumParams%c0), (/M*N/))
                DEALLOCATE (NeighbLocInd); DEALLOCATE (TimePoints)
            end if
        end if

        iMemAllocated = iMemAllocated - PRODUCT(SHAPE(P_scattered_Total))*sizeof(P_scattered_Total(1, 1))
        DEALLOCATE (P_scattered_Total)
        iMemAllocated = iMemAllocated - PRODUCT(shape(NeighbGPIndex))*i4bs
        iF (ALLOCATED(NeighbGPIndex)) DEALLOCATE (NeighbGPIndex)

        write (acTemp, "('Proceedure finished.')"); call PrintToLog(acTemp, 3)

    END SUBROUTINE CalculateCloudMonoEnergeticPressure

END MODULE ParnacPointSourceCloud
