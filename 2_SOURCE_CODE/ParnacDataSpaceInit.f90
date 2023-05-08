
MODULE ParnacDataSpaceInit

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
!   The module ParnacDataSpaceInit contains subroutines and functions used in
!   the initialization, destruction and basic use of the Space, Grid and
!   RRMSNorm structures.
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

#ifndef RUNESTIMATOR
    USE ParnacTransformInit, ONLY: &
        GridInitTransforms, GridDestroyTransforms
#endif

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
!   SpaceInfo       sub   Echo the descriptive parameters in a Space structure
!                         to the general log file
!   InitSpaceFromPhysicalDims   sub   Initialize a Space structure from the
!                                     supplied physical dimensions.
!   InitSpace       sub   Initialize a Space structure from the supplied grid
!                         dimensions.
!   InitConvolutionSpaces   sub   Initialize two Space structures used in the
!                                 evaluation of the convolution, in the
!                                 subroutine ConvolutionGreen
!   DestructSpace   sub   Destruct a Space structure
!   InitGrid        sub   Initialize a Grid structure
!   GridDistr0Allocate      sub   Allocate a Distribution 0 array in a Grid
!   GridDistr0Deallocate    sub   Deallocate a Distribution 0 array in a Grid
!   GridDistr0CreateEmpty   sub   Create a Distribution 0 array in a Grid and
!                                 set to zero
!   GridDistr0CreateTest_Int   sub   Create a Distribution 0 array in a Grid
!                                    with a certain content for testing.
!   GridDistr0CreateTest_Cos   sub   Create a Distribution 0 array in a Grid
!                                    with a certain content for testing.
!   GridDistr1Allocate      sub   Allocate a Distribution 1 array in a Grid
!   GridDistr1Deallocate    sub   Deallocate a Distribution 1 array in a Grid
!   GridDistr1CreateEmpty   sub   Create a Distribution 1 array in a Grid
!                                 and set to zero
!   GridDistr2Allocate      sub   Allocate a Distribution 2 array in a Grid
!   GridDistr2Deallocate    sub   Deallocate a Distribution 2 array in a Grid
!   GridDistr2CreateEmpty   sub   Create a Distribution 2 array in a Grid
!                                 and set to zero
!   InitGridSpace   sub   Initialize the dimensions of the beam within the
!                         structure as a function of the dimensions in the
!                         Space structure. Internal routine used by InitGrid.
!   InitGridDistrLookupTables   sub   Initialize the internal parameters
!                                     related to the storage method in the
!                                     Grid structure. Internal routine used by
!                                     InitGrid.
!   iBeamOffsetT    fun   Returns the offset in t as a function of the beam
!                         index in z. Relevant for skew beams only.
!   iBeamOffsetX    fun   Returns the offset in x as a function of the beam
!                         index in z. Relevant for skew beams only.
!   iBeamOffsetY    fun   Returns the offset in y as a function of the beam
!                         index in z. Relevant for skew beams only.
!   ObtainNextSmallPrime   fun   Returns the closest larger integer number that
!                                is a product of a number of small primes
!                                2,3,5,7...
!   Isprodsmallprimes   fun   Returns a true or false depending on the input
!                             being a product of small primes
!   InitSaveSlices  sub   Translates the x/y/z save positions of the slices to
!                         be output to file to beam numbers and beam indices.
!   InitRRMSNorm    sub   Initialize the RRMSNorm structure.
!   InitRRMSPreviousCorrection     sub   Initialize the PreviousCorrection data
!                                        array within the RRMSNorm structure.
!   DestroyRRMSPreviousCorrection  sub   Destroy the PreviousCorrection array
!                                        within the RRMSNorm structure.
!   DestroyRRMSNorm   sub   Destroy the RRMSNorm structure.
!
! =============================================================================
!
CONTAINS

    SUBROUTINE SpaceInfo(cSpace)

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
!   The subroutine SpaceInfo echoes the descriptive parameters in a Space
!   structure to the general log file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   i   type(space)  Space structure for which the information is to
!                             be echoed.
!
        type(Space), intent(in) ::                 cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::                acTemp; 
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
!
! =============================================================================

        write (acTemp, "('Debug information on a space')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Space Identifier:       ', I6, '.')") cSpace%iSpaceIdentifier; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Y-dimension Symmetry:   ', L6, '.')") cSpace%bYSymm; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Beam Tangents')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  TanT: ', F10.3, '.')") cSpace%dTanT; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  TanX: ', F10.3, '.')") cSpace%dTanX; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  TanY: ', F10.3, '.')") cSpace%dTanY; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Sampling parameters')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  FNyq:      ', F10.3, '.')") cSpace%dFnyq; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  delta t:   ', F10.3, '.')") cSpace%dDt; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  delta x:   ', F10.3, '.')") cSpace%dDx; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Number of points in each dimension')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  T: ', I8, '.')") cSpace%iDimT; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  X: ', I8, '.')") cSpace%iDimX; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Y: ', I8, '.')") cSpace%iDimY; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Z: ', I8, '.')") cSpace%iDimZ; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Beam Offset in points')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StartT: ', I8, '.')") cSpace%iStartT; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StartX: ', I8, '.')") cSpace%iStartX; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StartY: ', I8, '.')") cSpace%iStartY; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StartZ: ', I8, '.')") cSpace%iStartZ; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Beam index Start')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  X: ', I8, '.')") cSpace%iBeamIndexStartX; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Y: ', I8, '.')") cSpace%iBeamIndexStartY; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Z: ', I8, '.')") cSpace%iBeamIndexStartZ; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        ! Size of the physical problem
        write (acTemp, "('Starting positions in mus/mm')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StT: ', F10.3, '.')") &
            real(cSpace%iStartT)*cSpace%dDt/cModelParams%freq0*1.0E6; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StX: ', F10.3, '.')") &
            real(cSpace%iStartX)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StY: ', F10.3, '.')") &
            real(cSpace%iStartY)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  StZ: ', F10.3, '.')") &
            real(cSpace%iStartZ)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        write (acTemp, "('Beam lengths in mus/mm')"); call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Lt: ', F10.3, '.')") &
            real(cSpace%iDimT - 1)*cSpace%dDt/cModelParams%freq0*1.0E6; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Lx: ', F10.3, '.')") &
            real(cSpace%iDimX - 1)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Ly: ', F10.3, '.')") &
            real(cSpace%iDimY - 1)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        write (acTemp, "('  Lz: ', F10.3, '.')") &
            real(cSpace%iDimZ - 1)*cSpace%dDx*cMediumParams%c0/cModelParams%freq0*1.0E3; call PrintToLog(acTemp, 1); 
        call PrintToLog(' ', 1); 
        call PrintToLog(' ', 1); 
    END SUBROUTINE SpaceInfo

    SUBROUTINE InitSpaceFromPhysicalDims(cSpace, iSpaceIdentifier, bYSymm, dLx, dLy, dLt, dSz, dLz, dthetaX, dthetaY, dthetaT, dFnyq_in)

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
!   The subroutine InitSpaceFromPhysicalDims initializes a Space structure from
!   physical starting points and sizes generally supplied by the configuration
!   file.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iSpaceIdentifier   i   i8b   Identifier for the type of the variable within
!                                the space
!   bYSymm     i   lgt   Whether the variable exhibits symmetry in y
!   dLx        i   dp    Length of the x-dimension in wavelengths rel. to freq0
!   dLy        i   dp    Length of the y-dimension in wavelengths rel. to freq0
!   dLt        i   dp    Length of the t-dimension in periods relative to freq0
!   dLz        i   dp    Length of the z-dimension in wavelengths rel. to freq0
!   dSz        i   dp    Start  of the z-dimension in wavelengths rel. to freq0
!   dThetaX    i   dp    skew beam angle in X in radians
!   dThetaY    i   dp    skew beam angle in Y in radians
!   dThetaT    i   dp    skew beam angle in T in radians
!   dFnyq_in   i   dp    Nyquist frequency relative to freq0
!   cSpace     o   type(space)  Space structure which is to be initialized
!
        integer(i8b), intent(in) ::           iSpaceIdentifier
        logical(lgt), intent(in) ::    bYSymm; 
        real(dp), INTENT(in) ::        dLx, dLy, dLt, dSz, dLz, dthetaX, dthetaY, dthetaT, dFnyq_in
        type(Space), INTENT(out) ::    cSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::                acTemp; 
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
!
! =============================================================================

        cSpace%dFnyq = dFnyq_in; 
        cSpace%iSpaceIdentifier = iSpaceIdentifier; 
        cSpace%bYSymm = bYSymm; 
        write (acTemp, "('InitSpaceFromPhysicalDims')"); call PrintToLog(acTemp, 1); 
        ! Initialize the angles
        !if (dthetaX==half_pi) then
        if (abs(dthetaX - half_pi) < 1e-5) then                   ! Changed by Agis M. 23/03/2020, the comment above was the initial
            cSpace%dTanX = 0.0_dp
        else
            cSpace%dTanX = 1.0_dp/tan(dthetaX)
        end if
        if (abs(dthetaY - half_pi) < 1e-5) then
            cSpace%dTanY = 0.0_dp
        else
            cSpace%dTanY = 1.0_dp/tan(dthetaY)
        end if
        if (abs(dthetaT - half_pi) < 1e-5) then
            cSpace%dTanT = 0.0_dp
        else
            cSpace%dTanT = 1.0_dp/tan(dthetaT)
        end if

        ! We obtain the step sizes. The factors are due to numerical errors introduced in the algorithm
        cSpace%dDt = 1.0_dp/(2.0_dp*cSpace%dFnyq)                                !time  step
        cSpace%dDx = 1.0_dp/(2.0_dp*cSpace%dFnyq)                                !Step in x, y and z

        ! We obtain the dimensions of the problem
        ! The lengths are always more than sufficient with the ceiling(... -1e-5) construction
        cSpace%iDimT = ceiling(dLt/cSpace%dDt - 1e-5); 
        cSpace%iDimX = ceiling(dLx/cSpace%dDx - 1e-5); 
        cSpace%iDimY = ceiling(dLy/cSpace%dDx - 1e-5); 
        cSpace%iDimZ = ceiling(dLz/cSpace%dDx - 1e-5); 
        ! StartT, StartX, StartY and StartZ give the array starting positions in the grid
        ! ensure that the midpoint in X and Y is at ceiling((Dim+1)/2)
        cSpace%iStartZ = nint(dSz/cSpace%dDx); 
        cSpace%iStartT = 0 + nint(cSpace%iStartZ*cSpace%dTanT)
        cSpace%iStartX = -floor(cSpace%iDimX/2.0_dp) + nint(cSpace%iStartZ*cSpace%dTanX); 
        if (bYSymm .EQV. .false.) then
            cSpace%iStartY = -floor(cSpace%iDimY/2.0_dp) + nint(cSpace%iStartZ*cSpace%dTanY); 
        else
            cSpace%iStartY = 0
        end if

        !iBeamIndexStartX/Y/Z are of importance when the space contains a convolution slice
        ! in that case, we could have beam indices <0, which are stored at the end of the array.
        ! For other types of spaces, these are zero
        cSpace%iBeamIndexStartX = 0
        cSpace%iBeamIndexStartY = 0
        cSpace%iBeamIndexStartZ = 0

    END SUBROUTINE InitSpaceFromPhysicalDims

    SUBROUTINE InitSpace(cSpace, iSpaceIdentifier, bYSymm, iDimT, iDimX, iDimY, iDimZ, iStartT, iStartX, iStartY, iStartZ, &
                         iBeamIndexStartX, iBeamIndexStartY, iBeamIndexStartZ, dFnyq_in, dTanX, dTanY, dTanT)

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
!   The subroutine InitSpace initializes a Space structure from starting
!   indices and integer sizes in the global grid. As the space structures
!   always concern discretized spaces, this is a much safer method for Space
!   definition withinthe program - no roundoff problems. Moreover, it has more
!   flexibility with respect to the starting positions in t/x/y/z.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iSpaceIdentifier   i   i8b   Identifier for the type of the variable within
!                                the space
!   bYSymm     i   lgt   Whether the variable exhibits symmetry in y
!   iDimT      i   i8b   Length of the t-dimension in number of points
!   iDimX      i   i8b   Length of the x-dimension in number of points
!   iDimY      i   i8b   Length of the y-dimension in number of points
!   iDimZ      i   i8b   Length of the z-dimension in number of points
!   iStartT    i   i8b   Start of the t-dimension in the global grid
!   iStartX    i   i8b   Start of the x-dimension in the global grid
!   iStartY    i   i8b   Start of the y-dimension in the global grid
!       BEWARE: The iStartX and iStartY values are not always copied
!               into the Space structure. Often, standard values are
!               used instead!!
!   iStartZ    i   i8b   Start of the y-dimension in the global grid
!   iBeamIndexStartX   i   i8b   see ParnacDataDef
!   iBeamIndexStartY   i   i8b   see ParnacDataDef
!   iBeamIndexStartZ   i   i8b   see ParnacDataDef
!   dTanX    i   dp    tangent related to the skew beam angle in X in radians
!   dTanY    i   dp    tangent related to the skew beam angle in Y in radians
!   dTanT    i   dp    tangent related to the skew beam angle in T in radians
!   dFnyq_in   i   dp    Nyquist frequency relative to freq0
!   cSpace     o   type(space)  Space structure which is to be initialized
!
        integer(i8b), intent(in) ::     iSpaceIdentifier
        logical(lgt), intent(in) ::     bYSymm; 
        integer(i8b), INTENT(in) ::     iDimT, iDimX, iDimY, iDimZ
        integer(i8b), INTENT(in) ::     iStartT, iStartX, iStartY, iStartZ
        integer(i8b), INTENT(in) ::     iBeamIndexStartX, iBeamIndexStartY, iBeamIndexStartZ
        real(dp), INTENT(in) ::         dFnyq_in, dTanX, dTanY, dTanT; 
        type(Space), INTENT(inout) ::   cSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::                acTemp; 
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
!
! =============================================================================

        write (acTemp, "('InitSpace')"); call PrintToLog(acTemp, 2); 
        cSpace.dTanX = dTanX; 
        cSpace.dTanY = dTanY; 
        cSpace.dTanT = dTanT; 
        cSpace%iSpaceIdentifier = iSpaceIdentifier; 
        cSpace%bYSymm = bYSymm; 
        ! We obtain the step sizes. The factors are due to numerical errors introduced in the algorithm
        cSpace%dFnyq = dFnyq_in; 
        cSpace%dDt = 1/(2.0*cSpace%dFnyq)    !time  step
        cSpace%dDx = 1/(2.0*cSpace%dFnyq)    !Step in x, y and z

        ! We obtain the dimensions of the problem
        cSpace%iDimT = iDimT
        cSpace%iDimX = iDimX
        cSpace%iDimY = iDimY
        cSpace%iDimZ = iDimZ

!        ! STATIC VALUES
!        cSpace%iDimT        = 0.74224021e-6
!        cSpace%iDimX        = 110e-6
!        cSpace%iDimY        = 110e-6
!        cSpace%iDimZ        = 110e-6

        ! StartT, StartX, StartY and StartZ give the array starting positions in the grid
        !
        ! If iSpaceIdentifier=iSI_XYZSLICE, then use the given iStartX and iStartY
        ! Otherwise ensure that the midpoint in X and Y is at ceiling((Dim+1)/2)
        cSpace%iStartT = iStartT; 
        if (iSpaceIdentifier == iSI_XYZSLICE) then
            cSpace%iStartX = iStartX
            cSpace%iStartY = iStartY
        else
            cSpace%iStartX = -floor(iDimX/2.0_dp) + nint(iStartZ*dTanX); 
            if (bYSymm .EQV. .false.) then
                cSpace%iStartY = -floor(iDimY/2.0_dp) + nint(iStartZ*dTanY); 
            else
                cSpace%iStartY = 0
            end if
        end if
        cSpace%iStartZ = iStartZ; 
        ! iBeamIndexStartX/Y/Z are of importance when the space contains a convolution
        ! slice, i.e. iSpaceIdentifier = iSI_XYZSLICE. In this case, we can have beam
        ! indices <0, which are stored at the end of the array. For other types of
        ! spaces, these are zero.
        cSpace%iBeamIndexStartX = iBeamIndexStartX
        cSpace%iBeamIndexStartY = iBeamIndexStartY
        cSpace%iBeamIndexStartZ = iBeamIndexStartZ

    END SUBROUTINE InitSpace

    SUBROUTINE InitConvolutionSpaces(cConvSrcSpace, cConvDestSpace, cSrcSpace, cDestSpace)

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
!   The subroutine InitConvolutionSpaces initializes the Space structures used
!   in the evaluation of the convolution - in the subroutine ConvolutionGreen.
!   These Space structures are XYZ convolution slices. Each processor has one
!   XYZ, and therefore iDimT = iProcN-1 which has no physical meaning. The
!   other dimensions and iStart and iBeamStart depend heavily on the source
!   and destination spaces.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSrcSpace       i   type(Space)   Source space for the convolution
!   cDestSpace      i   type(Space)   Destination space for the convolution
!   cConvSrcSpace   o   type(Space)   Source XYZ slice used in the convolution
!   cConvDestSpace  o   type(Space)   Destination XYZ slice used in the
!                                     convolution
!
        type(Space), INTENT(in) ::      cSrcSpace, cDestSpace
        type(Space), INTENT(out) ::     cConvSrcSpace, cConvDestSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iDimX     i8b   Size of Src and Dest Slice in the X-dimension
!   iDimY     i8b   Size of Src and Dest Slice in the Y-dimension
!   iDimZ     i8b   Size of Src and Dest Slice in the Z-dimension
!   iSkewXR   i8b   Extra skewness margin in X on the right side of the beam
!   iSkewXL   i8b   Extra skewness margin in X on the left side of the beam
!   iSkewYR   i8b   Extra skewness margin in Y on the right side of the beam
!   iSkewYL   i8b   Extra skewness margin in Y on the left side of the beam
!   iStartX   i8b   Global start index in X
!   iStartY   i8b   Global start index in Y
!   iStartZ   i8b   Global start index in Z
!   iBeamIndexStartXSrc    i8b   iBeamIndexStartX of the Src Slice
!   iBeamIndexStartYSrc    i8b   iBeamIndexStartY of the Src Slice
!   iBeamIndexStartZSrc    i8b   iBeamIndexStartZ of the Src Slice
!   iBeamIndexStartXDest   i8b   iBeamIndexStartX of the Dest Slice
!   iBeamIndexStartYDest   i8b   iBeamIndexStartY of the Dest Slice
!   iBeamIndexStartZDest   i8b   iBeamIndexStartZ of the Dest Slice
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b) ::                                iDimX, iDimY, iDimZ
        integer(i8b) ::                                iSkewXR, iSkewXL, iSkewYR, iSkewYL, iStartX, iStartY, iStartZ
        integer(i8b) ::                                iBeamIndexStartXSrc, iBeamIndexStartYSrc, iBeamIndexStartZSrc
        integer(i8b) ::                                iBeamIndexStartXDest, iBeamIndexStartYDest, iBeamIndexStartZDest
        character(LEN=2048) ::                    acTemp; 
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
!   ObtainNextSmallPrime
!   InitSpace
!
! =============================================================================

        write (acTemp, "('InitConvolutionSpaces')"); call PrintToLog(acTemp, 3); 
        ! Determine the dimensions, starts and beam starts in the T/X/Y/Z axes
        ! For the T axis, this is simple, since iDimT=1 for an X/Y/Z slice on each processor
        ! For the X, Y and Z axes, this is a more complex issue, which needs a good insight in
        ! the ranges/wraparound regions needed for convolving skewed support domains.
        !
        ! iBeamIndexStartX/Y/Z are of importance when the space contains a Green"s function
        ! in that case, we can have beam indices <0, which are stored at the end of the array
        ! for other types of spaces, these are zero.

        iDimZ = ObtainNextSmallPrime(cSrcSpace%iDimZ + cDestSpace%iDimZ - 1)
        iBeamIndexStartZSrc = -(cDestSpace%iDimZ - 1)
        iBeamIndexStartZDest = -(cSrcSpace%iDimZ - 1)

        if (abs(cSrcSpace%dTanX) < 1.0D-5 .and. abs(cDestSpace%dTanX) < 1.0D-5) then
            iDimX = ObtainNextSmallPrime(cSrcSpace%iDimX + cDestSpace%iDimX - 1)
            iBeamIndexStartXSrc = -(cDestSpace%iDimX - 1)
            iBeamIndexStartXDest = -(cSrcSpace%iDimX - 1)
        else
            iSkewXL = ceiling(max(0.0_dp, cSrcSpace%iDimZ*(cSrcSpace%dTanX - cDestSpace%dTanX)) - 1.0D-5)
            iSkewXR = ceiling(max(0.0_dp, cSrcSpace%iDimZ*(cDestSpace%dTanX - cSrcSpace%dTanX) - cDestSpace%iDimX) - 1.0D-5)
            iDimX = ObtainNextSmallPrime(cSrcSpace%iDimX + cDestSpace%iDimX - 1 + iSkewXL + iSkewXR + 2)
            iBeamIndexStartXSrc = -(cDestSpace%iDimX - 1) - iSkewXL - 1
            iBeamIndexStartXDest = -(cSrcSpace%iDimX - 1) - iSkewXL - 1
        end if

        if (abs(cSrcSpace%dTanY) < 1.0D-5 .and. abs(cDestSpace%dTanY) < 1.0D-5) then
            iDimY = ObtainNextSmallPrime(cSrcSpace%iDimY + cDestSpace%iDimY - 1)
            iBeamIndexStartYSrc = -(cDestSpace%iDimY - 1)
            iBeamIndexStartYDest = -(cSrcSpace%iDimY - 1)
        else
            iSkewYL = ceiling(max(0.0_dp, cSrcSpace%iDimZ*(cSrcSpace%dTanY - cDestSpace%dTanY)) - 1.0D-5)
            iSkewYR = ceiling(max(0.0_dp, cSrcSpace%iDimZ*(cDestSpace%dTanY - cSrcSpace%dTanY) - cDestSpace%iDimY) - 1.0D-5)
            iDimY = ObtainNextSmallPrime(cSrcSpace%iDimY + cDestSpace%iDimY - 1 + iSkewYL + iSkewYR + 2)
            iBeamIndexStartYSrc = -(cDestSpace%iDimY - 1) - iSkewYL - 1
            iBeamIndexStartYDest = -(cSrcSpace%iDimY - 1) - iSkewYL - 1
        end if

        if (cSrcSpace%iSpaceIdentifier == iSI_PRIMARYSOURCE .or. cSrcSpace%iSpaceIdentifier == iSI_FIELDSLICE) then
            ! In this case, iDimZ=1 is sufficient because of the special efficient
            ! implementation for these types in ComputeConvolution()

            call InitSpace(cConvSrcSpace, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
                           cSrcSpace%cGrid%iProcN - 1, iDimX, iDimY, 1_i8b, &
                           0_i8b, cSrcSpace%iStartX, cSrcSpace%iStartY, cSrcSpace%iStartZ, &
                           iBeamIndexStartXSrc, iBeamIndexStartYSrc, 0_i8b, &
                           cSrcSpace%dFnyq, cSrcSpace%dTanX, cSrcSpace%dTanY, cSrcSpace%dTanT)
            ! iDimT has no physical meaning, but we need one slice on each processor

        else

            call InitSpace(cConvSrcSpace, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
                           cSrcSpace%cGrid%iProcN - 1, iDimX, iDimY, iDimZ, &
                           0_i8b, cSrcSpace%iStartX, cSrcSpace%iStartY, cSrcSpace%iStartZ, &
                           iBeamIndexStartXSrc, iBeamIndexStartYSrc, iBeamIndexStartZSrc, &
                           cSrcSpace%dFnyq, cSrcSpace%dTanX, cSrcSpace%dTanY, cSrcSpace%dTanT)
            ! iDimT has no physical meaning, but we need one slice on each processor

        end if

        call InitSpace(cConvDestSpace, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
                       cDestSpace%cGrid%iProcN - 1, iDimX, iDimY, iDimZ, &
                       0_i8b, cDestSpace%iStartX, cDestSpace%iStartY, cDestSpace%iStartZ, &
                       iBeamIndexStartXDest, iBeamIndexStartYDest, iBeamIndexStartZDest, &
                       cDestSpace%dFnyq, cDestSpace%dTanX, cDestSpace%dTanY, cDestSpace%dTanT)
        ! iDimT has no physical meaning, but we need one slice on each processor

    END SUBROUTINE InitConvolutionSpaces

    SUBROUTINE DestructSpace(cSpace)

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
!   The subroutine DestructSpace Destroys the Space structure and deallocates
!   all arrays.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space to be destructed
!
        type(Space), INTENT(inout) ::                cSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::    acTemp; 
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
!   GridDestroyTransforms
!   GridDistr0DeAllocate
!   GridDistr1DeAllocate
!   GridDistr2DeAllocate
!
! =============================================================================

        write (acTemp, "('DestructSpace')"); call PrintToLog(acTemp, 2); 
#ifndef RUNESTIMATOR
        ! Destruct the derived parameters
        call GridDestroyTransforms(cSpace); 
#endif
        ! Deallocate the info arrays
        if (allocated(cSpace%cGrid%aiD0Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD0Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD0Loc); 
        end if
        if (allocated(cSpace%cGrid%aiD1Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD1Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD1Loc); 
        end if
        if (allocated(cSpace%cGrid%aiD2Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD2Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD2Loc); 
        end if

        ! Empty the data stored in this type
        call GridDistr0DeAllocate(cSpace%cGrid); 
        call GridDistr1DeAllocate(cSpace%cGrid); 
        call GridDistr2DeAllocate(cSpace%cGrid); 
    END SUBROUTINE DestructSpace

    SUBROUTINE InitGrid(cSpace, iProcN, iProcID, alResize)

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
!   The subroutine InitGrid initializes the Grid structure within a Space
!   structure. The Space structure has to be initialized first.
!
!   BEWARE:
!       Through the call to InitGridSpace, this function may actually change
!       the dimensions of the Space structure.. This is done to make sure that
!       the data can be distributed nicely over the different processors, and
!       that it has good dimensions for fast FFT's.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace     io   type(Space)   Space for which the Grid is to be initialized
!   iProcN     i    i8b           Total number of processors
!   iProcID    i    i8b           ID number of the current processor
!   alResize   i    lgt           Array of four logical values indicating
!                                 whether each of the dimensions T/X/Y/Z
!                                 may be resized (true) or not (false) in
!                                 InitGrid.
!
        type(Space), intent(inout) ::   cSpace; 
        integer(i8b), intent(in) ::    iProcN, iProcID
        logical, intent(in) ::          alResize(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::    acTemp; 
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
!   InitGridSpace
!   InitGridDistrLookupTables
!   GridInitTransforms
!
! =============================================================================

        write (acTemp, "('InitGrid')"); call PrintToLog(acTemp, 2); 
        ! Initialize the standard variables
        cSpace%cGrid%iProcN = iProcN; 
        cSpace%cGrid%iProcID = iProcID; 
        ! We initialize the space associated to this grid
        call InitGridSpace(cSpace, alResize); 
        ! Now we have the associated space, we can fill in the details for each distribution type
        cSpace%cGrid%iDistr = -1; !This grid is empty, so there is no distribution yet
        call InitGridDistrLookupTables(cSpace); 
!#ifndef RUNESTIMATOR
        ! The last step is to initialize the transforms
        call GridInitTransforms(cSpace); 
        ! These checks shouldn't be necessary, but it seems so due to some bug in Huygens...
        if (associated(cSpace%cGrid%parD0)) then
            call PrintToLog("******** WARNING: Space parD0 is already associated. Use Nullify.", 2)
            nullify (cSpace%cGrid%parD0)
        end if
        if (associated(cSpace%cGrid%pacD1)) then
            call PrintToLog("******** WARNING: Space pacD1 is already associated. Use Nullify.", 2)
            nullify (cSpace%cGrid%pacD1)
        end if
        if (associated(cSpace%cGrid%pacD2)) then
            call PrintToLog("******** WARNING: Space pacD2 is already associated. Use Nullify.", 2)
            nullify (cSpace%cGrid%pacD2)
        end if
!#endif

    END SUBROUTINE InitGrid

    SUBROUTINE GridDistr0Allocate(cGrid)

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
!   The subroutine GridDistr0Allocate allocates the memory for a Distribution 0
!   array in the Grid - if it was not already allocated.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 0 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   err       i8b    error number
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b) :: err
        character(LEN=2048) ::    acTemp; 
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
!   alloctest
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('GridDistr0Allocate')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        if (.NOT. associated(cGrid%parD0)) then

            allocate (cGrid%parD0(cGrid%iD0LocSize), STAT=err); 
            call alloctest(err, 'Error in allocation of a Distr 0 array')
            iMemAllocated = iMemAllocated + cGrid%iD0LocSize*dpS

            if (cGrid%iDistr == -1) then
                cGrid%iDistr = 0; 
            end if

        else
            write (acTemp, "('******ERROR: Was already associated')"); call PrintToLog(acTemp, 4); 
        end if

    END SUBROUTINE GridDistr0Allocate

    SUBROUTINE GridDistr0DeAllocate(cGrid)

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
!   The subroutine GridDistr0DeAllocate deallocates the memory for a
!   Distribution 0 array in the Grid - if it was allocated. Otherwise, it does
!   nothing.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 0 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::    acTemp; 
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
!
! =============================================================================

        write (acTemp, "('GridDistr0Deallocate')"); call PrintToLog(acTemp, 3); 
        if (associated(cGrid%parD0)) then

            iMemAllocated = iMemAllocated - cGrid%iD0LocSize*dpS
            deallocate (cGrid%parD0)
            nullify (cGrid%parD0); 
            if (cGrid%iDistr == 0) then
                cGrid%iDistr = -1; 
            end if

        end if

        write (acTemp, "('Deallocated')"); call PrintToLog(acTemp, 4); 
    END SUBROUTINE GridDistr0DeAllocate

    SUBROUTINE GridDistr0CreateEmpty(cGrid)

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
!   The subroutine GridDistr0CreateEmpty allocates the memory for a
!   Distribution 0 array in the Grid - if it was not already allocated - and
!   puts it to zero.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 0 array needs
!                              to be created.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b    loop counter over all array positions
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b)::                                iLi
        character(LEN=2048) ::                acTemp; 
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
!   GridDistr0Allocate
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('GridDistr0CreateEmpty')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        call GridDistr0Allocate(cGrid)
        cGrid%parD0 = 0.0D0

    END SUBROUTINE GridDistr0CreateEmpty

    SUBROUTINE GridDistr0CreateTest_Int(cGrid)

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
!   The subroutine GridDistr0CreateTest_Int allocates the memory for a
!   Distribution 0 array in the Grid - if it was not already allocated - and
!   fills it with a test input.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 0 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b    loop counter over all xyz positions
!   iLt       i8b    loop counter over all time instants
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b)::                                iLi, iLt
        character(LEN=2048) ::                acTemp; 
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
!
! =============================================================================

        write (acTemp, "('GridDistr0CreateTest_Int')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        call GridDistr0Allocate(cGrid)

        do iLi = 0, cGrid%iD0LocN - 1
            do iLt = 0, cGrid%iD0TL - 1
                cGrid%parD0(1 + iLi*cGrid%iD0IS + iLt*cGrid%iD0TS) = iLt*1000000 + sum(cGrid%aiD0Loc(1 + iLi, :)*(/10000, 100, 1/)); 
            end do
        end do

    END SUBROUTINE GridDistr0CreateTest_Int

    SUBROUTINE GridDistr0CreateTest_Cos(cSpace)

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
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine GridDistr0CreateTest_Cos allocates the memory for a
!   Distribution 0 array in the Grid - if it was not already allocated - and
!   fills it with a test input.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace    io   type(Space)   Space for which the Distribution 0 array needs
!                                to be allocated.
!
        type(Space), INTENT(inout) ::                cSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   dAlpha    dp     angular frequency of test signal
!   iX        i8b    global index in x
!   iY        i8b    global index in y
!   iZ        i8b    global index in z
!   iT        i8b    global index in t
!   iIndex    i8b    array index
!   iLi       i8b    loop counter over all xyz positions
!   iLt       i8b    loop counter over all time instants
!   acTemp    char   temporary char array for output log messages
!
        real(dp) ::                 dAlpha = 0.25; 
        integer(i8b) ::             iX, iY, iZ, iT, iIndex; 
        integer(i8b) ::             iLi, iLt
        character(LEN=2048) ::    acTemp; 
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
!   iBeamOffsetX
!   iBeamOffsetY
!   iBeamOffsetT
!
! =============================================================================

        ! Local variables

        write (acTemp, "('GridDistr0CreateEmptyTest_Cos')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        call GridDistr0Allocate(cSpace.cGrid)

        do iLi = 0, cSpace.cGrid.iD0LocN - 1
            do iLt = 0, cSpace.cGrid.iD0TL - 1
                iZ = cSpace.cGrid.aiD0Loc(1 + iLi, 3); 
                iX = cSpace.cGrid.aiD0Loc(1 + iLi, 1) + cSpace.iStartX + iBeamOffsetX(cSpace, iZ); 
                iY = cSpace.cGrid.aiD0Loc(1 + iLi, 2) + cSpace.iStartY + iBeamOffsetY(cSpace, iZ); 
                iT = iLt + cSpace.iStartT + iBeamOffsetT(cSpace, iZ); 
                iZ = cSpace.cGrid.aiD0Loc(1 + iLi, 3) + cSpace.iStartZ; 
                iIndex = 1 + iLi*cSpace.cGrid.iD0IS + iLt; 
                cSpace.cGrid.parD0(iIndex) = cos(dAlpha*(iZ + iX)); 
            end do
        end do
    END SUBROUTINE GridDistr0CreateTest_Cos

    SUBROUTINE GridDistr1Allocate(cGrid)

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
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine GridDistr1Allocate allocates the memory for a Distribution 1
!   array in the Grid - if it was not already allocated.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 1 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   err       i8b    error number
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b) :: err
        character(LEN=2048) ::    acTemp; 
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
!   alloctest
!
! =============================================================================

        write (acTemp, "('GridDistr1Allocate')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        if (.NOT. associated(cGrid%pacD1)) then

            allocate (cGrid%pacD1(cGrid%iD1LocSize), STAT=err); 
            call alloctest(err, 'Error in allocation of a Distr 1 array')
            iMemAllocated = iMemAllocated + cGrid%iD1LocSize*dpcS

            if (cGrid%iDistr == -1) then
                cGrid%iDistr = 1; 
            end if

        else
            write (acTemp, "('******ERROR: Was already associated')"); call PrintToLog(acTemp, 4); 
        end if

    END SUBROUTINE GridDistr1Allocate

    SUBROUTINE GridDistr1DeAllocate(cGrid)

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
!   The subroutine GridDistr1DeAllocate deallocates the memory for a
!   Distribution 1 array in the Grid - if it was allocated. Otherwise, it does
!   nothing.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   Grid    io   type(cGrid)   Grid for which the Distribution 1 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::    acTemp; 
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
!
! =============================================================================

        write (acTemp, "('GridDistr1DeAllocate')"); call PrintToLog(acTemp, 3); 
        if (associated(cGrid%pacD1)) then

            iMemAllocated = iMemAllocated - cGrid%iD1LocSize*dpcS
            deallocate (cGrid%pacD1)
            nullify (cGrid%pacD1); 
            if (cGrid%iDistr == 1) then
                cGrid%iDistr = -1; 
            end if

        end if

        write (acTemp, "('Deallocated')"); call PrintToLog(acTemp, 4); 
    END SUBROUTINE GridDistr1DeAllocate

    SUBROUTINE GridDistr1CreateEmpty(cGrid)

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
!   The subroutine GridDistr1CreateEmpty allocates the memory for a
!   Distribution 1 array in the Grid - if it was not already allocated - and
!   puts it to zero.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 1 array needs
!                              to be created.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b    loop counter over all array positions
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b)::                                iLi
        character(LEN=2048) ::                acTemp; 
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
!   GridDistr1Allocate
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('GridDistr1CreateEmpty')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        call GridDistr1Allocate(cGrid); 
        cGrid%pacD1 = 0; 
    END SUBROUTINE GridDistr1CreateEmpty

    SUBROUTINE GridDistr2Allocate(cGrid)

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
!   The subroutine GridDistr2Allocate allocates the memory for a Distribution 2
!   array in the Grid - if it was not already allocated.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   Grid    io   type(cGrid)   Grid for which the Distribution 2 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   err       i8b    error number
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b) :: err
        character(LEN=2048) ::    acTemp; 
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
!   alloctest
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('GridDistr2Allocate')"); call PrintToLog(acTemp, 3); 
        ! Allocate the array (BIG array)
        if (.NOT. associated(cGrid%pacD2)) then

            allocate (cGrid%pacD2(cGrid%iD2LocSize), STAT=err)
            call alloctest(err, 'Error in allocation of a Distr 2 array')

            iMemAllocated = iMemAllocated + cGrid%iD2LocSize*dpcS

            if (cGrid%iDistr == -1) then
                cGrid%iDistr = 2; 
            end if

        else
            write (acTemp, "('******ERROR: Was already associated')"); call PrintToLog(acTemp, 4); 
        end if

    END SUBROUTINE GridDistr2Allocate

    SUBROUTINE GridDistr2DeAllocate(cGrid)

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
!   The subroutine GridDistr2DeAllocate deallocates the memory for a
!   Distribution 2 array in the Grid - if it was allocated. Otherwise, it does
!   nothing.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   Grid    io   type(cGrid)   Grid for which the Distribution 2 array needs
!                              to be allocated.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   acTemp    char   temporary char array for output log messages
!
        character(LEN=2048) ::    acTemp; 
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
!
! =============================================================================

        write (acTemp, "('GridDistr2Deallocate')"); call PrintToLog(acTemp, 3); 
        if (associated(cGrid%pacD2)) then

            deallocate (cGrid%pacD2); 
            iMemAllocated = iMemAllocated - cGrid%iD2LocSize*dpcS
            nullify (cGrid%pacD2); 
            if (cGrid%iDistr == 2) then
                cGrid%iDistr = -1; 
            end if

        end if

        write (acTemp, "('Deallocated')"); call PrintToLog(acTemp, 4); 
    END SUBROUTINE GridDistr2DeAllocate

    SUBROUTINE GridDistr2CreateEmpty(cGrid)

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
!   The subroutine GridDistr2CreateEmpty allocates the memory for a
!   Distribution 2 array in the Grid - if it was not already allocated - and
!   puts it to zero.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid    io   type(Grid)   Grid for which the Distribution 2 array needs
!                              to be created.
!
        type(Grid), INTENT(inout) ::                cGrid; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b    loop counter over all array positions
!   acTemp    char   temporary char array for output log messages
!
        integer(i8b)::              iLi
        character(LEN=2048) ::    acTemp; 
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
!   GridDistr0Allocate
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('GridDistr2CreateEmpty')"); call PrintToLog(acTemp, 3); 
        ! Empty the data currently stored in this type
        call GridDistr2Allocate(cGrid); 
        cGrid%pacD2 = 0.0D0; 
    END SUBROUTINE GridDistr2CreateEmpty

    SUBROUTINE InitGridSpace(cSpace, alResize)

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
!   The subroutine InitGridSpace initializes the dimensions of the beam in the
!   grid. The dimensions given by the space are taken as basis, and if allowed
!   and necessary, dimensions are increased to better fit the
!   following requirements:
!   1) Mod(cSpace%iDimX * cSpace%iDimY * cSpace%iDimZ, iProcN) = 0
!   2) Mod(cSpace%iDimT+1, iProcN ) = 0.
!   3) The separate dimensions iDimX, iDimY, iDimZ and iDimT should be
!       products of small primes.
!   The requirements 1) and 2) are necessary in order to distribute the data
!   nicely over the processors in all three distributions. The iDimT+1 in
!   requirement 2)is necessary because of the length of the frequency axes in
!   Distr.1 and 2 being iDimT+1. Note that it is an extremely bad idea to take
!   iProcN = 31 (i.e. a large prime). The requirement 3) ensures that the
!   dimensions are nice for efficient FFT's. Whether or not it is allowed to
!   change a dimension is determined from alResize.
!
!   BEWARE:
!       This function may actually change the dimensions of the Space
!       structure..
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace     io   type(Space)   Space on which InitGridSpace operates.
!   alResize   i    lgt           Array of four logical values indicating
!                                 whether each of the dimensions T/X/Y/Z
!                                 may be resized (true) or not (false) in
!                                 InitGrid.
!
        type(Space), intent(inout) ::                cSpace; 
        logical, intent(in) ::                       alResize(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iDimT    i8b    size of T-dimension
!   iDimX    i8b    size of X-dimension
!   iDimY    i8b    size of Y-dimension
!   iDimZ    i8b    size of Z-dimension
!   i        i8b    temporary variable for finding the optimal values for X/Y/Z
!   j        i8b    temporary variable for finding the optimal values for X/Y/Z
!   k        i8b    temporary variable for finding the optimal values for X/Y/Z
!   iPa      i8b    product of small primes for X
!   iPb      i8b    product of small primes for Y
!   iPc      i8b    product of small primes for Z
!   iK       i8b    multiplication factor for an appropriate T-dimension
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::                iDimT, iDimX, iDimY, iDimZ
        integer(i8b) ::                i, j, k, iPa, iPb, iPc, iK; 
        character(LEN=2048) ::                acTemp; 
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
!   ObtainNextSmallPrime
!   isProdSmallPrimes
!   InitSpace
!   PrintToLog
!
! =============================================================================

        write (acTemp, "('InitGridSpace')"); call PrintToLog(acTemp, 3); 
        ! Obtain the sizes of each dimension
        iDimT = cSpace%iDimT
        iDimX = cSpace%iDimX
        iDimY = cSpace%iDimY
        iDimZ = cSpace%iDimZ

        ! May we resize the T-dimension?
        if (alResize(1)) then
            ! First we will handle the easiest: iDimT
            !        we need iDimT+1 <= \alpha * p = iDimTNEW, with \alpha minimal and iDimTNEW a product of small primes, so:
            iK = ceiling((iDimT + 1)/real(cSpace%cGrid%iProcN, dp)); 
            iDimT = iK*cSpace%cGrid%iProcN - 1; 
!KH How about small primes for T?? Now we simply assume
!KH that both iK and iProcN are a product of small primes...
!KH e.g. with iProcN=4 and iDimT=171 -> iK = 43, which is NOT a product of small primes!
        else

            ! If we may not resize it, then it should be good already.
            ! If not, issue a warning but continue running
            iDimT = cSpace%iDimT; 
            if ((mod(iDimT + 1, cSpace%cGrid%iProcN) /= 0) .and. (cSpace%iSpaceIdentifier /= iSI_INHOMCONTRAST)) then
                write (acTemp, "('Warning while initializing the grid: time dimension was not allowed to be changed.')")
                call PrintToLog(acTemp, -1)
                write (acTemp, "('However, the current choice for iDimT is not correct.')")
                call PrintToLog(acTemp, -1)
                write (acTemp, "('iDimT = ', I5, ' the number of processors is: ', I4, '.')") iDimT, cSpace%cGrid%iProcN
                call PrintToLog(acTemp, -1)
            end if
        end if

        ! The demands on x,y,z are significantly harder to satisfy. We decided to implement it like the case for T
        ! we will first try to find xyz = abc p, with p = p_a p_b p_c; x = a p_a; y = b p_b; z = c p_c.
        ! It is clear that for optimal choices for p_a, p_b and p_c this expression will give an optimal
        ! solution.
        ! We will settle for a heuristic and take p = Prod(p_i), with i \in [0..N], with p_i <= p_{i+1}, primes,
        ! we take:         p_a = Prod(p_i), with i \in [0..N] and i mod 3 == 0
        !                  p_b = Prod(p_i), with i \in [0..N] and i mod 3 == 1
        !                  p_b = Prod(p_i), with i \in [0..N] and i mod 3 == 2                (This is for the case that all three dimensions may be resized)
        ! However, if the conditions are already met for these values of x,y,z, we will just use them
        !KH Alternative approach: make arrays of the first N numbers of small
        !KH primes larger than or equal to the required dimension sizes. Multiply each of
        !KH them with each other to get N**3 numbers and take the smallest one
        !KH that is also a multiple of the number of processors.
        !KH Much easier, much more robust, and negligibly slower.
        if (alResize(2) .or. alResize(3) .or. alResize(4)) then
            if (.NOT. (isProdSmallPrimes(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ) .AND. &
                       (mod(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ, cSpace%cGrid%iProcN) == 0))) then

!KH write (acTemp, '("iDimX, iDimY, iDimZ before:    ",3I4)') cSpace%iDimX,cSpace%iDimY,cSpace%iDimZ
!KH call mexPrintf(acTemp//achar(13))

                iPa = 1; 
                iPb = 1; 
                iPc = 1; 
                i = cSpace%cGrid%iProcN; 
                j = 2; ! The first prime number
                k = 0; ! a counter
                do while (i /= 1)

                    ! If not j | i, we should try the next number
                    do while (mod(i, j) /= 0)
                        j = j + 1; 
                    end do

                    ! Select the first dimension we are allowed to change:
                    do while (.NOT. alResize(mod(k, 3_i8b) + 2))
                        k = k + 1; 
                    end do

                    select case (mod(k, 3_i8b))
                    case (0)
                        iPa = iPa*j; 
                    case (1)
                        iPb = iPb*j; 
                    case (2)
                        iPc = iPc*j; 
                    end select

!KH write (acTemp, '("i, j, k:    ",3I4)') i,j,k
!KH call mexPrintf(acTemp//achar(13))

                    i = i/j; 
                    k = k + 1; 
                end do        ! Note that this method will only yield primes
                ! BEWARE: this works only if ALL alResize(2:4) are .true.
                ! glitch...
                iDimX = ObtainNextSmallPrime(int(ceiling(iDimX/real(iPa)), i8b))*iPa; 
                iDimY = ObtainNextSmallPrime(int(ceiling(iDimY/real(iPb)), i8b))*iPb; 
                iDimZ = ObtainNextSmallPrime(int(ceiling(iDimZ/real(iPc, dp)), i8b))*iPc; 
            else

                iDimX = cSpace%iDimX; 
                iDimY = cSpace%iDimY; 
                iDimZ = cSpace%iDimZ; 
!KH write (acTemp, '("iDimX, iDimY, iDimZ after:    ",3I4)') iDimX,iDimY,iDimZ
!KH call mexPrintf(acTemp//achar(13))

            end if

        else

            ! If we are not allowed to change any dimension, we assume that
            ! the dimensions are correct. If not, issue a warning, but continue running
            iDimX = cSpace%iDimX; 
            iDimY = cSpace%iDimY; 
            iDimZ = cSpace%iDimZ; 
            if (.NOT. (isProdSmallPrimes(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ) .AND. &
                       (mod(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ, cSpace%cGrid%iProcN) == 0))) then
                write (acTemp, "('Warning while initializing the grid: no dimension was allowed to be changed.')")
                call PrintToLog(acTemp, -1)
                write (acTemp, "('However, the current choice for iDimX, iDimY, iDimZ,  is not correct.')")
                call PrintToLog(acTemp, -1)
                write (acTemp, "('(iDimX, iDimY, iDimZ) = (', I5, ', ', I5, ', ', I5, ') the number of processors is: ', I4, '.')") iDimX, iDimY, iDimZ, cSpace%cGrid%iProcN
                call PrintToLog(acTemp, -1)
            end if

        end if

        write (acTemp, "('Original problem size: ', I4, ' * ', I4, ' * ', I4, ' * ', I4)") cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ; 
        call PrintToLog(acTemp, 3); 
        write (acTemp, "('The resized problem has dimensions: ', I4, ' * ', I4, ' * ', I4, ' * ', I4)") iDimT, iDimX, iDimY, iDimZ; 
        call PrintToLog(acTemp, 3); 
        ! Now we have obtained values which satisfy our demands, we can calculate the new start
        ! and end-point of each dimension, but now with good values for the dimensions.
        call InitSpace(cSpace, cSpace%iSpaceIdentifier, cSpace%bYSymm, &
                       iDimT, iDimX, iDimY, iDimZ, &
                       cSpace.iStartT, cSpace.iStartX, cSpace.iStartY, cSpace.iStartZ, &
                       cSpace.iBeamIndexStartX, cSpace.iBeamIndexStartY, cSpace.iBeamIndexStartZ, &
                       cSpace%dFnyq, cSpace.dTanX, cSpace.dTanY, cSpace.dTanT)

    END SUBROUTINE InitGridSpace

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %                InitGridDistrLookupTables
! %
! % This function initializes all information, concerning the number of items we store locally, the size
! %  of this data etc etc. It will also calculate which t-blocks are stored locally in Distribution 0 and 1
! %  and which xyz blocks are stored locally in distribution 2. (Note that we distribute both the t and xyz
! %  blocks over the different processors, using a block-distribution)
! %
! %  Input:
! %          cGrid- a Grid for which the associated Space is already initialized.
! %
! %  Output:
! %          cGrid- a Grid for which all the information concerning the distributions is computed
! %
! %
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    SUBROUTINE InitGridDistrLookupTables(cSpace)

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
!   The subroutine InitGridDistrLookupTables initializes the parameters within
!   the Grid structure concerning the dimension sizes, array strides and beam
!   indices of the data array stored in the form of each of the three
!   distributions.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace     io   type(Space)   Space on which InitGridSpace operates.
!
        type(Space), INTENT(inout) ::                cSpace; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i        i8b    loop counter
!   j        i8b    loop counter
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::                                i, j
        character(LEN=2048) ::                acTemp; 
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
!
! =============================================================================

        write (acTemp, "('InitGridDistrLookupTables')"); call PrintToLog(acTemp, 3); 
        ! Before we start, we will check whether the lookup-arrays are still not-allocated. If this is
        !  not the case, we have to deallocate them ourselves.
        if (allocated(cSpace%cGrid%aiD0Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD0Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD0Loc); 
        end if
        if (allocated(cSpace%cGrid%aiD1Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD1Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD1Loc); 
        end if
        if (allocated(cSpace%cGrid%aiD2Loc)) then
            iMemAllocated = iMemAllocated - (PRODUCT(SHAPE(cSpace%cGrid%aiD2Loc)))*i4bS
            deallocate (cSpace%cGrid%aiD2Loc); 
        end if

        ! We will start the initialization for distribution 0 and continue till distribution 2

        ! First we initialize all the values for Distribution 1, note that it is block-distributed
        cSpace%cGrid%iD0TL = cSpace%iDimT; 
        cSpace%cGrid%iD0GlobN = cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ; 
        cSpace%cGrid%iD0LocN = cSpace%cGrid%iD0GlobN/cSpace%cGrid%iProcN; 
        cSpace%cGrid%iD0LocSize = cSpace%cGrid%iD0LocN*cSpace%cGrid%iD0TL; 
        cSpace%cGrid%iD0TS = 1; 
        cSpace%cGrid%iD0IS = cSpace%cGrid%iD0TL; 
#ifndef RUNESTIMATOR
        !        We initialize the lookuptable, which contains the indices of the t-blocks present
        !        locally on this processor
        allocate (cSpace%cGrid%aiD0Loc(cSpace%cGrid%iD0LocN, 3))
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(cSpace%cGrid%aiD0Loc)))*i4bS
        do i = 0, cSpace%cGrid%iD0LocN - 1
            ! We have i = x + y * iDimX + z * iDimX * iDimY
            j = i + (cSpace%cGrid%iProcID)*cSpace%cGrid%iD0LocN; !THIS CHANGED BY A.M 21/09/2020
            cSpace%cGrid%aiD0Loc(1 + i, 1) = mod(j, cSpace%iDimX); ; 
            j = (j - mod(j, cSpace%iDimX))/cSpace%iDimX; 
            cSpace%cGrid%aiD0Loc(1 + i, 2) = mod(j, cSpace%iDimY); 
            j = (j - mod(j, cSpace%iDimY))/cSpace%iDimY; 
            cSpace%cGrid%aiD0Loc(1 + i, 3) = j; 
        end do
        
#endif

        ! We initialize the lookup-table for Distribution 1, note that it is block-distributed
!KH r2c stuff: is it always plus 1, also with odd iDimT?? no, don"t think so...
!KH Does this assume a double length T axis already? yes, do think so...
!KH Change for periodical T
        cSpace%cGrid%iD1TL = cSpace%iDimT + 1; ! Added by A.M 28/06/2020
        cSpace%cGrid%iD1GlobN = cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ; 
        cSpace%cGrid%iD1LocN = cSpace%cGrid%iD1GlobN/cSpace%cGrid%iProcN; 
        cSpace%cGrid%iD1LocSize = cSpace%cGrid%iD1LocN*cSpace%cGrid%iD1TL; 
        cSpace%cGrid%iD1TS = 1; 
        cSpace%cGrid%iD1IS = cSpace%cGrid%iD1TL; 
#ifndef RUNESTIMATOR
        !         We initialize the lookuptable, which contains the indices of the t-blocks present
        !        locally on this processor
        allocate (cSpace%cGrid%aiD1Loc(cSpace%cGrid%iD1LocN, 3))
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(cSpace%cGrid%aiD1Loc)))*i4bS
        do i = 0, cSpace%cGrid%iD1LocN - 1
            ! We have i = x + y * iDimX + z * iDimX * iDimY
            j = i + (cSpace%cGrid%iProcID)*cSpace%cGrid%iD0LocN; 
            cSpace%cGrid%aiD1Loc(1 + i, 1) = mod(j, cSpace%iDimX); ; 
            j = (j - mod(j, cSpace%iDimX))/cSpace%iDimX; 
            cSpace%cGrid%aiD1Loc(1 + i, 2) = mod(j, cSpace%iDimY); 
            j = (j - mod(j, cSpace%iDimY))/cSpace%iDimY; 
            cSpace%cGrid%aiD1Loc(1 + i, 3) = j; 
        end do
#endif

        ! We initialize the lookup-table for Distribution 1, note that it is block-distributed
        cSpace%cGrid%iD2XL = cSpace%iDimX; 
        cSpace%cGrid%iD2YL = cSpace%iDimY; 
        cSpace%cGrid%iD2ZL = cSpace%iDimZ; 
        cSpace%cGrid%iD2GlobN = cSpace%iDimT + 1; ! Added by A.M 28/06/2020
        cSpace%cGrid%iD2LocN = cSpace%cGrid%iD2GlobN/cSpace%cGrid%iProcN; 
        cSpace%cGrid%iD2LocSize = cSpace%cGrid%iD2LocN*cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ; 
        !KH We have i = x + y * iDimX + z * iDimX * iDimY
        cSpace%cGrid%iD2XS = 1; 
        cSpace%cGrid%iD2YS = cSpace%iDimX; 
        cSpace%cGrid%iD2ZS = cSpace%iDimX*cSpace%iDimY; 
        cSpace%cGrid%iD2IS = cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ; 
#ifndef RUNESTIMATOR
        ! We initialize the lookuptable, which contains the indices of the t-blocks present
        ! locally on this processor
        allocate (cSpace%cGrid%aiD2Loc(cSpace%cGrid%iD2LocN)); 
        iMemAllocated = iMemAllocated + (PRODUCT(SHAPE(cSpace%cGrid%aiD2Loc)))*i4bS
        do i = 0, cSpace%cGrid%iD2LocN - 1
            cSpace%cGrid%aiD2Loc(1 + i) = i + cSpace%cGrid%iProcID*cSpace%cGrid%iD2LocN; 
        end do
#endif

    END SUBROUTINE InitGridDistrLookupTables

    INTEGER(i8b) FUNCTION iBeamOffsetT(cSpace, iBeamIndexZ)

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
!   The function iBeamOffsetT returns the offset in t as a function of the beam
!   index in z. With this offset, a beam index may be translated to a global
!   index and vice versa. Relevant for skew beams only. See also the
!   ParnacDataDef module.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iBeamOffsetT  o   i8b         BeamOffset in T
!   cSpace        i   type(Space) Space for which the BeamOffset needs to be
!                                 given
!   iBeamIndexZ   i   i8b         Beam Index in z for which the BeamOffset
!                                 needs to be given
!
        type(Space), INTENT(in) ::    cSpace; 
        integer(i8b), INTENT(in) :: iBeamIndexZ

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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

        iBeamOffsetT = nint((iBeamIndexZ + cSpace%iStartZ)*cSpace%dTanT) - nint(cSpace%iStartZ*cSpace%dTanT); 
    END FUNCTION iBeamOffsetT

    INTEGER(i8b) FUNCTION iBeamOffsetX(cSpace, iBeamIndexZ)

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
!   The function iBeamOffsetX returns the offset in x as a function of the beam
!   index in z. With this offset, a beam index may be translated to a global
!   index and vice versa. Relevant for skew beams only. See also the
!   ParnacDataDef module.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iBeamOffsetX  o   i8b         BeamOffset in X
!   cSpace        i   type(Space) Space for which the BeamOffset needs to be
!                                 given
!   iBeamIndexZ   i   i8b         Beam Index in z for which the BeamOffset
!                                 needs to be given
!
        type(Space), INTENT(in) ::    cSpace; 
        integer(i8b), INTENT(in) :: iBeamIndexZ

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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

        iBeamOffsetX = nint((iBeamIndexZ + cSpace%iStartZ)*cSpace%dTanX) - nint(cSpace%iStartZ*cSpace%dTanX); 
    END FUNCTION iBeamOffsetX

    INTEGER(i8b) FUNCTION iBeamOffsetY(cSpace, iBeamIndexZ)

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
!   The function iBeamOffsetY returns the offset in y as a function of the beam
!   index in z. With this offset, a beam index may be translated to a global
!   index and vice versa. Relevant for skew beams only. See also the
!   ParnacDataDef module.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iBeamOffsetY  o   i8b         BeamOffset in Y
!   cSpace        i   type(Space) Space for which the BeamOffset needs to be
!                                 given
!   iBeamIndexZ   i   i8b         Beam Index in z for which the BeamOffset
!                                 needs to be given
!
        type(Space), INTENT(in) ::    cSpace; 
        integer(i8b), INTENT(in) :: iBeamIndexZ

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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

        iBeamOffsetY = nint((iBeamIndexZ + cSpace%iStartZ)*cSpace%dTanY) - nint(cSpace%iStartZ*cSpace%dTanY); 
    END FUNCTION iBeamOffsetY

    FUNCTION ObtainNextSmallPrime(iN)

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
!   The function ObtainNextSmallPrime returns the product of small primes 2, 3,
!   5 and 7 that is equal or next larger than the number iN.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   ObtainNextSmallPrime   o   i8b   Next product of small primes
!   iN                     i   i8b   Number for which the next product of small
!                                    primes needs to be determined
!
        integer(i8b) :: ObtainNextSmallPrime
        integer(i8b), INTENT(in) ::                 iN; 
! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i   i8b   temporary variable
!
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
!   IsProdSmallPrimes
!
! =============================================================================

        i = iN; 
        do while (IsProdSmallPrimes(i) .EQV. .false.)
            i = i + 1; 
        end do

        ObtainNExtSmallPrime = i; 
    END FUNCTION ObtainNextSmallPrime

    FUNCTION Isprodsmallprimes(Number)

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
!   The function Isprodsmallprimes tests whether a number is a product of
!   small primes (2,3,5 and 7)
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   isProdSmallPrimes  o   lgt    .true. if Number is a product of small primes
!                                 .false. if not
!   Number             i   i8b    number which is tested
!
        logical(lgt) :: Isprodsmallprimes
        integer(i8b), intent(in) :: Number

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   Nhelp   i8b   Number divided over and over by the small primes
!   loop    i8b   loop counter
!   primes  i8b   array containing the small primes
!
        integer(i8b) :: Nhelp, loop
        integer(i8b), dimension(4) :: primes = (/2, 3, 5, 7/)

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

        Nhelp = Number

        if (Nhelp == 0) then
            Isprodsmallprimes = .false.
            return
        end if

        do loop = 1, size(primes)
            do
                if (mod(Nhelp, primes(loop)) == 0) then
                    Nhelp = Nhelp/primes(loop)
                else
                    exit
                end if
            end do
        end do

        if (Nhelp == 1) then
            Isprodsmallprimes = .true.
        else
            Isprodsmallprimes = .false.
        end if

    END FUNCTION Isprodsmallprimes

    SUBROUTINE InitSaveSlices(cSpace)

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
!   The subroutine InitSaveSlices initializes the Save slices variables. It
!   translates the x/y/z save positions of the slices that need to be output
!   to file to beam numbers and beam indices.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   i   type(Space)  Space for which the Save Slices are determined
!
        type(Space), intent(in) :: cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i           i8b    loop counter
!   iBeamIndex  i8b    beam index for a save slice coordinate
!   err         i8b    error number
!   acTemp      char   temporary char array for output log messages
!
        integer(i8b) :: i
        integer(i8b) :: iBeamIndex
        character(LEN=1024) :: acTemp

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
!   none
!
! =============================================================================

        call PrintToLog("InitializeSlices", 1)

        !Translate the x/y/z save positions to indices
        ! t/x/y: all beams (set cModelParams%beam = -1)
        ! z: special case: only located in one specific beam

        do i = 1, cModelParams%numslices
            if (cModelParams%xyzslicedim(i) == 't') then
                iBeamIndex = nint(cModelParams%xyzslicepos(i)/cSpace%dDt)
                cModelParams%xyzsliceindex(i) = iBeamIndex - cSpace%iStartT
                cModelParams%xyzslicebeam(i) = -1
            elseif (cModelParams%xyzslicedim(i) == 'x') then
                iBeamIndex = nint(cModelParams%xyzslicepos(i)/cSpace%dDx)
                cModelParams%xyzsliceindex(i) = iBeamIndex - cSpace%iStartx
                cModelParams%xyzslicebeam(i) = -1
            elseif (cModelParams%xyzslicedim(i) == 'y') then
                if (cSpace%bYSymm) then
                    cModelParams%xyzslicepos(i) = abs(cModelParams%xyzslicepos(i))
                end if
                iBeamIndex = nint(cModelParams%xyzslicepos(i)/cSpace%dDx)
                cModelParams%xyzsliceindex(i) = iBeamIndex - cSpace%iStarty
                cModelParams%xyzslicebeam(i) = -1
            elseif (cModelParams%xyzslicedim(i) == 'z') then
                iBeamIndex = nint(cModelParams%xyzslicepos(i)/cSpace%dDx)
                cModelParams%xyzslicebeam(i) = int((iBeamIndex - cSpace%iStartz)/cSpace%iDimZ)
                cModelParams%xyzsliceindex(i) = iBeamIndex - cSpace%iStartz - cSpace%iDimZ*cModelParams%xyzslicebeam(i)
            end if
        end do

    END SUBROUTINE InitSaveSlices

    SUBROUTINE InitRRMSNorm(cRRMSNorm, iProcID, iNumBeams, iNumIterations)

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
!   The subroutine InitRRMSNorm initializes the RRMSNorm structure.
!   Memory is allocated, and all is set to zero.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cRRMSNorm   io   type(RRMSNormType)  RRMSNorm to be initialized
!   iProcID     i    i8b   ID number of the current processor
!   iNumBeams   i    i8b   Number of beams
!   iNumIterations   i    i8b   Total number of iterations
!
        type(RRMSNormtype), intent(inout) :: cRRMSNorm
        integer(i8b) :: iProcID, iNumBeams, iNumIterations

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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
!   none
!
! =============================================================================

        call PrintToLog("InitRRMSNorm", 3)

        if (iProcID == 0) then
            allocate (cRRMSNorm%arRRMSNorm(1:iNumBeams, 1:iNumIterations - 1))

            cRRMSNorm%arRRMSNorm = 0.0_dp

        end if

    END SUBROUTINE InitRRMSNorm

    SUBROUTINE InitRRMSPreviousCorrection(cRRMSNorm, cSpace)

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
!   The subroutine InitRRMSPreviousCorrection initializes the
!   PreviousCorrection data array within the RRMSNorm structure.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cRRMSNorm   io   type(RRMSNormType)  RRMSNorm structure in which the
!                                        PreviousCorrection array is to be
!                                        initialized
!   cSpace      i    type(Space)         Space for which the RRMSNorm is to be
!                                        determined
!
        type(RRMSNormtype), INTENT(inout) :: cRRMSNorm
        type(Space), INTENT(in) :: cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi   i8b   loop counter over all xyz-positions in the beam
!
        integer(i8b) :: iLi

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
!   none
!
! =============================================================================

        call PrintToLog("InitRRMSPreviousCorrection", 3)

        !check how many points there are located in this beam
        !that lie on the slice at the last z position
        !search for the first xy position and include rest of the trace in
        !the RRMS calculation (assuming that the z dimension is the largest
        !running index as in distribution 0).
        do iLi = 0, cSpace.cGrid.iD0LocN - 1
            if (cSpace.cGrid.aiD0Loc(1 + iLi, 3) == cSpace.iDimZ - 1) then
                exit
            end if
        end do

        if (iLi < cSpace.cGrid.iD0LocN) then
            !set cRRMSNorm%iRRMSstart, cRRMSNorm%iRRMSlength
            cRRMSNorm%iStart = iLi*cSpace%cGrid%iD0TL
            !allocate and initialize cRRMSNorm%arPreviouscorrection
            allocate (cRRMSNorm%arPreviouscorrection(1:(cSpace.cGrid.iD0LocN - iLi)*cSpace%cGrid%iD0TL))
            cRRMSNorm%arPreviouscorrection = 0.0_dp
        else
            cRRMSNorm%iStart = -1
        end if

    END SUBROUTINE InitRRMSPreviousCorrection

    SUBROUTINE DestroyRRMSPreviousCorrection(cRRMSNorm)

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
!   The subroutine DestroyRRMSPreviousCorrection destructs the
!   PreviousCorrection data array within the RRMSNorm structure.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cRRMSNorm   io   type(RRMSNormType)  RRMSNorm structure in which the
!                                        PreviousCorrection array is destructed
!
        type(RRMSNormtype), INTENT(inout) :: cRRMSNorm

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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
!   none
!
! =============================================================================

        call PrintToLog("DestroyRRMSPreviousCorrection", 3)

        if (allocated(cRRMSNorm%arPreviouscorrection)) then
            deallocate (cRRMSNorm%arPreviouscorrection)
        end if

    END SUBROUTINE DestroyRRMSPreviousCorrection

    SUBROUTINE DestroyRRMSNorm(cRRMSNorm)

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
!   The subroutine DestroyRRMSNorm destructs the RRMSNorm structure.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cRRMSNorm   io   type(RRMSNormType)  RRMSNorm to be destructed
!
        type(RRMSNormtype), INTENT(inout) :: cRRMSNorm

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   none
!
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
!   none
!
! =============================================================================

        call PrintToLog("DestroyRRMSNorm", 3)

        if (allocated(cRRMSNorm%arRRMSNorm)) then
            deallocate (cRRMSNorm%arRRMSNorm)
        end if
        if (allocated(cRRMSNorm%arPreviouscorrection)) then
            deallocate (cRRMSNorm%arPreviouscorrection)
        end if

    END SUBROUTINE DestroyRRMSNorm

END MODULE ParnacDataSpaceInit
