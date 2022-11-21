MODULE ParnacDataMap

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The module ParnacDataMap contains functions to map and inverse map the data
!   inside a beam to include the offsets due to the skewness of a beam in the
!   T, X and/or Y dimensions. The offsets are applied to shift the data, and if
!   the data then overflows the upper boundary, it is wrapped and shows up at
!   the lower boundary. (idem for overflow at the lower boundary). For the
!   mappings in the different dimensions we first define generic forward and
!   inverse mapping functions for an arbitrary dimension and for real or
!   complex data.
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types
    USE Constants
    USE ParnacGeneral
    USE ParnacDataDef

    USE ParnacDataSpaceInit, ONLY: &
        iBeamOffsetT, iBeamOffsetX, iBeamOffsetY

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
!   MapVReal            sub   Map a space of real numbers in one dimension
!   MapVInvReal         sub   Inverse map a space of real numbers in one dim.
!   MapVComplex         sub   Map a space of complex numbers in one dimension
!   MapVInvComplex      sub   Inverse map a space of complex numbers in one dim
!   MapVRealComplex     sub   Map a space of real numbers, stored in complex
!                               data format, in one dimension
!   MapVInvRealComplex  sub   Inverse map a space of real numbers, stored in
!                               complex data format, in one dimension
!   MapVt               sub   Map a space in its temporal dimension
!   MapVtInv            sub   Inverse map a space in its temporal dimension
!   MapVx               sub   Map a space in its x-dimension
!   MapVy               sub   Map a space in its y-dimension
!   MapVx_Small         sub   Map a space in its x-dimension over half
!                             the dimension length
!   MapVy_Small         sub   Map a space in its y-dimension over half
!                             the dimension length
!   MapVxInv            sub   Inverse map a space in its x-dimension
!   MapVyInv            sub   Inverse map a space in its y-dimension
!   MapVxInv_Small      sub   Inverse map a space in its x-dimension over half
!                             the dimension length
!   MapVyInv_Small      sub   Inverse map a space in its y-dimension over half
!                             the dimension length
!
! =============================================================================
!
CONTAINS

    SUBROUTINE MapVReal(arData, arWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVReal applies a mapping to the data in a beam to account
!   for the offsets due to the skewness of a beam in a certain dimension. The
!   offsets are applied to shift the data over a number of positions, and if
!   the data then overflows the upper boundary, it is wrapped and shows up at
!   the lower boundary. (idem for overflow at the lower boundary for a negative
!   shift). The routine wraps the dimension that is the first one in the
!   aiBlockL and aiBlockS arrays; this dimension is here referred to as the
!   row.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   arData    i  dp    The original data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data
!   iDim      i  i8b   number of dimensions of the data in arData and arWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions in arData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides of the dimensions in arData
!                      The first dimension is the dimension to be wrapped.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in arWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   arWrap    o  dp    The wrapped data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        real(dp), intent(in)::       arData(:)
        integer(i8b), intent(in)::   iDim
        integer(i8b), intent(in)::   aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::   iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::   aiOffsets(:)
        real(dp), intent(out)::      arWrap(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   arBuffer   dp   temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        real(dp), allocatable::         arBuffer(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! *****************************************************************************
!
!   PSEUDOCODE
!
!   initialize buffer, loop variables and indices
!   for each row
!       update array indices for arData and arWrap to start of the current row
!       update aiLoop indices to the current row
!       get offset for this row and
!       obtain the modulo of the offset with the wrap length
!       copy the current row from arData to a buffer
!       copy first part of the buffer to the shifted position in the row in
!         arWrap
!       copy second part of the buffer to the start of the row in arWrap
!   deinitialize buffer, loop variable
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVReal", 6)

        ! Initialize the buffer
        call PrintToLog("Initialize MapVReal", 7)
        allocate (arBuffer(aiBlockL(1))); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 7)
        do iLi = 0, iIndMax - 1

            ! Update the location indices to the current row in arData and arWrap,
            ! and update the loop index variable aiLoop(2:4) to the current position.
            ! If the index in aiLoop(2:4) overruns its dimension length, reset it to
            ! 0, increase the dimension index iLD, increase the index of the next
            ! higher dimension aiLoop(iLD) and adjust the indices iIndData and
            ! iIndWrap accordingly. In other words, this results in a chique multiple
            ! do loop, without knowing the number of dimensions in advance:
            ! do aiLoop(2) = 0, aiBlockL(2)-1
            !     iIndData = iIndData + aiBlockS(2)
            !     iIndWrap = iIndWrap + aiWrapS(2)
            !     do aiLoop(3) = 0, aiBlockL(3)-1
            !         iIndData = iIndData + aiBlockS(3)
            !         iIndWrap = iIndWrap + aiWrapS(3)
            !
            !         ...
            !
            !     end do
            ! end do
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            ! Store the block temporarily in a buffer-variable
            arBuffer(:) = arData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)); 
            arWrap(iIndWrap:iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)) = 0; 
            !Put the elements back, but now mapped
            !KH The condition is never met, unless iOffset=0 and aiBlockL(1) = iWrapL
            if (iOffset + aiBlockL(1) <= iWrapL) then
                arWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iOffset + aiBlockL(1) - 1)*aiWrapS(1):aiWrapS(1)) = arBuffer(:); 
            else
                iBreak = iWrapL - iOffset; 
                arWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)) = arBuffer(1:iBreak); 
                arWrap(iIndWrap:iIndWrap + (aiBlockL(1) - iBreak - 1)*aiWrapS(1):aiWrapS(1)) = arBuffer(iBreak + 1:aiBlockL(1))
            end if
        end do

        ! Stopwatches, stop
        deallocate (arBuffer); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVReal

    SUBROUTINE MapVInvReal(arData, arWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVInvReal applies the inverse mapping of MapVReal.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   arWrap    i  dp    The wrapped data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data.
!   iDim      i  i8b   number of dimensions of the data in arData and arWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions in arData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides of the dimensions in arData
!                      The first dimension is the dimension to be wrapped.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in arWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   arData    o  dp    The original data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        real(dp), intent(in)::    arWrap(:)
        integer(i8b), intent(in)::      iDim
        integer(i8b), intent(in)::      aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::      iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::      aiOffsets(:)
        real(dp), intent(out)::   arData(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   arBuffer   dp   temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        real(dp), allocatable::         arBuffer(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVInvReal", 6)

        ! Initialize a handy pointer and the buffer
        call PrintToLog("Initialize MapVInvReal", 7)
        allocate (arBuffer(aiBlockL(1))); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 7)
        do iLi = 0, iIndMax - 1
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            ! Store the block temporary in a buffer-variable
            if (iOffset + aiBlockL(1) <= iWrapL) then
                arBuffer(:) = arWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + iOffset + (aiBlockL(1) - 1)*aiWrapS(1):aiWrapS(1)); 
            else
                iBreak = iWrapL - iOffset; 
                arBuffer(1:iBreak) = arWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)); 
                arBuffer(iBreak + 1:aiBlockL(1)) = arWrap(iIndWrap:iIndWrap + (aiBlockL(1) - iBreak - 1)*aiWrapS(1):aiWrapS(1)); 
            end if

            ! Put the elements back, but now unmapped
            arData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)) = 0; 
            arData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)) = arBuffer(:); 
        end do

        ! Stopwatches, stop
        deallocate (arBuffer); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVInvReal

    SUBROUTINE MapVComplex(acData, acWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVComplex performs the same operation as MapVReal, but
!   on a data array containing a complex variable.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   acData    i  dpc   The original data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data
!   iDim      i  i8b   number of dimensions of the data in acData and acWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions in acData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides of the dimensions in acData
!                      The first dimension is the dimension to be wrapped.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in acWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   acWrap    o  dpc   The wrapped data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        complex(dpc), intent(in)::       acData(:)
        integer(i8b), intent(in)::   iDim
        integer(i8b), intent(in)::   aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::   iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::   aiOffsets(:)
        complex(dpc), intent(out)::      acWrap(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   acBuffer  dpc   temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        complex(dpc), allocatable::       acBuffer(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVComplex", 6)

        ! Initialize a handy pointer and the buffer
        call PrintToLog("Initialize MapVComplex", 7)
        allocate (acBuffer(aiBlockL(1))); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 7)
        do iLi = 0, iIndMax - 1
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo iWrapL
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            ! Store the block temporary in a buffer-variable
            acBuffer(:) = acData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)); 
            acWrap(iIndWrap:iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)) = 0; 
            !Put the elements back, but now mapped
            if (iOffset + aiBlockL(1) <= iWrapL) then
                acWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iOffset + aiBlockL(1) - 1)*aiWrapS(1):aiWrapS(1)) = acBuffer(:); 
            else
                iBreak = iWrapL - iOffset; 
                acWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)) = acBuffer(1:iBreak); 
                acWrap(iIndWrap:iIndWrap + (aiBlockL(1) - iBreak - 1)*aiWrapS(1):aiWrapS(1)) = acBuffer(iBreak + 1:aiBlockL(1))
            end if
        end do

        ! Stopwatches, stop
        deallocate (acBuffer); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVComplex

    SUBROUTINE MapVInvComplex(acData, acWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVInvComplex applies the inverse mapping of MapVComplex.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   acWrap    i  dpc   The wrapped data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data.
!   iDim      i  i8b   number of dimensions of the data in acData and acWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions in acData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides of the dimensions in acData
!                      The first dimension is the dimension to be wrapped.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in acWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   acData    o  dpc   The original data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        complex(dpc), intent(in)::    acWrap(:)
        integer(i8b), intent(in)::      iDim
        integer(i8b), intent(in)::      aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::      iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::      aiOffsets(:)
        complex(dpc), intent(out)::   acData(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   acBuffer  dpc   temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        complex(dpc), allocatable::         acBuffer(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVInvComplex", 6)

        ! Initialize a handy pointer and the buffer
        call PrintToLog("Initialize MapVinvComplex", 7)
        allocate (acBuffer(aiBlockL(1))); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 6)
        do iLi = 0, iIndMax - 1
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo iWrapL
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            ! Store the block temporary in a buffer-variable
            if (iOffset + aiBlockL(1) <= iWrapL) then
                acBuffer(:) = acWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + iOffset + (aiBlockL(1) - 1)*aiWrapS(1):aiWrapS(1)); 
            else
                iBreak = iWrapL - iOffset; 
                acBuffer(1:iBreak) = acWrap(iIndWrap + iOffset*aiWrapS(1):iIndWrap + (iWrapL - 1)*aiWrapS(1):aiWrapS(1)); 
                acBuffer(iBreak + 1:aiBlockL(1)) = acWrap(iIndWrap:iIndWrap + (aiBlockL(1) - iBreak - 1)*aiWrapS(1):aiWrapS(1)); 
            end if

            ! Put the elements back, but now mapped
            acData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)) = 0; 
            acData(iIndData:iIndData + (aiBlockL(1) - 1)*aiBlockS(1):aiBlockS(1)) = acBuffer(:); 
        end do

        ! Stopwatches, stop
        deallocate (acBuffer); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVInvComplex

    SUBROUTINE MapVRealComplex(acData, acWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVRealComplex performs the same operation as MapVReal, but
!   on a data array containing a real variable stored in the complex datatype,
!   i.e. [arData(1) arData(2) arData(3) ...] = [real(acData(1)) imag(acData(1))
!   real(acData(2)) ...]. The length parameters aiBlockL and iWrapL refer to
!   the lengths of the REAL variable, the stride parameters aiBlockS and
!   aiWrapS refer to the strides in terms of the COMPLEX variable.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   acData    i  dpc   The original data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data.
!   iDim      i  i8b   number of dimensions of the data in acData and acWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions of the real
!                      variable stored inside acData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides for each of the dimensions
!                      in acData. The first dimension is the dimension to be
!                      wrapped. The strides are in terms of the complex
!                      datatype.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in acWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   acWrap    o  dpc   The wrapped data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        complex(dpc), intent(in)::   acData(:)
        integer(i8b), intent(in)::   iDim
        integer(i8b), intent(in)::   aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::   iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::   aiOffsets(:)
        complex(dpc), intent(out)::  acWrap(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   arBufferIn  dp  temporary buffer for the shift
!   arBufferOut dp  temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        real(dp), allocatable::         arBufferIn(:), arBufferOut(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVRealComplex", 6)

        ! Initialize a handy pointer and the buffer
        call PrintToLog("Initialize MapVRealComplex", 7)
        allocate (arBufferIn(aiBlockL(1)), arBufferOut(iWrapL)); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 7)
        do iLi = 0, iIndMax - 1
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            !fill real input buffer from complex array
            if (mod(aiBlockL(1), 2_i8b) == 1) then
                arBufferIn(1:aiBlockL(1):2) = real(acData(iIndData:iIndData + ((aiBlockL(1) + 1)/2 - 1)*aiBlockS(1):aiBlockS(1)), dp); 
            else
                arBufferIn(1:aiBlockL(1):2) = real(acData(iIndData:iIndData + (aiBlockL(1)/2 - 1)*aiBlockS(1):aiBlockS(1)), dp); 
            end if
            arBufferIn(2:aiBlockL(1):2) = dimag(acData(iIndData:iIndData + (aiBlockL(1)/2 - 1)*aiBlockS(1):aiBlockS(1))); 
            !clear output buffer
            arBufferOut = 0; 
            !Put the elements back, but now mapped
            !KH With the changes I made, the condition is never met...
            if (iOffset + aiBlockL(1) <= iWrapL) then
                arBufferOut(1 + iOffset:aiBlockL(1) + iOffset) = arBufferIn
            else
                iBreak = iWrapL - iOffset; 
                arBufferOut(1 + iOffset:iWrapL) = arBufferIn(1:iBreak); 
                arBufferOut(1:aiBlockL(1) - iBreak) = arBufferIn(iBreak + 1:aiBlockL(1)); 
            end if

            !store real output buffer in complex array
            if (mod(iWrapL, 2_i8b) == 1) then
                acWrap(iIndWrap:iIndWrap + (iWrapL/2 - 1)*aiWrapS(1):aiWrapS(1)) = arBufferOut(1:iWrapL - 1:2) + im*arBufferOut(2:iWrapL - 1:2); 
                acWrap(iIndWrap + ((iWrapL + 1)/2 - 1)*aiWrapS(1)) = arBufferOut(iWrapL); 
            else
                acWrap(iIndWrap:iIndWrap + (iWrapL/2 - 1)*aiWrapS(1):aiWrapS(1)) = arBufferOut(1:iWrapL:2) + im*arBufferOut(2:iWrapL:2); 
            end if

        end do

        ! Stopwatches, stop
        deallocate (arBufferIn, arBufferOut); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVRealComplex

    SUBROUTINE MapVInvRealComplex(acData, acWrap, iDim, aiBlockL, aiBlockS, iWrapL, aiWrapS, aiOffsets)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVInvRealComplex applies the inverse mapping of
!   MapVRealComplex.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   acWrap    i  dpc   The wrapped data, passed as a 1D array. The parameters
!                      iDim, aiBlockL and aiBlockS determine the shape of the
!                      data.
!   iDim      i  i8b   number of dimensions of the data in acData and acWrap
!   aiBlockL  i  i8b   array with the lengths of the dimensions of the real
!                      variable stored inside acData.
!                      The first dimension is the dimension to be wrapped.
!   aiBlockS  i  i8b   array with the strides for each of the dimensions
!                      in acData. The first dimension is the dimension to be
!                      wrapped. The strides are in terms of the complex
!                      datatype.
!   iWrapL    i  i8b   The wrapping length over which the first dimension
!                      is wrapped.
!   aiWrapS   i  i8b   The stride of each dimension in acWrap
!   aiOffsets i  i8b   The offset of each of the the rows.
!                      The length of aiOffsets is PRODUCT(iDim(2:))
!   acData    o  dpc   The original data, passed as a 1D array. In the calling
!                      environment, this array may refer to the same array ,
!                      i.e. this implementation supports in-place evaluation.
!
        complex(dpc), intent(in)::    acWrap(:)
        integer(i8b), intent(in)::      iDim
        integer(i8b), intent(in)::      aiBlockL(:), aiBlockS(:)
        integer(i8b), intent(in)::      iWrapL, aiWrapS(:)
        integer(i8b), intent(in)::      aiOffsets(:)
        complex(dpc), intent(out)::   acData(:)

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   iLi       i8b   loop counter, temporary or pointing to the rows
!   lLD       i8b   loop counter, pointing to the dimensions
!   iOffset   i8b   Offset modulo iWrapL
!   iBreak    i8b   Index where to break the current row, shift the first
!                   part to a higher index and the rest of the row to index 0.
!   iIndData  i8b   Index to the current location in arData
!   iIndWrap  i8b   Index to the current location in arWrap
!   iIndMax   i8b   Total number of rows in arData
!   aiLoop    i8b   1D array of size iDim, running index over the positions
!                   in the dimensions 2:iDim. aiLoop(1) is not used.
!   arBufferIn  dp  temporary buffer for the shift
!   arBufferOut dp  temporary buffer for the shift
!   acTemp   char   temporary char array for output log messages
!
        integer(i8b) ::            iLi, iLD, iOffset, iBreak, iIndData, iIndWrap, iIndMax
        integer(i8b), allocatable ::      aiLoop(:); 
        real(dp), allocatable::         arBufferIn(:), arBufferOut(:); 
        CHARACTER(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!
! =============================================================================

        ! Stopwatches, start
        call PrintToLog("Start MapVInvRealComplex", 6)

        ! Initialize a handy pointer and the buffer
        call PrintToLog("Initialize MapVInvRealComplex", 7)
        allocate (arBufferIn(iWrapL), arBufferOut(aiBlockL(1))); 
        ! Initialize the loopvariables
        allocate (aiLoop(iDim)); 
        do iLi = 1, iDim
            aiLoop(iLi) = 0; 
        end do
        aiLoop(2) = -1; 
        iIndData = 1 - aiBlockS(2); 
        iIndWrap = 1 - aiWrapS(2); 
        iIndMax = 1; 
        do iLi = 2, iDim
            iIndMax = iIndMax*aiBlockL(iLi); 
        end do

        ! Apply the map for each block
        call PrintToLog("Start the loop", 7)
        do iLi = 0, iIndMax - 1
            iIndData = iIndData + aiBlockS(2); 
            iIndWrap = iIndWrap + aiWrapS(2); 
            aiLoop(2) = aiLoop(2) + 1; 
            iLD = 2; 
            do while (aiLoop(iLD) >= aiBlockL(iLD))
                iIndData = iIndData - aiBlockS(iLD)*aiBlockL(iLD); 
                iIndWrap = iIndWrap - aiWrapS(iLD)*aiBlockL(iLD); 
                aiLoop(iLD) = 0; 
                iLD = iLD + 1; 
                aiLoop(iLD) = aiLoop(iLD) + 1; 
                iIndData = iIndData + aiBlockS(iLD); 
                iIndWrap = iIndWrap + aiWrapS(iLD); 
            end do

            ! Obtain the absolute offset of this t-block
            iOffset = aiOffsets(1 + iLi); 
            ! Obtain the offset modulo
            iOffset = modulo(iOffset, iWrapL); 
            ! For a normal modulo, the condition iOffset<0 is never met...
            if (iOffset < 0) then
                iOffset = iOffset + iWrapL
            end if

            !fill real input buffer from complex array
            if (mod(iWrapL, 2_i8b) == 1) then
                arBufferIn(1:iWrapL:2) = real(acWrap(iIndWrap:iIndWrap + ((iWrapL + 1)/2 - 1)*aiWrapS(1):aiWrapS(1)), dp); 
            else
                arBufferIn(1:iWrapL:2) = real(acWrap(iIndWrap:iIndWrap + (iWrapL/2 - 1)*aiWrapS(1):aiWrapS(1)), dp); 
            end if
            arBufferIn(2:iWrapL:2) = dimag(acWrap(iIndWrap:iIndWrap + (iWrapL/2 - 1)*aiWrapS(1):aiWrapS(1))); 
            !clear output buffer
            arBufferOut = 0; 
            ! Store the block temporary in a buffer-variable
            if (iOffset + aiBlockL(1) <= iWrapL) then
                arBufferOut(:) = arBufferIn(1 + iOffset:aiBlockL(1) + iOffset)
            else
                iBreak = iWrapL - iOffset; 
                arBufferOut(1:iBreak) = arBufferIn(1 + iOffset:iWrapL)
                arBufferOut(iBreak + 1:aiBlockL(1)) = arBufferIn(1:aiBlockL(1) - iBreak)
            end if

            !store real output buffer in complex array
            if (mod(aiBlockL(1), 2_i8b) == 1) then
                acData(iIndData:iIndData + (aiBlockL(1)/2 - 1)*aiBlockS(1):aiBlockS(1)) = arBufferOut(1:aiBlockL(1) - 1:2) + im*arBufferOut(2:aiBlockL(1) - 1:2); 
                acData(iIndData + ((aiBlockL(1) + 1)/2 - 1)*aiBlockS(1)) = arBufferOut(aiBlockL(1)); 
            else
                acData(iIndData:iIndData + (aiBlockL(1)/2 - 1)*aiBlockS(1):aiBlockS(1)) = arBufferOut(1:aiBlockL(1):2) + im*arBufferOut(2:aiBlockL(1):2); 
            end if

        end do

        ! Stopwatches, stop
        deallocate (arBufferIn, arBufferOut); 
        deallocate (aiLoop); 
    END SUBROUTINE MapVInvRealComplex

    SUBROUTINE MapVt(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVt applies a mapping to the data in a beam to account
!   for the offsets due to the skewness of a beam in the temporal dimension.
!   The grid should be in distribution 1.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset      i8b   1D array containing all the offsets as a function of
!                       the z-position
!   iLi           i8b   loop counter over all xyz locations in the space
!   iIndZ         i8b   beam index in z
!   iWrapTL       i8b   wrap length, in this case equal to 2*iDimT
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::      pcGrid; 
        integer(i8b), allocatable::   aiOffset(:)
        integer(i8b)::            iLi, iWrapTL, iIndZ; 
        CHARACTER(LEN=2048)::       acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetT
!   MapVRealComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVt); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVt", 5)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 1) then
            write (acTemp, '("Error in MapVt. This function requires the grid to be in distribution 1.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! times two because we have a wraparound region in dist. 1
        iWrapTL = 2*cSpace%iDimT

        ! Initialize the table with offsets
        allocate (aiOffset(pcGrid.iD1LocN)); 
        do iLi = 0, pcGrid.iD1LocN - 1
            iIndZ = mod(pcGrid.aiD1Loc(1 + iLi, 3) - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
            aiOffset(1 + iLi) = iBeamOffsetT(cSpace, iIndZ); 
        end do

        ! Call MapVRealComplex with the correct parameters
        call MapVRealComplex(pcGrid%pacD1, pcGrid%pacD1, 2_i8b, (/cSpace.iDimT, pcGrid.iD1LocN/), &
                             (/pcGrid.iD1TS, pcGrid.iD1IS/), &
                             iWrapTL, &
                             (/pcGrid.iD1TS, pcGrid.iD1IS/), &
                             aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVt); 
    END SUBROUTINE MapVt

    SUBROUTINE MapVtInv(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVtInv applies the inverse mapping to the data in a beam
!   in the temporal dimension. The grid should be in distribution 1.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset      i8b   1D array containing all the offsets as a function of
!                       the z-position
!   iLi           i8b   loop counter over all xyz locations in the space
!   iIndZ         i8b   beam index in z
!   iWrapTL       i8b   wrap length, in this case equal to 2*iDimT
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::      pcGrid; 
        integer(i8b), allocatable::   aiOffset(:)
        integer(i8b)::            iLi, iWrapTL, iIndZ
        CHARACTER(LEN=2048)::       acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetT
!   MapVInvRealComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVtInv); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVtInv", 5)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 1) then
            write (acTemp, '("Error in MapVtInv. This function requires the grid to be in distribution 1.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        iWrapTL = 2*cSpace%iDimT

        ! Initialize the table with offsets
        allocate (aiOffset(pcGrid.iD1LocN)); 
        do iLi = 0, pcGrid.iD1LocN - 1
            iIndZ = mod(pcGrid.aiD1Loc(1 + iLi, 3) - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
            aiOffset(1 + iLi) = iBeamOffsetT(cSpace, iIndZ); 
        end do

        ! Call MapVInvReal with the correct parameters
        call MapVInvRealComplex(pcGrid%pacD1, pcGrid%pacD1, 2_i8b, (/2*cSpace.iDimT, cSpace.cGrid.iD1LocN/), &
                                (/cSpace.cGrid.iD1TS, cSpace.cGrid.iD1IS/), &
                                iWrapTL, &
                                (/cSpace.cGrid.iD1TS, cSpace.cGrid.iD1IS/), &
                                aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVtInv); 
    END SUBROUTINE MapVtInv

    SUBROUTINE MapVx(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVx applies a mapping to the data in a beam to account
!   for the offsets due to the skewness of a beam in the x-dimension.
!   The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLy           i8b   loop counter over all y-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLy, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetX
!   MapVComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxy); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVx", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVx. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimY*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLy = 0, cSpace.iDimY - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetX(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimX, cSpace.iDimY, cSpace.iDimZ/), &
                         (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                         cSpace.iDimX, &
                         (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                         aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxy); 
    END SUBROUTINE MapVx

    SUBROUTINE MapVy(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVy applies a mapping to the data in a beam to account
!   for the offsets due to the skewness of a beam in the y-dimension.
!   The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLx           i8b   loop counter over all x-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLx, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetY
!   MapVComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxy); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVy", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVy. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimX*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLx = 0, cSpace.iDimX - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetY(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimY, cSpace.iDimX, cSpace.iDimZ/), &
                         (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                         cSpace.iDimY, &
                         (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                         aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxy); 
    END SUBROUTINE MapVy

    SUBROUTINE MapVx_Small(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVx_Small applies a mapping to the data in a beam to
!   account for the offsets due to the skewness of a beam in the x-dimension.
!   As the wrapping length it takes half the length of the x-dimension, and all
!   other dimensions are taken half as well. This is appropriate for XYZ slices
!   which have a wraparound region but which is not used as such (currently
!   only in InhomContrastOperator_Ali). The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLy           i8b   loop counter over all y-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLy, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetX
!   MapVComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxy); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVx_Small", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVx_Small. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimY/2*cSpace.iDimZ/2)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ/2 - 1
            do iLy = 0, cSpace.iDimY/2 - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ/2) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetX(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimX/2, cSpace.iDimY/2, cSpace.iDimZ/2/), &
                         (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                         cSpace.iDimX/2, &
                         (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                         aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxy); 
    END SUBROUTINE MapVx_Small

    SUBROUTINE MapVy_Small(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVy_Small applies a mapping to the data in a beam to
!   account for the offsets due to the skewness of a beam in the y-dimension.
!   As the wrapping length it takes half the length of the y-dimension, and all
!   other dimensions are taken half as well. This is appropriate for XYZ slices
!   which do have a wraparound region but which is not used as such (currently
!   only in InhomContrastOperator_Ali). The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLx           i8b   loop counter over all y-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLx, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetY
!   MapVComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxy); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVy_Small", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVy_Small. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimX/2*cSpace.iDimZ/2)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ/2 - 1
            do iLx = 0, cSpace.iDimX/2 - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ/2) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetY(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimY/2, cSpace.iDimX/2, cSpace.iDimZ/2/), &
                         (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                         cSpace.iDimY/2, &
                         (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                         aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxy); 
    END SUBROUTINE MapVy_Small

    SUBROUTINE MapVxInv(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVxInv applies the inverse mapping operation of MapVx.
!   The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLy           i8b   loop counter over all y-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLy, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetX
!   MapVInvComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxyInv); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVxInv", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVxInv. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimY*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLy = 0, cSpace.iDimY - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetX(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVInvComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimX, cSpace.iDimY, cSpace.iDimZ/), &
                            (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                            cSpace.iDimX, &
                            (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                            aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxyInv); 
    END SUBROUTINE MapVxInv

    SUBROUTINE MapVyInv(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVyInv applies the inverse mapping operation of MapVy.
!   The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLx           i8b   loop counter over all x-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLx, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetY
!   MapVInvComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxyInv); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVyInv", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVyInv. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimX*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLx = 0, cSpace.iDimX - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetY(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVInvComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimY, cSpace.iDimX, cSpace.iDimZ/), &
                            (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                            cSpace.iDimY, &
                            (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                            aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxyInv); 
    END SUBROUTINE MapVyInv

    SUBROUTINE MapVxInv_Small(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVxInv_Small applies the inverse mapping operation of
!   MapVx_Small. The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLy           i8b   loop counter over all y-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLy, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetX
!   MapVInvComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxyInv); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVxInv_Small", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVxInv_Small. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimY*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLy = 0, cSpace.iDimY - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetX(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVInvComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimX, cSpace.iDimY, cSpace.iDimZ/), &
                            (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                            cSpace.iDimX, &
                            (/pcGrid.iD2XS, pcGrid.iD2YS, pcGrid.iD2ZS/), &
                            aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxyInv); 
    END SUBROUTINE MapVxInv_Small

    SUBROUTINE MapVyInv_Small(cSpace)

! =============================================================================
!
!   Programmer: Jasper de Koning / Koos Huijssen
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
!   The subroutine MapVyInv_Small applies the inverse mapping operation of
!   MapVy_Small. The grid should be in distribution 2.
!
!   WARNING: It is assumed that we have only one frequency, i.e. iDimT=1.
!            This is true for XYZ slices
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space for which the data needs to be mapped
!
        type(Space), target, intent(inout)::      cSpace

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   pcGrid   type(Grid) pointer to the grid in cSpace
!   aiOffset   i8b   1D array containing all the offsets as a function of
!                    the z-position
!   iLx           i8b   loop counter over all x-positions
!   iLz           i8b   loop counter over all z-positions
!   iIndz         i8b   beam index in z
!   iIndex        i8b   index to aiOffset
!   acTemp       char   temporary char array for output log messages
!
        type(Grid), pointer ::         pcGrid; 
        integer(i8b), allocatable ::      aiOffset(:); 
        integer(i8b) ::            iLx, iLz, iIndZ, iIndex
        character(LEN=2048)::          acTemp

! *****************************************************************************
!
!   I/O
!
!   Log file entries
!
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   PrintToLog
!   SwStartAndCount
!   SwStop
!   iBeamOffsetX
!   MapVInvComplex
!
! =============================================================================

        call SwStartAndCount(cswMapVxyInv); call SwStartAndCount(cswMaps); 
        call PrintToLog("Start MapVyInv_Small", 6)

        pcGrid => cSpace.cGrid; 
        ! Check if we are in the correct distribution
        if (pcGrid.iDistr /= 2) then
            write (acTemp, '("Error in MapVyInv_Small. This function requires the grid to be in distribution 2.")')
            call PrintToLog(acTemp, -1)
            write (acTemp, '("However, the given grid was in distribution ", I3, ".")') pcGrid.iDistr
            call PrintToLog(acTemp, -1)
            stop
        end if

        ! Initialization
        allocate (aiOffset(cSpace.iDimX*cSpace.iDimZ)); 
        iIndex = 0; 
        do iLz = 0, cSpace.iDimZ - 1
            do iLx = 0, cSpace.iDimX - 1
                iIndZ = mod(iLz - cSpace%iBeamIndexStartZ, cSpace%iDimZ) + cSpace%iBeamIndexStartZ
                aiOffset(1 + iIndex) = iBeamOffsetY(cSpace, iIndZ); 
                iIndex = iIndex + 1; 
            end do
        end do

        ! Call MapVComplex with the correct parameters
        call MapVInvComplex(pcGrid.pacD2, pcGrid.pacD2, 3_i8b, (/cSpace.iDimY, cSpace.iDimX, cSpace.iDimZ/), &
                            (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                            cSpace.iDimY, &
                            (/pcGrid.iD2YS, pcGrid.iD2XS, pcGrid.iD2ZS/), &
                            aiOffset); 
        deallocate (aiOffset); 
        ! Stopwatches, stop
        call SWStop(cswMaps); call SWStop(cswMapVxyInv); 
    END SUBROUTINE MapVyInv_Small

END MODULE ParnacDataMap
