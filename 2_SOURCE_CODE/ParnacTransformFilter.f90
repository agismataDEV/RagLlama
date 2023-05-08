
MODULE ParnacTransformFilter

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
!   The module ParnacTransformFilter contains subroutines that perform forward
!   and inverse FFT's in all dimensions, and that apply the spatial filters.
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types
    USE Constants
    USE FFTW_cons
    USE ParnacGeneral
    USE ParnacTransformDef
    USE ParnacDataDef
    USE ParnacParamDef
    USE SpecialFun, ONLY: &
        cisi4, cisi, dSinc, ExpInt_Prime
    USE ParnacDataSpaceInit, ONLY: &
        iBeamOffsetT, iBeamOffsetX, iBeamOffsetY
    USE ParnacDataRedist

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   none - if necessary, the different FFTW routines need aliases defined in
!   fftw_wrappers.inc for a correct calling of the routines - this depends
!   on the used FFTW library.
!
!#ifdef FFTW_C_NAMING
    include "fftw_wrappers.inc"
!#endif

    IMPLICIT NONE

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   TransformT            sub   Apply the Forward FFT to the T-dimension
!   TransformTInv         sub   Apply the Inverse FFT to the T-dimension
!   TransformXY           sub   Apply the Forward FFT to the XY-dimensions
!   TransformXYInv        sub   Apply the Inverse FFT to the XY-dimensions
!   TransformXYZ          sub   Apply the Forward FFT to the XYZ-dimensions
!   TransformXYZInv       sub   Apply the Inverse FFT to the XYZ-dimensions
!   TransformT_sml        sub   Apply the Forward FFT to the T-dimension on a
!                               reduced range
!   TransformTInv_sml     sub   Apply the Inverse FFT to the T-dimension on a
!                               reduced range
!   TransformXYZ_sml      sub   Apply the Forward FFT to the XYZ-dimensions on
!                               a reduced range
!   TransformXYZInv_sml   sub   Apply the Inverse FFT to the XYZ-dimensions on
!                               a reduced range
!   FilterSpatial1D       sub   Apply an ideal rectangular filter to a 2D array
!                               in the second dimension - array is real
!   FilterSpatial2D       sub   Apply an ideal spherical filter to a 2D array
!                               in both dimensions - array is real
!   FFilterSpatial1D      sub   Apply an ideal rectangular filter to a 2D array
!                               in the second dimension - array is complex
!   FFilterSpatial2D      sub   Apply an ideal spherical filter to a 2D array
!                               in both dimensions - array is complex
!   FFilterSpaceSpatial3D sub   Apply an ideal spherical filter to a space
!                               in all dimensions
!
! =============================================================================
!
CONTAINS

    SUBROUTINE TransformT(cGrid)

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
!   The subroutine TransformT applies the forward FFT in the T-dimension on a
!   Grid. The data array should be stored in Distr. 1.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::                        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft_r2c
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransT)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 1) then
            write (*, '("Attempted to do a transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, " (should be 1)")') cGrid%iDistr
            write (30, '("Attempted to do a transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, " (should be 1)")') cGrid%iDistr
            stop
        end if

        ! Perform the FFT
!        print *, "Number fft"
!        print *, cGrid%cTransforms%iPlanTransformT
!        print *, "Size in-out"
!        print *,size(cGrid%pacD1(:))
        !write (6,*)'I m here 0'
        call dfftw_execute_dft_r2c(cGrid%cTransforms%iPlanTransformT, cGrid%pacD1(:), cGrid%pacD1(:))
        !write (6,*)'I m here 1'
        !print *, "After dtttw_execute_dft-r2c"

        call SWStop(cswTrans); call SWStop(cswTransT)
    END SUBROUTINE TransformT

    SUBROUTINE TransformTInv(cGrid)

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
!   The subroutine TransformTInv applies the inverse FFT in the T-dimension on a
!   Grid. The data array should be stored in Distr. 1.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::                        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft_c2r
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransTInv)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 1) then
            write (*, '("Attempted to do an inverse transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, "(should be 1)")') cGrid%iProcID
            write (30, '("Attempted to do an inverse transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, "(should be 1)")') cGrid%iProcID
            stop
        end if

        ! Perform the FFT
!        print *, "Number ifft"
!        print *, cGrid%cTransforms%iPlanTransformT_inv
!        print *, "Size in-out"
!        print *,size(cGrid%pacD1(:))

        call dfftw_execute_dft_c2r(cGrid%cTransforms%iPlanTransformT_inv, cGrid%pacD1(:), cGrid%pacD1(:))

        call SWStop(cswTrans); call SWStop(cswTransTInv)
    END SUBROUTINE TransformTInv

    SUBROUTINE TransformXY(cGrid)

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
!   The subroutine TransformXY applies the forward FFT in the X- and Y-
!   dimensions on a Grid. The data array should be stored in Distr. 2. It is
!   assumed that iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::                        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXY)
        call PrintToLog("TransformXY", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XY-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XY-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT
        call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXY, cGrid%pacD2, cGrid%pacD2); 
        call SWStop(cswTrans); call SWStop(cswTransXY)
    END SUBROUTINE TransformXY

    SUBROUTINE TransformXYInv(cGrid)

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
!   The subroutine TransformXY applies the inverse FFT in the X- and Y-
!   dimensions on a Grid. The data array should be stored in Distr. 2. It is
!   assumed that iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::                        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! ============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXYInv)
        call PrintToLog("TransformXYInv", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XY-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XY-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT
        call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXY_inv, cGrid%pacD2, cGrid%pacD2); 
        call SWStop(cswTrans); call SWStop(cswTransXYInv)
    END SUBROUTINE TransformXYInv

    SUBROUTINE TransformXYZ(cGrid, bFullZ)

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
!   The subroutine TransformXYZ applies the forward FFT in the X-, Y- and Z-
!   dimensions on a Grid. The data array should be stored in Distr. 2. It is
!   assumed that iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!   bFullZ  i      lgt        Whether or not to apply the FFT's in X and Y over
!                             all Z-positions, or only over the relevant ones.
!                             The FFT in Z is always applied over all
!                             Z-positions.
!
        type(Grid), intent(inout) ::        cGrid; 
        logical, intent(in)::                        bFullZ; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXYZ)
        call PrintToLog("TransformXYZ", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT, if bFullZ==.false. then only the first half z is XY-transformed
        if (bFullZ) then
!            print *, "Sono qui"
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXYZ, cGrid%pacD2(:), cGrid%pacD2(:)); 
        else
!            print *, "Sono qui 2"
!            print *, size(cGrid%pacD2(:))
!            print *, "Number"
!            print *, cGrid%cTransforms%iPlanTransformXY_relevantZ
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXY_relevantZ, cGrid%pacD2(:), cGrid%pacD2(:)); 
!                print *, "Sono qui 3"
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformZ, cGrid%pacD2(:), cGrid%pacD2(:)); 
        end if

        call SWStop(cswTrans); call SWStop(cswTransXYZ)
    END SUBROUTINE TransformXYZ

    SUBROUTINE TransformXYZInv(cGrid, bFullZ)

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
!   The subroutine TransformXYZ applies the inverse FFT in the X-, Y- and Z-
!   dimensions on a Grid. The data array should be stored in Distr. 2. It is
!   assumed that iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!   bFullZ  i      lgt        Whether or not to apply the FFT's in X and Y over
!                             all Z-positions, or only over the relevant ones.
!                             The FFT in Z is always applied over all
!                             Z-positions.
!
        type(Grid), intent(inout) ::        cGrid; 
        logical, intent(in)::                        bFullZ; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXYZInv)
        call PrintToLog("TransformXYZInv", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT, if bFullZ==.false. then only the first half z is XY-inv-transformed
        if (bFullZ) then
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXYZ_inv, cGrid%pacD2(:), cGrid%pacD2(:)); 
        else
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformZ_inv, cGrid%pacD2(:), cGrid%pacD2(:)); 
            call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXY_relevantZ_inv, cGrid%pacD2(:), cGrid%pacD2(:)); 
        end if
        call SWStop(cswTrans); call SWStop(cswTransXYZInv)
    END SUBROUTINE TransformXYZInv

    SUBROUTINE TransformT_sml(cGrid)

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
!   The subroutine TransformT_sml applies the forward FFT in the T-dimension
!   on a Grid. The data array should be stored in Distr. 1. As the length of
!   the T-axis, we use iD0TL, that is the time axis without the wraparound
!   region.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft_r2c
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransT)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 1) then
            write (*, '("Attempted to do a transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, " (should be 1)")') cGrid%iDistr
            write (30, '("Attempted to do a transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, " (should be 1)")') cGrid%iDistr
            stop
        end if

        ! Perform the FFT
!        print *, "Number fft_small"
!        print *, cGrid%cTransforms%iPlanTransformT_sml
!        print *, "Size in-out"
!        print *,size(cGrid%pacD1(:))
        call dfftw_execute_dft_r2c(cGrid%cTransforms%iPlanTransformT_sml, cGrid%pacD1(:), cGrid%pacD1(:))

        call SWStop(cswTrans); call SWStop(cswTransT)
    END SUBROUTINE TransformT_sml

    SUBROUTINE TransformTInv_sml(cGrid)

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
!   The subroutine TransformTInv_sml applies the inverse FFT in the T-dimension
!   on a Grid. The data array should be stored in Distr. 1. As the length of
!   the T-axis, we use iD0TL, that is the time axis without the wraparound
!   region.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft_c2r
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransTInv)
        call PrintToLog("AttenuationEval6", 2);  
        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 1) then
            write (*, '("Attempted to do an inverse transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, "(should be 1)")') cGrid%iProcID
            write (30, '("Attempted to do an inverse transform on a grid, with respect to the T-axis, however the grid is in distribution ", I3, "(should be 1)")') cGrid%iProcID
            stop
        end if

        ! Perform the FFT
!        print *, "Number ifft_small"
!        print *, cGrid%cTransforms%iPlanTransformT_sml_inv
!        print *, "Size in-out"
!        print *,size(cGrid%pacD1(:))
        call dfftw_execute_dft_c2r(cGrid%cTransforms%iPlanTransformT_sml_inv, cGrid%pacD1(:), cGrid%pacD1(:))

        call SWStop(cswTrans); call SWStop(cswTransTInv)
    END SUBROUTINE TransformTInv_sml

    SUBROUTINE TransformXYZ_sml(cGrid)

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
!   The subroutine TransformXYZ_sml applies the forward FFT in the X-, Y- and
!   Z-dimensions on a Grid. The data array should be stored in Distr. 2. As the
!   length of the X-, Y- and Z- axes, we use iD2XL/2, iD2YL/2 and iD2ZL/2, that
!   is the XYZ-axes without the wraparound regions. It is assumed that
!   iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXYZ)
        call PrintToLog("TransformXYZ_sml", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT
        call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXYZ_sml, cGrid%pacD2(:), cGrid%pacD2(:)); 
        call SWStop(cswTrans); call SWStop(cswTransXYZ)
    END SUBROUTINE TransformXYZ_sml

    SUBROUTINE TransformXYZInv_sml(cGrid)

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
!   The subroutine TransformXYZInv_sml applies the inverse FFT in the X-, Y-
!   and Z-dimensions on a Grid. The data array should be stored in Distr. 2. As
!   the length of the X-, Y- and Z- axes, we use iD2XL/2, iD2YL/2 and iD2ZL/2,
!   that is the XYZ-axes without the wraparound regions. It is assumed that
!   iD2LocN=1, i.e. the Grid belongs to a convolution slice space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGrid   io   type(Grid)   Grid containing the data array to which the
!                             FFT is applied
!
        type(Grid), intent(inout) ::        cGrid; 
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
!   SwStartAndCount
!   SwStop
!   dfftw_execute_dft
!
! =============================================================================

        call SwStartAndCount(cswTrans); call SwStartAndCount(cswTransXYZInv)
        call PrintToLog("TransformXYZInv_sml", 5)

        ! Standard test, to make sure that we are in the correct distribution
        if (cGrid%iDistr /= 2) then
            write (*, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            write (30, '("Attempted to do a transform on a grid, with respect to the XYZ-axis, however the grid is in distribution ", I3, " (should be 2)")') cGrid%iDistr, cGrid%iProcID
            stop
        end if

        ! Perform the FFT, if bFullZ==.false. then only the first half z is XY-inv-transformed
        call dfftw_execute_dft(cGrid%cTransforms%iPlanTransformXYZ_sml_inv, cGrid%pacD2(:), cGrid%pacD2(:)); 
        call SWStop(cswTrans); call SWStop(cswTransXYZInv)
    END SUBROUTINE TransformXYZInv_sml

    SUBROUTINE FilterSpatial1D(shapein, shapeout, dxin, factor, Kc)

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
!   The subroutine FilterSpatial1D applies an ideal rectangular filter with
!   angular cutoff frequency Kc to the 2nd dimension of the real, 2D input
!   array. It also reduces the sampling in that dimension by the factor given
!   as argument.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   shapein   i   dp   2D input array
!   shapeout  o   dp   2D output array; the second dimension has been filtered
!                      and has reduced sampling
!   dxin      i   dp   step size in the second dimension
!   factor    i   i8b  oversampling factor in shapein, or reduction factor in
!                      the sampling of shapeout
!   Kc        i   dp   angular cutoff frequency of the filter
!
        real(dp), intent(IN) :: shapein(:, :)
        real(dp), intent(OUT) :: shapeout(:, :)
        real(dp), intent(IN) :: dxin
        integer(i8b), intent(IN) :: factor
        real(dp), intent(IN) :: Kc

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i         i8b   loop counter
!   j         i8b   loop counter
!   Kfin      i8b   number of points in the first dimension in the input array
!   Kfout     i8b   number of points in the first dimension in the output array
!   Mfin      i8b   number of points in the second dimension in the input array
!   Mfout     i8b   number of points in the second dimension in the output array
!   error     i8b   error number
!   filter    dp    1D array with the filter in the spectral domain
!   shapeoutlarge  dp  2D large filtered output array (resampling has not yet
!                      been applied)
!   fshapein  dpc   input array transformed in the second dimension
!   k2        dp    k-axis, i.e. axis containing all spatial frequencies
!   dk2       dp    step size in the k-axis
!   kc2pi     dp    filter cutoff frequency
!   plan_1D_transform   i8b   FFTW transform plan
!
        integer(i8b) :: i, j, Kfin, Mfin, Kfout, Mfout, error
        real(dp), allocatable :: filter(:), shapeoutlarge(:, :)
        complex(dpc), allocatable :: fshapein(:, :)
        real(dp), allocatable :: k2(:)
        real(dp) :: dk2, kc2pi
!        real(dp) :: kmaxfinal,band
        integer(i8b) :: plan_1D_transform

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
!   alloctest
!   dfftw_plan_many_dft_r2c
!   dfftw_execute
!   dfftw_plan_many_dft_c2r
!   dfftw_destroy_plan
!
! =============================================================================

        Kfin = size(shapein, 1)
        Mfin = size(shapein, 2)

        Kfout = size(shapeout, 1)
        Mfout = size(shapeout, 2)

        if ((Kfin /= Kfout) .OR. (Mfin /= Mfout*factor)) then
            write (*, "('Error in FilterSpatial1D in init_source')")
            stop
        end if

        allocate (filter(1:(Mfin/2 + 1)), fshapein(1:Kfin, 1:(Mfin/2 + 1)), &
                  shapeoutlarge(1:Kfin, 1:Mfin), k2(1:(Mfin/2 + 1)), STAT=error)
        call alloctest(error, 'filter,fshape.. in FilterSpatial1D')

        !Define spatial frequency k2
        dk2 = 1/(Mfin*dxin)

        k2 = real((/(i, i=0, Mfin/2)/), dp)*dk2

        !use smooth filter between Knyq and Kc
!        filter=0.0_dp
!        kc2pi=Kc/(two_pi)
!        kmaxfinal=1.0_dp/(2.0_dp*dxin*real(factor,dp))
!        band=kc2pi-kmaxfinal
        ! kmaxfinal is the maximum k in the final, normally sampled shape;
        ! band gives the fall band between kmaxfinal and kc, in which we have a
        ! smooth decay of the filter; it is assumed that Kc is given in wavenumber K,
        ! not in spectral parameter k = K / 2 pi (same as f = omega / 2 pi)
        ! for 2D, we need sqrt(2)*kmaxfinal instead of kmaxfinal
        ! for 3D, we need sqrt(3)*kmaxfinal
!        do j=1,Mfin/2
!                if (k2(j)<kc2pi-band) then
!                        filter(j)=1.0_dp
!                else if (k2(j)<kc2pi) then
!                        filter(j)=(0.5_dp+9.0_dp/16.0_dp*cos(pi*((k2(j)-(kc2pi-band))/band))                &
!                                                         -1.0_dp/16.0_dp*cos(3.0_dp*pi*((k2(j)-(kc2pi-band))/band)))
!                end if
!        end do

        !use sharp filter if Knyq=Kc
        filter = 0.0_dp
        kc2pi = Kc/(two_pi)
        do j = 1, Mfin/2
            if (k2(j) < kc2pi + 1e-6) then
                filter(j) = 1.0_dp
            end if
        end do

        !transform fshape and filter it
        call dfftw_plan_many_dft_r2c(plan_1D_transform, 1, int(Mfin), int(Kfin), &
                                     shapein, 0, int(Kfin), 1, &
                                     fshapein, 0, int(Kfin), 1, FFTW_ESTIMATE)
        call dfftw_execute(plan_1D_transform)

        do j = 1, Kfin
            fshapein(j, :) = fshapein(j, :)*filter
        end do

        call dfftw_plan_many_dft_c2r(plan_1D_transform, 1, int(Mfin), int(Kfin), &
                                     fshapein, 0, int(Kfin), 1, &
                                     shapeoutlarge, 0, int(Kfin), 1, FFTW_ESTIMATE)
        call dfftw_execute(plan_1D_transform)

        call dfftw_destroy_plan(plan_1D_transform)

        !renormalize output by 1/(Mfin)
        shapeoutlarge = shapeoutlarge/(Mfin)

        !reduce number of samples by striding with step factor through shapein
        shapeout = shapeoutlarge(:, ::factor)

        deallocate (filter, fshapein, shapeoutlarge, k2)

    END SUBROUTINE FilterSpatial1D

    SUBROUTINE FilterSpatial2D(shapein, shapeout, dxin, factor, Kc)

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
!   The subroutine FilterSpatial2D applies an ideal spherical filter with
!   angular cutoff frequency Kc to both dimensions of the real, 2D input
!   array. It also reduces the sampling in both dimensions by the factor given
!   as argument.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   shapein   i   dp   2D input array
!   shapeout  o   dp   2D output array; both dimensions have been filtered and
!                      have reduced sampling
!   dxin      i   dp   step size in both dimensions
!   factor    i   i8b  oversampling factor in shapein, or reduction factor in
!                      the sampling of shapeout
!   Kc        i   dp   angular cutoff frequency of the filter
!
        real(dp), intent(IN) :: shapein(:, :)
        real(dp), intent(OUT) :: shapeout(:, :)
        real(dp), intent(IN) :: dxin
        integer(i8b), intent(IN) :: factor
        real(dp), intent(IN) :: Kc

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i         i8b   loop counter
!   j         i8b   loop counter
!   Lfin      i8b   number of points in the first dimension in the input array
!   Lfout     i8b   number of points in the first dimension in the output array
!   Mfin      i8b   number of points in the second dimension in the input array
!   Mfout     i8b   number of points in the second dimension in the output array
!   error     i8b   error number
!   filter    dp    2D array with the filter in the spectral domain
!   shapeoutlarge  dp  2D large filtered output array (resampling has not yet
!                      been applied)
!   fshapein  dpc   input array transformed in both dimensions
!   k1        dp    k1-axis, i.e. axis containing all spatial frequencies of the
!                   first dimension
!   k2        dp    k2-axis, i.e. axis containing all spatial frequencies of the
!                   second dimension
!   dk1       dp    step size in k1
!   dk2       dp    step size in k2
!   kc2pi     dp    filter cutoff frequency
!   k1k2      dp    sqrt(k1**2+k2**2)
!   plan_2D_transform   i8b   FFTW transform plan
!
        integer(i8b) :: i, j, Lfin, Mfin, Lfout, Mfout, error
        real(dp), allocatable :: filter(:, :), shapeoutlarge(:, :)
        complex(dpc), allocatable :: fshapein(:, :)
        real(dp), allocatable :: k1(:), k2(:)
        real(dp) :: K1max, K2max, dk1, dk2, kc2pi, k1k2
!        real(dp) :: kmaxfinal,band
        integer(i8b) :: plan_2D_transform

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
!   alloctest
!   dfftw_plan_many_dft_r2c_2d
!   dfftw_execute
!   dfftw_plan_many_dft_c2r_2d
!   dfftw_destroy_plan
!
! =============================================================================

        call PrintToLog("FilterSpatial2D", 3); 
        Lfin = size(shapein, 1)
        Mfin = size(shapein, 2)

        Lfout = size(shapeout, 1)
        Mfout = size(shapeout, 2)

        if ((Lfin /= Lfout*factor) .OR. (Mfin /= Mfout*factor)) then
            write (*, "('Error in FilterSpatial2D in init_source')")
            stop
        end if

        allocate (filter(1:(Lfin/2 + 1), 1:Mfin), fshapein(1:(Lfin/2 + 1), 1:Mfin), &
                  shapeoutlarge(1:Lfin, 1:Mfin), k1(1:(Lfin/2 + 1)), k2(1:Mfin), STAT=error)

        call alloctest(error, 'filter,fshape.. in FilterSpatial2D')

        !Define spatial frequency k1 and k2
        dk1 = 1/(Lfin*dxin)
        dk2 = 1/(Mfin*dxin)

        k1 = real((/(i, i=0, Lfin/2)/), dp)*dk1
        k2 = real((/(i, i=0, (Mfin + 1)/2 - 1), (i, i=-Mfin/2, -1)/), dp)*dk2

        !use smooth filter between Knyq and Kc
!        filter=0.0_dp
!        kc2pi=Kc/(two_pi)
!        kmaxfinal=1.0_dp/(2.0_dp*dxin*real(factor,dp))
!        band=kc2pi-sqrt(2.0_dp)*kmaxfinal                !
        ! kmaxfinal is the maximum k in the final, normally sampled shape;
        ! band gives the fall band between kmaxfinal and kc, in which we have a
        ! smooth decay of the filter; it is assumed that Kc is given in wavenumber K,
        ! not in spectral parameter k = K / 2 pi (same as f = omega / 2 pi)
        ! for 2D, we need sqrt(2)*kmaxfinal instead of kmaxfinal
        ! for 3D, we need sqrt(3)*kmaxfinal
!        do i=1,Lfin/2
!                do j=1,Mfin
!                        k1k2=sqrt(k1(i)**2+k2(j)**2)
!                        if (k1k2<kc2pi-band) then
!                                filter(i,j)=1.0_dp
!                        else if (k1k2<kc2pi) then
!                                filter(i,j)=(0.5_dp+9.0_dp/16.0_dp*cos(pi*((k1k2-(kc2pi-band))/band))                &
!                                                                   -1.0_dp/16.0_dp*cos(3.0_dp*pi*((k1k2-(kc2pi-band))/band)))
!                        end if
!                end do
!        end do

        !use sharp filter if Knyq=Kc
        filter = 0.0_dp
        kc2pi = Kc/(two_pi)
        do i = 1, Lfin/2
            do j = 1, Mfin/2
                k1k2 = sqrt(k1(i)**2 + k2(j)**2)
                if (k1k2 < kc2pi + 1e-6) then
                    filter(i, j) = 1.0_dp
                    if (j > 1) then
                        filter(i, Mfin + 2 - j) = 1.0_dp
                    end if
                end if
            end do
        end do

        !transform fshape and filter it
        call dfftw_plan_dft_r2c_2d(plan_2D_transform, int(Lfin), int(Mfin), shapein, fshapein, FFTW_ESTIMATE)
        call dfftw_execute(plan_2D_transform)
        call dfftw_destroy_plan(plan_2D_transform)

        fshapein = fshapein*filter
        fshapein = fshapein*filter

        !reduce number of samples by striding with step factor through shapein
        call dfftw_plan_dft_c2r_2d(plan_2D_transform, int(Lfin), int(Mfin), fshapein, shapeoutlarge, FFTW_ESTIMATE)
        call dfftw_execute(plan_2D_transform)
        call dfftw_destroy_plan(plan_2D_transform)
        shapeout = shapeoutlarge(::factor, ::factor)

!        !alternative procedure to reduce number of samples: by selecting low frequency terms only
!        fshapeout(:,1:(Mfout+1)/2)=fshapein(1:(Lfout/2+1),1:(Mfout+1)/2)
!        fshapeout(:,Mfout-Mfout/2+1:Mfout)=fshapein(1:(Lfout/2+1),Mfin-Mfout/2+1:Mfin)
!        call dfftw_plan_dft_c2r_2d(plan_2D_transform,int(Lfout),int(Mfout),fshapeout,shapeout,FFTW_ESTIMATE)
!        call dfftw_execute(plan_2D_transform)
!        call dfftw_destroy_plan(plan_2D_transform)

        !renormalize output by 1/(Lfin*Mfin)
        shapeout = shapeout/(Lfin*Mfin)

        deallocate (filter, fshapein, shapeoutlarge, k1, k2)

        call PrintToLog("End FilterSpatial2D", 5); 
    END SUBROUTINE FilterSpatial2D

    SUBROUTINE FFilterSpatial1D(fshapein, fshapeout, dxin, factor, Kc)

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
!   The subroutine FFilterSpatial1D applies an ideal rectangular filter with
!   angular cutoff frequency Kc to the 2nd dimension of the complex, 2D input
!   array. It also reduces the sampling in that dimension by the factor given
!   as argument.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   fshapein   i   dpc  2D input array
!   fshapeout  o   dpc  2D output array; the second dimension has been filtered
!                       and has reduced sampling
!   dxin       i   dp   step size in the second dimension
!   factor     i   i8b  oversampling factor in shapein, or reduction factor in
!                       the sampling of shapeout
!   Kc         i   dp   angular cutoff frequency of the filter
!
        complex(dpc), intent(IN) :: fshapein(:, :)
        complex(dpc), intent(OUT) :: fshapeout(:, :)
        real(dp), intent(IN) :: dxin
        integer(i8b), intent(IN) :: factor
        real(dp), intent(IN) :: Kc

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i         i8b   loop counter
!   j         i8b   loop counter
!   Kfin      i8b   number of points in the first dimension in the input array
!   Kfout     i8b   number of points in the first dimension in the output array
!   Mfin      i8b   number of points in the second dimension in the input array
!   Mfout     i8b   number of points in the second dimension in the output array
!   error     i8b   error number
!   filter    dp    1D array with the filter in the spectral domain
!   fshapeoutlarge dpc  2D large filtered output array (resampling has not yet
!                       been applied)
!   ffshapein dpc   input array transformed in the second dimension
!   k2        dp    k-axis, i.e. axis containing all spatial frequencies
!   dk2       dp    step size in the k-axis
!   kc2pi     dp    filter cutoff frequency
!   plan_1D_transform   i8b   FFTW transform plan
!
        integer(i8b) :: i, j, Kfin, Mfin, Kfout, Mfout, error
        real(dp), allocatable :: filter(:)
        complex(dpc), allocatable :: ffshapein(:, :), fshapeoutlarge(:, :)
        real(dp), allocatable :: k2(:)
        real(dp) :: dk2, kc2pi
!        real(dp) :: kmaxfinal,band
        integer(i8b) :: plan_1D_transform

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
!   alloctest
!   dfftw_plan_many_dft
!   dfftw_execute
!   dfftw_plan_many_dft
!   dfftw_destroy_plan
!
! =============================================================================

        Kfin = size(fshapein, 1)
        Mfin = size(fshapein, 2)

        Kfout = size(fshapeout, 1)
        Mfout = size(fshapeout, 2)

        if ((Kfin /= Kfout) .OR. (Mfin /= Mfout*factor)) then
            write (*, "('Error in FFilterSpatial1D in init_source')")
            stop
        end if

        allocate (filter(1:(Mfin)), ffshapein(1:Kfin, 1:(Mfin)), &
                  fshapeoutlarge(1:Kfin, 1:Mfin), k2(1:(Mfin)), STAT=error)
        call alloctest(error, 'filter,ffshape.. in FFilterSpatial2D')

        !Define spatial frequency k2
        dk2 = 1/(Mfin*dxin)
        k2 = real((/(i, i=0, (Mfin - 1)/2), (i, i=-Mfin/2, -1)/), dp)*dk2

        !use sharp filter if Knyq=Kc
        filter = 0.0_dp
        kc2pi = Kc/two_pi
        do j = 1, Mfin
            if (abs(k2(j)) < kc2pi + 1e-6) then
                filter(j) = 1.0_dp
            end if
        end do

        !transform fshape and filter it
        call dfftw_plan_many_dft(plan_1D_transform, 1, int(Mfin), int(Kfin), &
                                 fshapein, 0, int(Kfin), 1, &
                                 ffshapein, 0, int(Kfin), 1, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan_1D_transform)

        do j = 1, Kfin
            ffshapein(j, :) = ffshapein(j, :)*filter
        end do

        call dfftw_plan_many_dft(plan_1D_transform, 1, int(Mfin), int(Kfin), &
                                 ffshapein, 0, int(Kfin), 1, &
                                 fshapeoutlarge, 0, int(Kfin), 1, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan_1D_transform)

        call dfftw_destroy_plan(plan_1D_transform)

        !reduce number of samples by striding with step factor through shapein
        fshapeout = fshapeoutlarge(:, ::factor)

        !renormalize output by 1/(Mfin)
        fshapeout = fshapeout/(Mfin)

        deallocate (filter, ffshapein, fshapeoutlarge, k2)

    END SUBROUTINE FFilterSpatial1D

    SUBROUTINE FFilterSpatial2D(fshapein, fshapeout, dxin, factor, Kc)

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
!   The subroutine FFilterSpatial2D applies an ideal spherical filter with
!   angular cutoff frequency Kc to both dimensions of the complex, 2D input
!   array. It also reduces the sampling in both dimensions by the factor given
!   as argument.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   fshapein   i   dpc  2D input array
!   fshapeout  o   dpc  2D output array; both dimensions have been filtered and
!                       have reduced sampling
!   dxin       i   dp   step size in both dimensions
!   factor     i   i8b  oversampling factor in shapein, or reduction factor in
!                       the sampling of shapeout
!   Kc         i   dp   angular cutoff frequency of the filter
!
        complex(dpc), intent(IN) :: fshapein(:, :)
        complex(dpc), intent(OUT) :: fshapeout(:, :)
        real(dp), intent(IN) :: dxin
        integer(i8b), intent(IN) :: factor
        real(dp), intent(IN) :: Kc

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i         i8b   loop counter
!   j         i8b   loop counter
!   Lfin      i8b   number of points in the first dimension in the input array
!   Lfout     i8b   number of points in the first dimension in the output array
!   Mfin      i8b   number of points in the second dimension in the input array
!   Mfout     i8b   number of points in the second dimension in the output array
!   error     i8b   error number
!   filter    dp    2D array with the filter in the spectral domain
!   fshapeoutlarge dpc  2D large filtered output array (resampling has not yet
!                       been applied)
!   ffshapein dpc   input array transformed in both dimensions
!   k1        dp    k1-axis, i.e. axis containing all spatial frequencies of the
!                   first dimension
!   k2        dp    k2-axis, i.e. axis containing all spatial frequencies of the
!                   second dimension
!   dk1       dp    step size in k1
!   dk2       dp    step size in k2
!   kc2pi     dp    filter cutoff frequency
!   k1k2      dp    sqrt(k1**2+k2**2)
!   plan_2D_transform   i8b   FFTW transform plan
!
        integer(i8b) :: i, j, Lfin, Mfin, Lfout, Mfout, error
        real(dp), allocatable :: filter(:, :)
        complex(dpc), allocatable :: ffshapein(:, :), fshapeoutlarge(:, :)
        real(dp), allocatable :: k1(:), k2(:)
        real(dp) :: K1max, K2max, dk1, dk2, kc2pi, k1k2
!        real(dp) :: kmaxfinal,band
        integer(i8b) :: plan_2D_transform

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
!   alloctest
!   dfftw_plan_many_dft_2d
!   dfftw_execute
!   dfftw_plan_many_dft_2d
!   dfftw_destroy_plan
!
! =============================================================================

        Lfin = size(fshapein, 1)
        Mfin = size(fshapein, 2)

        Lfout = size(fshapeout, 1)
        Mfout = size(fshapeout, 2)

        if ((Lfin /= Lfout*factor) .OR. (Mfin /= Mfout*factor)) then
            write (*, "('Error in FFilterSpatial1D in init_source')")
            stop
        end if

        allocate (filter(1:Lfin, 1:Mfin), ffshapein(1:Lfin, 1:Mfin), &
                  fshapeoutlarge(1:Lfin, 1:Mfin), k1(1:Lfin), k2(1:Mfin), STAT=error)
        call alloctest(error, 'filter,fshape.. in FFilterSpatial1D')

        !Define spatial frequency k1 and k2
        K1max = 1/(2*dxin)
        K2max = 1/(2*dxin)
        dk1 = 1/(Lfin*dxin)
        dk2 = 1/(Mfin*dxin)
        k1 = real((/(i, i=0, (Lfin - 1)/2), (i, i=-Lfin/2, -1)/), dp)*dk1
        k2 = real((/(i, i=0, (Mfin - 1)/2), (i, i=-Mfin/2, -1)/), dp)*dk2

        !use smooth filter between Knyq and Kc
!        filter=0.0_dp
!        kc2pi=Kc/(two_pi)
!        kmaxfinal=1.0_dp/(2.0_dp*dxin*real(factor,dp))
!        band=kc2pi-sqrt(2.0_dp)*kmaxfinal                !
        ! kmaxfinal is the maximum k in the final, normally sampled shape;
        ! band gives the fall band between kmaxfinal and kc, in which we have a
        ! smooth decay of the filter; it is assumed that Kc is given in wavenumber K,
        ! not in spectral parameter k = K / 2 pi (same as f = omega / 2 pi)
        ! for 2D, we need sqrt(2)*kmaxfinal instead of kmaxfinal
        ! for 3D, we need sqrt(3)*kmaxfinal
!        do i=1,Lfin
!                do j=1,Mfin
!                        k1k2=sqrt(k1(i)**2+k2(j)**2)
!                        if (k1k2<kc2pi-band) then
!                                filter(i,j)=1.0_dp
!                        else if (k1k2<kc2pi) then
!                                filter(i,j)=(0.5_dp+9.0_dp/16.0_dp*cos(pi*((k1k2-(kc2pi-band))/band))                &
!                                                                   -1.0_dp/16.0_dp*cos(3.0_dp*pi*((k1k2-(kc2pi-band))/band)))
!                        end if
!                end do
!        end do

        !use sharp filter if Knyq=Kc
        filter = 0.0_dp
        kc2pi = Kc/(two_pi)
        do i = 1, Lfin
            do j = 1, Mfin
                k1k2 = sqrt(k1(i)**2 + k2(j)**2)
                if (k1k2 < kc2pi + 1e-6) then
                    filter(i, j) = 1.0_dp
                end if
            end do
        end do

        !transform fshape and filter it
        call dfftw_plan_dft_2d(plan_2D_transform, int(Lfin), int(Mfin), fshapein, ffshapein, FFTW_FORWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan_2D_transform)

        ffshapein = ffshapein*filter

        call dfftw_plan_dft_2d(plan_2D_transform, int(Lfin), int(Mfin), ffshapein, fshapeoutlarge, FFTW_BACKWARD, FFTW_ESTIMATE)
        call dfftw_execute(plan_2D_transform)

        call dfftw_destroy_plan(plan_2D_transform)

        !reduce number of samples by striding with step factor through shapein
        fshapeout = fshapeoutlarge(::factor, ::factor)

!        !alternative procedure to reduce number of samples: by selecting low frequency terms only
!        fshapeout(:,1:(Mfout+1)/2)=fshapein(1:(Lfout/2+1),1:(Mfout+1)/2)
!        fshapeout(:,Mfout-Mfout/2+1:Mfout)=fshapein(1:(Lfout/2+1),Mfin-Mfout/2+1:Mfin)
!!! and of course the parts Lfin-Lfout/2+1:Lfin etc...
!!! is not implemented...
!
!        call dfftw_plan_dft_2d(plan_2D_transform,int(Lfout),int(Mfout),fshapeout,shapeout,FFTW_BACKWARD,FFTW_ESTIMATE)
!        call dfftw_execute(plan_2D_transform)
!
!        call dfftw_destroy_plan(plan_2D_transform)

        !renormalize output by 1/(Lfin*Mfin)
        fshapeout = fshapeout/(Lfin*Mfin)

        deallocate (filter, ffshapein, fshapeoutlarge, k1, k2)

    END SUBROUTINE FFilterSpatial2D

    SUBROUTINE FFilterSpaceSpatial3D(cSpace, iOversampling)

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
!   The subroutine FFilterSpaceSpatial3D applies an ideal spherical filter to
!   the three spatial dimensions of the Grid in the space. The grid is assumed
!   to be in Distr. 2 and to be oversampled by a factor iOversampling. The
!   cutoff frequency of the spherical filter is the Nyquist frequency of the
!   original sampling (i.e. without the oversampling factor). The output is
!   a filtered but still oversampled space.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space containing the grid to which the filter
!                               is applied
!   iOversampling   i   i8b     Oversampling factor
!
        type(Space), intent(inout), target :: cSpace
        integer(i8b), intent(in) :: iOversampling

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   error   i8b   error number
!   dK1     dp    K1-axis containing all frequencies of the 1st dimension
!   dK2     dp    K2-axis containing all frequencies of the 2nd dimension
!   dK3     dp    K3-axis containing all frequencies of the 3rd dimension
!   dDk1    dp    step size in the K1-axis
!   dDk2    dp    step size in the K2-axis
!   dDk3    dp    step size in the K3-axis
!   dk1k2k3   dp  sqrt(k1**2+k2**2+k3**2)
!   dKcutoff  dp  cutoff frequency of the filter
!   dHalfbandstart dp  start of the half band filter
!   dHalfbandend   dp  end of the half band filter
!   iIndex   i8b   index into the data array
!   iIndex1  i8b   loop counter in the X-dimension
!   iIndex2  i8b   loop counter in the Y-dimension
!   iIndex3  i8b   loop counter in the Z-dimension
!   iFilterselect i8b  whether to use an ideal (1) or half band (2) filter
!   i        i8b   temporary variable
!   pcGrid   type(Grid)   pointer to the Grid inside the space
!
        integer(i8b) :: error
        real(dp), allocatable :: dK1(:), dK2(:), dK3(:)
        real(dp) :: dDk1, dDk2, dDk3, dk1k2k3, dKcutoff, dHalfbandstart, dHalfbandend
        integer(i8b) :: iIndex, iIndex1, iIndex2, iIndex3, iFilterselect, i
        type(Grid), pointer :: pcGrid

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
!   alloctest
!   transformXYZ
!   transformXYZInv
!
! =============================================================================

        pcGrid => cSpace%cGrid

        allocate (dK1(cSpace%iDimX), dK2(cSpace%iDimY), dK3(cSpace%iDimZ), STAT=error)
        call alloctest(error, 'FFilterSpaceSpatial3D')

        iFilterselect = 1
        ! 1 = ideal filter
        ! 2 = halfband filter

        !Define spatial frequency axes
        dKcutoff = two_pi*cSpace%dFnyq/iOversampling
        dDk1 = two_pi/(cSpace%iDimX*cSpace%dDx)
        dDk2 = two_pi/(cSpace%iDimY*cSpace%dDx)
        dDk3 = two_pi/(cSpace%iDimZ*cSpace%dDx)
        dK1 = real((/(i, i=0, (cSpace%iDimX - 1)/2), (i, i=-cSpace%iDimX/2, -1)/), dp)*dDk1
        dK2 = real((/(i, i=0, (cSpace%iDimY - 1)/2), (i, i=-cSpace%iDimY/2, -1)/), dp)*dDk2
        dK3 = real((/(i, i=0, (cSpace%iDimZ - 1)/2), (i, i=-cSpace%iDimZ/2, -1)/), dp)*dDk3

        !transform cSpace and filter it
!        print *, "Punto di Interesse Odierno 5"
        call TransformXYZ(pcGrid, .true.)

        if (iFilterselect == 1) then
            !Use Ideal filter
            do iIndex1 = 0, cSpace%iDimX - 1
                do iIndex2 = 0, cSpace%iDimY - 1
                    do iIndex3 = 0, cSpace%iDimZ - 1
                        dk1k2k3 = sqrt(dK1(1 + iIndex1)**2 + dK2(1 + iIndex2)**2 + dK3(1 + iIndex3)**2)
                        iIndex = 1 + iIndex1*pcGrid%iD2XS + iIndex2*pcGrid%iD2YS + iIndex3*pcGrid%iD2ZS
                        if (dk1k2k3 > dKcutoff) then
                            pcGrid%pacD2(iIndex) = 0.0_dp
                        end if
                    end do
                end do
            end do
        else
            !Use Halfband filter
            dHalfbandstart = 0.5_dp*dKcutoff
            dHalfbandend = 1.5_dp*dKcutoff
            do iIndex1 = 0, cSpace%iDimX - 1
                do iIndex2 = 0, cSpace%iDimY - 1
                    do iIndex3 = 0, cSpace%iDimZ - 1
                        dk1k2k3 = sqrt(dK1(1 + iIndex1)**2 + dK2(1 + iIndex2)**2 + dK3(1 + iIndex3)**2)
                        iIndex = 1 + iIndex1*pcGrid%iD2XS + iIndex2*pcGrid%iD2YS + iIndex3*pcGrid%iD2ZS
                        if (dk1k2k3 > dHalfbandend) then
                            pcGrid%pacD2(iIndex) = 0.0_dp
                        else if (dk1k2k3 >= dHalfbandstart) then
                            pcGrid%pacD2(iIndex) = pcGrid%pacD2(iIndex)* &
                                                   (0.5_dp + 9.0_dp/16.0_dp*cos(pi*((dk1k2k3 - dHalfbandstart)/dKcutoff)) &
                                                    - 1.0_dp/16.0_dp*cos(3.0_dp*pi*((dk1k2k3 - dHalfbandstart)/dKcutoff)))
                        end if
                    end do
                end do
            end do
        end if

        !inverse transform and normalize cSpace
        call TransformXYZinv(pcGrid, .true.)
        pcGrid%pacD2 = pcGrid%pacD2/(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ)

        deallocate (dK1, dK2, dK3)

    END SUBROUTINE FFilterSpaceSpatial3D

    FUNCTION dTaperingWindow(iDimT, dDt, dLeftBand, dRightBand)

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
!   The function dTaperingWindow returns a tapering window with cosined
!   tapering at start and end of the range.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iDimT            i   i8b   length of the required tapering window
!   dDt              i   dp    step size
!   dLeftBand        i   dp   width of the part to be tapered on the left side
!   dRightBand       i   dp   width of the part to be tapered on the right side
!   dTaperingWindow  o   dp   requested tapering window
!
        integer(i8b), intent(in) :: iDimT
        real(dp), intent(in) :: dDt
        real(dp), intent(in) :: dLeftBand
        real(dp), intent(in) :: dRightBand
        real(dp), dimension(iDimT) :: dTaperingWindow

! *****************************************************************************
!
!   LOCAL PARAMETERS
!
!   i       i8b   loop counter
!   ctr     i8b   loop counter
!   dTaxis   dp   axis of the variable
!
        integer(i8b) :: i, ctr
        real(dp) :: dTaxis(iDimT)

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

        dTaxis = real((/(i, i=0, iDimT - 1)/), dp)*dDt

        dTaperingWindow = 1.0_dp
        do ctr = 1, iDimT 
            if (dTaxis(ctr) < dLeftBand) then
                dTaperingWindow(ctr) = 0.5_dp + &
                                       9.0_dp/16.0_dp*cos(pi*((dLeftBand - dTaxis(ctr))/dLeftBand)) &
                                       - 1.0_dp/16.0_dp*cos(3.0_dp*pi*((dLeftBand - dTaxis(ctr))/dLeftBand)); 
            elseif (dTaxis(ctr) > (dTaxis(iDimT) - dRightBand)) then
                dTaperingWindow(ctr) = 0.5_dp + &
                                       9.0_dp/16.0_dp*cos(pi*((dTaxis(ctr) - (dTaxis(iDimT) - dRightBand))/dRightBand)) &
                                       - 1.0_dp/16.0_dp*cos(3.0_dp*pi*((dTaxis(ctr) - (dTaxis(iDimT) - dRightBand))/dRightBand)); 
            end if
        end do  

    END FUNCTION dTaperingWindow

    SUBROUTINE INTERP1DFREQ(InputVal, OutputVal, cSpace, Order)

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
        real(dp), intent(in)                    ::                  InputVal(:)
        real(dp), intent(out)                   ::                  OutputVal(:)
        integer(i8b), intent(in)                ::                  Order
        type(Space), intent(in)                 ::                  cSpace

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

        complex(dpc)                                ::                      InputValW(size(InputVal)/2 + 1)
        complex(dpc), allocatable                   ::                      OutputValW(:)
        complex(dpc), allocatable                   ::                      dMultFactor(:)
        real(dp), allocatable                       ::                      dTaperMaxFreqWindow(:)
        integer(i8b)                                ::                      InputLen, OutputLen, EvenOrOdd, iNumPlanTransform, iNumPlanTransform_inv, iDimW

        !Filtering and windowing parameters
        real(dp)                                    ::                      dLeftBand, dRightBand, dOmega(size(InputVal)/2 + 1)
        integer                                     ::                      i, iErr

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
        iDimW = min(InputLen/2 + 1, OutputLen/2 + 1)

        ! Do a FFT of InputLen-point (Initial Variable Length) with complex numbers
        call dfftw_plan_dft_r2c_1d(iNumPlanTransform, InputLen, InputVal, InputValW, FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(iNumPlanTransform, InputVal, InputValW)
        call dfftw_destroy_plan(iNumPlanTransform)

        !Allocate to upsample or decimate . OutputValW has to have more or equal points to the initial signal
        ALLOCATE (OutputValW(max(InputLen/2 + 1, OutputLen/2 + 1)))
        ALLOCATE (dTaperMaxFreqWindow(iDimW))
        ! Tapering of the max frequency because of artifact created due to higher frequencies
        dLeftBand = 0.0
        dRightBand = 0.1_dp
        dTaperMaxFreqWindow = dTaperingWindow(iDimW, 1.0_dp/(min(InputLen, OutputLen)*cSpace%dDt), dLeftBand, dRightBand)

        ! Zero-padding after the values which will be used for interpolation
        ! This is done by initializing OutputValW with 0 ( Really small number due to precision)
        ! Replace the first half part with the values of initial variable.
        OutputValW = 0.0D0
        OutputValW(1:iDimW) = InputValW(1:iDimW)/InputLen*dTaperMaxFreqWindow(1:iDimW)

        ! Inverse FFT of OutputLen-point
        call dfftw_plan_dft_c2r_1d(iNumPlanTransform_inv, OutputLen, OutputValW, OutputVal, FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(iNumPlanTransform_inv, OutputValW, OutputVal)
        call dfftw_destroy_plan(iNumPlanTransform_inv)

        DEALLOCATE (OutputValW)
        DEALLOCATE (dTaperMaxFreqWindow)
    END SUBROUTINE INTERP1DFREQ

    SUBROUTINE GREENS1DFREQ(cSpace, cSpaceTemp, SrcPos, DestPos)

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
        !   The function GREENS1DFREQ returns the convolution between the input signal
        !   (InputVal) with the Green's function.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   InputVal         r   dp   a vector of the (time) values of the initial signal
        !   Bubble_diff      r   dp   a vector of the difference between the scatterer
        !                             and the Target where we want to compute the pressure
        !   OutputVal        r   dp   a vector of the  output signal
        !   iTargetIndex     i   i8b  a vector that has the location of the target gridpoint
        !                             where we want to compute the pressure 
        !
        !
        real(dp), intent(in)            ::                  SrcPos(3), DestPos(3)
        type(Space), intent(in)         ::                  cSpace

        type(Space), intent(inout)      ::                  cSpaceTemp

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
        integer(i8b)                                ::                        iDimT, iDimW, iTargetIndex(3)
        complex(dpc)                                ::                        Green(cSpaceTemp%iDimT + 1)
        integer                                     ::                        i

        !Green's function parameters:
        integer(i8b)                                ::                        iIndZ, iOmega, iDifft
        real(dp)                                    ::                        dStartT, dEndT, dRad
        real(dp)                                    ::                        dSi1, dSi2, dCi1, dCi2; 
        real(dp)                                    ::                        dOmega,  dKcutoff, dTaperMaxFreqWindow(cSpaceTemp%iDimT+1)
        complex(dpc)                                ::                        dKangular, E1mm, E1mp, E1pm, E1pp
        integer(i8b)                                ::                        error

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

        iDimT = cSpaceTemp%iDimT
        iDimW = cSpaceTemp%iDimT/2 + 1

        iDebugLvl=0
        call ReorderDistr0ToDistr1(cSpaceTemp%cGrid)
        call TransformT(cSpaceTemp%cGrid)
        iDebugLvl=5

        dKcutoff = two_pi*(1.0_dp + 0.5_dp/real(cSpace%iDimX, dp))*cSpace%dFnyq
        iDiffT = 0; 
        dTaperMaxFreqWindow = dTaperingWindow(iDimT+1, 1.0_dp/(2*iDimT*cSpace%dDt), 0.0D0, 0.1D0)

        iTargetIndex(3) = nint(DestPos(3)/cSpace%dDx)
        iTargetIndex(1) = nint(DestPos(1)/cSpace%dDx) - iBeamOffsetX(cSpace, iTargetIndex(3))
        iTargetIndex(2) = nint(DestPos(2)/cSpace%dDx) - iBeamOffsetY(cSpace, iTargetIndex(3))

        iIndZ = mod(iTargetIndex(3) + cSpace%iDimZ - 1, cSpace%iDimZ) + cSpace%iDimZ - 1; 
        dEndT = (iDiffT + iBeamOffSetT(cSpace, iIndZ) + cSpace%iDimT + 0.5)*cSpace%dDt; 
        dStartT = (iDiffT + iBeamOffSetT(cSpace, iIndZ) - cSpace%iDimT + 0.5)*cSpace%dDt; 
       
        dRad = SQRT(SUM((SrcPos - DestPos)**2) )
        Green = 0.0; 
        if (cMediumParams%a_att == 0.0_dp) then
                do iOmega = 0, iDimT
                    dOmega = two_pi*iOmega*cSpace%dFnyq/real(iDimT, dp)
                    call cisi4((dKcutoff - dOmega)*dRad, dCi1, dSi1)
                    call cisi4((dKcutoff + dOmega)*dRad, dCi2, dSi2)
                    Green(iOmega + 1) = 1.0_dp/(dRad*4.0_dp*pi**2)*(cos(dOmega*dRad)*(dSi1 + dSi2) + sin(dOmega*dRad)*(dCi1 - dCi2 - im*pi))
                    cSpaceTemp%cGrid%pacD1(1+iOmega) = cSpaceTemp%cGrid%pacD1(1+iOmega)*Green(iOmega + 1)*cSpace%dDx**3/real(2*iDimT, dp) &
                                                       *dTaperMaxFreqWindow(iOmega+1)
                end do
        else
                do iOmega = 0, iDimT
                    dOmega = two_pi*iOmega*cSpace%dFnyq/iDimT

                    dKangular = (dOmega*cModelParams%freq0/cMediumParams%c0 &
                                + cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*dOmega*cModelParams%freq0 &
                                *(abs(dOmega*cModelParams%freq0)**(cMediumParams%b_att - 1) - &
                                (two_pi*cMediumParams%fc0)**(cMediumParams%b_att - 1)) &
                                - im*cMediumParams%alpha0_att*abs(dOmega*cModelParams%freq0)**cMediumParams%b_att)*cMediumParams%c0/cModelParams%freq0

                    call Expint_prime(-im*(dKcutoff - dKangular)*dRad, E1mm, error)
                    call Expint_prime(im*(dKcutoff - dKangular)*dRad, E1pm, error)
                    call Expint_prime(-im*(dKcutoff + dKangular)*dRad, E1mp, error)
                    call Expint_prime(im*(dKcutoff + dKangular)*dRad, E1pp, error)

                    Green(iOmega + 1) = 1.0_dp/(4.0_dp*pi*dRad)*(exp(-im*dKangular*dRad) &
                                                                + (-exp(im*dKcutoff*dRad)*(E1mm + E1mp) &
                                                                    + exp(-im*dKcutoff*dRad)*(E1Pm + E1pp))/(im*two_pi))*cSpace%dDx**3/real(2*iDimT, dp)
                    cSpaceTemp%cGrid%pacD1(1+iOmega) = cSpaceTemp%cGrid%pacD1(1+iOmega)*Green(iOmega + 1)*cSpace%dDx**3/real(2*iDimT, dp) &
                                                       *dTaperMaxFreqWindow(iOmega+1)
                end do
        end if

        iDebugLvl=0
        call TransformTInv(cSpaceTemp%cGrid)
        call ReorderDistr1ToDistr0(cSpaceTemp%cGrid)
        iDebugLvl=5 

    END SUBROUTINE GREENS1DFREQ

END MODULE ParnacTransformFilter
