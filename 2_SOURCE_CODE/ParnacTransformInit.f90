MODULE ParnacTransformInit

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
    !   The module ParnacTransformInit contains subroutines that initialize and
    !   destroy the FFTW transform plans that define the Discrete Fourier
    !   Transforms for the different dimensions utilized in the Space structure.
    !
    ! *****************************************************************************
    !
    !   MODULES USED
    !
    USE Types
    USE FFTW_cons
    USE ParnacParamDef
    USE ParnacDataDef
    USE ParnacGeneral

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
    !   GridDestroyTransforms    sub   Destroys the FFTW transform plans in a grid
    !                                  structure
    !   GridInitTransforms       sub   Initializes the FFTW transform plans in a
    !                                  grid structure
    !
    ! =============================================================================
    !
CONTAINS

    SUBROUTINE GridDestroyTransforms(cSpace)

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
        !   The subroutine GridDestroyTransforms destroys all FFTW transform plans
        !   within a grid structure - there is no check whether the plans have actually
        !   been initialized.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace   io   type(Space)   The space that contains the grid that contains
        !                               the transform structure with all the transform
        !                               plans
        !
        type(Space), intent(inout), target ::                cSpace; 
        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   pcTransforms   type(Transforms)   pointer to the Transforms structure
        !                                     within the grid within the space
        !
        type(Transforms), pointer ::        pcTransforms

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
        !   dfftw_destroy_plan
        !   PrintToLog
        !
        ! =============================================================================

        call PrintToLog("GridDestroyTransforms", 3)

        pcTransforms => cSpace%cGrid%cTransforms

        call dfftw_destroy_plan(pcTransforms%iPlanTransformT)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformT_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXY)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXY_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXYZ)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXYZ_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXY_relevantZ)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformXY_relevantZ_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformZ)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformZ_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformT_sml)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformT_sml_inv)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformT_lrg)
        call dfftw_destroy_plan(pcTransforms%iPlanTransformT_lrg_inv)

    END SUBROUTINE GridDestroyTransforms

    SUBROUTINE GridInitTransforms(cSpace)

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
        !   The subroutine GridInitTransforms initializes all FFTW transform plans
        !   within a grid structure. The grid should be initialized before this
        !   subroutine is called.
        !
        !   We use the FFTW_ESTIMATE option for planning the FFTW transforms, and
        !   therefore we don't need to allocate the data arrays in the grid for the
        !   planning.
        !
        ! *****************************************************************************
        !
        !   INPUT/OUTPUT PARAMETERS
        !
        !   cSpace   io   type(Space)   The space that contains the grid that contains
        !                               the transform structure with all the transform
        !                               plans
        !
        type(Space), intent(inout), target ::                cSpace; 
        ! *****************************************************************************
        !
        !   LOCAL PARAMETERS
        !
        !   pcGrid           type(Grid)   pointer to the Grid structure within cSpace
        !   cFFTWtestarray   dpc          test array that is passed as dummy argument
        !                                 to the planning routines
        !
        type(Grid), pointer :: pcGrid
        complex(dpc) :: cFFTWtestarray(1)
        integer      :: iret

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
        !   dfftw_plan_guru_dft_r2c
        !   dfftw_plan_guru_dft_c2r
        !   dfftw_plan_guru_dft
        !   PrintToLog
        !
        ! =============================================================================

        call PrintToLog("GridInitTransforms", 3)

        pcGrid => cSpace%cGrid

        ! The forward and inverse FFT on the temporal dimension for a Distr. 1
        ! data array, i.e. a real data array that has been stored in a
        ! complex data type - this is somewhat awkward but it is done to enable
        ! the usage of in-place transforms only whil benefiting from the
        ! memory reduction in using real-to-complex transforms

        !if (pcGrid%iD0TL==0) then !L.D. 11-11-2011
        !pcGrid%iD0TL=1
        !end if

        !KH For periodical T, change here..
        !   print *, "pcGrid%iD0TL"
        !   print *, pcGrid%iD0TL
        !   print *, int((/ 2*pcGrid%iD0TL /))
        !   print *, int((/ 1+2*pcGrid%iD0TL /))
        !   print *, int((/ 2*pcGrid%iD1TL /))
        !   !print *, int((/ 1+2*pcGrid%iD0TL /))
        !   !print *, int((/ 2+2*pcGrid%iD0TL /))
        !   print *, "pcGrid%iD1TS"
        !   print *, pcGrid%iD1TS
        !   print *, "pcGrid%iD1LocN"
        !   print *, pcGrid%iD1LocN
        !   print *, "pcGrid%iD1TL"
        !   print *, pcGrid%iD1TL

        call dfftw_init_threads(iret)
        call dfftw_plan_with_nthreads(mkl_get_max_threads())
        
        ! Change this to 2*pcGrid%iD0TL both in T_inv for  2PPW
        call dfftw_plan_guru_dft_r2c(pcGrid%cTransforms%iPlanTransformT, &
                                     1, &
                                     int((/2*pcGrid%iD0TL/)), &        ! We take the length of the T-axis in distribution 1, the factor 2 is due to the fact that the input is real and the array in distribution 2 is complex
                                     int((/pcGrid%iD1TS/)), &        !KH no: the factor 2 is due to the wraparound region
                                     int((/pcGrid%iD1TS/)), &        ! in complex, we have floor((2*iD0TL)/2)+1 numbers
                                     1, &        ! for periodical T, we would need floor(iD0TL/2)+1 numbers in dist. 1
                                     int((/pcGrid%iD1LocN/)), &
                                     int((/2*pcGrid%iD1TL/)), &        !this factor two is merely because of storing real values in complex variables
                                     int((/pcGrid%iD1TL/)), &
                                     cFFTWtestarray, &        !empty array of size 1; for FFTW_ESTIMATE, the array isn't used,
                                     cFFTWtestarray, &        !but we have to throw the subroutine a bone...
                                     FFTW_ESTIMATE);
 
        !KH For periodical T, change here..
        call dfftw_plan_guru_dft_c2r(pcGrid%cTransforms%iPlanTransformT_inv, &
                                     1, &
                                     int((/2*pcGrid%iD0TL/)), &
                                     int((/pcGrid%iD1TS/)), &
                                     int((/pcGrid%iD1TS/)), &
                                     1, &
                                     int((/pcGrid%iD1LocN/)), &
                                     int((/pcGrid%iD1TL/)), &
                                     int((/2*pcGrid%iD1TL/)), &
                                     cFFTWtestarray, &
                                     cFFTWtestarray, &
                                     FFTW_ESTIMATE); 
        ! The forward and inverse transforms on the spatial dimensions for a
        ! Distr. 2 data array. All transforms in X/Y/Z are assumed to be
        ! performed on one XYZ-slice only on each processor, so iDimT=1
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXY, &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 int((/pcGrid%iD2ZL/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_FORWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXY_Inv, &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 int((/pcGrid%iD2ZL/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_BACKWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXYZ, &
                                 3, &
                                 int((/pcGrid%iD2ZL, pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 1, &
                                 1, &
                                 1, &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_FORWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXYZ_Inv, &
                                 3, &
                                 int((/pcGrid%iD2ZL, pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 1, &
                                 1, &
                                 1, &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_BACKWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformZ, &
                                 1, &
                                 int((/pcGrid%iD2ZL/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_FORWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformZ_Inv, &
                                 1, &
                                 int((/pcGrid%iD2ZL/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_BACKWARD, FFTW_ESTIMATE); 
        ! for ConvolutionSlice spaces, the length of the relevant range of
        ! positions in the Z-dimension is equal to the iDimZ+iBeamindexstartZ,
        ! which is much smaller than pcGrid%iD2ZL - saves some unnecessary
        !transforms
        !    print *, "pcGrid%iD2XL"
        !    print *, pcGrid%iD2XL
        !    print *, "pcGrid%iD2YL"
        !    print *, pcGrid%iD2YL
        !    print *, "pcGrid%iD2XS"
        !    print *, pcGrid%iD2XS
        !    print *, "pcGrid%iD2YS"
        !    print *, pcGrid%iD2YS
        !    print *, "cSpace%iDimZ + cSpace%iBeamindexstartZ"
        !    print *, cSpace%iDimZ + cSpace%iBeamindexstartZ
        !    print *, "pcGrid%iD2ZS"
        !    print *, pcGrid%iD2ZS

        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXY_relevantZ, &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 int((/cSpace%iDimZ + cSpace%iBeamindexstartZ/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_FORWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXY_relevantZ_Inv, &
                                 2, &
                                 int((/pcGrid%iD2XL, pcGrid%iD2YL/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 int((/cSpace%iDimZ + cSpace%iBeamindexstartZ/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 int((/pcGrid%iD2ZS/)), &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_BACKWARD, FFTW_ESTIMATE); 
        ! Forward and inverse transforms on the temporal dimension. The length of
        ! the time axis is smaller than the normal temporal FFT's - used in the
        ! nonlinear operator

        !KH For periodical T, change here..
        call dfftw_plan_guru_dft_r2c(pcGrid%cTransforms%iPlanTransformT_sml, &
                                     1, &
                                     int((/pcGrid%iD0TL/)), &        ! We take the length of the T-axis in distribution 0
                                     int((/pcGrid%iD1TS/)), &        ! That is about half the total length of reals in dist. 1
                                     int((/pcGrid%iD1TS/)), &        ! in complex, we have floor(iD0TL/2)+1 numbers
                                     1, &        ! however, pacD1 is larger, dim in T is iD1TL = iD0TL+1
                                     int((/pcGrid%iD1LocN/)), &
                                     int((/2*pcGrid%iD1TL/)), &        !this factor two is merely because of storing real values in complex variables
                                     int((/pcGrid%iD1TL/)), &
                                     cFFTWtestarray, &
                                     cFFTWtestarray, &
                                     FFTW_ESTIMATE);
 
        !KH For periodical T, change here..
        call dfftw_plan_guru_dft_c2r(pcGrid%cTransforms%iPlanTransformT_sml_inv, &
                                     1, &
                                     int((/pcGrid%iD0TL/)), &
                                     int((/pcGrid%iD1TS/)), &
                                     int((/pcGrid%iD1TS/)), &
                                     1, &
                                     int((/pcGrid%iD1LocN/)), &
                                     int((/pcGrid%iD1TL/)), &
                                     int((/2*pcGrid%iD1TL/)), &
                                     cFFTWtestarray, &
                                     cFFTWtestarray, &
                                     FFTW_ESTIMATE); 
        ! Forward and inverse FFT's on the spatial dimensions, employing smaller
        ! axes than normal - used in the inhomogeneity contrast operator

        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXYZ_sml, &
                                 3, &
                                 int((/pcGrid%iD2ZL/2, pcGrid%iD2XL/2, pcGrid%iD2YL/2/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 1, &
                                 1, &
                                 1, &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_FORWARD, FFTW_ESTIMATE); 
        call dfftw_plan_guru_dft(pcGrid%cTransforms%iPlanTransformXYZ_sml_Inv, &
                                 3, &
                                 int((/pcGrid%iD2ZL/2, pcGrid%iD2XL/2, pcGrid%iD2YL/2/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 int((/pcGrid%iD2ZS, pcGrid%iD2XS, pcGrid%iD2YS/)), &
                                 1, &
                                 1, &
                                 1, &
                                 1, &
                                 cFFTWtestarray, &
                                 cFFTWtestarray, &
                                 FFTW_BACKWARD, FFTW_ESTIMATE); 
    END SUBROUTINE GridInitTransforms

END MODULE ParnacTransformInit
