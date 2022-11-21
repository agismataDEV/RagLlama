
MODULE ParnacDataDef

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
!   The module ParnacDataDef contains the type definition of the data structure
!   that forms the basis of the Parnac program: the Space structure. The Space
!   structure describes and contains a quantity in a 4-dimensional spatio-
!   temporal domain (t,x,y,z). This quantity may amongst others be an acoustic
!   pressure, a contrast source, a primary source excitation...
!
!   The spatiotemporal domain may be parallel-shaped in the XZ-, YZ- and/or
!   TZ-plane. The skewness is described by the angles dThetaT, dThetaX and
!   dThetaY, being the angle in radians between the T/X/Y axis and that side
!   of the parallel that for a rectangular-shaped domain would lie along the
!   Z-axis. E.g. an angle dThetaX=pi/2 means that in the XZ-plane we have a
!   normal rectangular-shaped domain. All angles are between 0 and pi.
!   In the Space structure, the angles are represented with the values
!   dTanT = 1/tan(dThetaT); 0 if dThetaT = pi/2
!   dTanX = 1/tan(dThetaX); 0 if dThetaX = pi/2
!   dTanY = 1/tan(dThetaY); 0 if dThetaY = pi/2
!   The domain is further determined by the lengths (dLt, dLx, dLy, dLz)
!   and the starting positions (dSt, dSx, dSy, dSz) in all dimensions. The
!   lengths are set in the configuration file. The starting position in z, dSz
!   may also be set in the configuration file, but it should normally be put
!   to zero. However, if dSz is zero, in the procedure
!   SpecialChecksConfigParameters it is reset to the first layer beyond the
!   source plane. The starting positions in t/x/y are determined from their
!   starting position in the plane Z=0, and then corrected with the angles
!   dThetaT, dThetaX and dThetaY. At Z=0, the T axis starts at T=0 and the
!   X and Y axes are centered around the origin. With dTanT, dTanX and dTanY,
!   we therefore have
!   dSt = 0 + dSz * dTanT
!   dSx = -dLx/2 + dSz * dTanX
!   dSy = -dLy/2 + dSz * dTanY
!
!   The spatiotemporal domain is discretized on an equidistant grid determined
!   by the step sizes dDt and dDx, which are derived from the Nyquist Frequency
!   Fnyq. In the grid, we now distinguish between two different sets of four
!   indices to denote a position:
!   - The global index (iGt, iGx, iGy, iGz): the spatiotemporal coordinate of
!       the point divided by the stepsize dDt or dDx.
!   - The beam index (iBt, iBx, iBy, iBz): index within the beam relative to
!       the first position in the beam, which has index 0. For indices in the
!       T/X/Y dimension, even with a skewed beam, the index always starts at 0
!       on the first position.
!   The dimension lengths are derived from the physical lengths with
!   e.g. iDimT = ceiling(dLt/dDt). The global indices of the starting positions
!   are obtained as
!   iStartZ = nint(dSz/dDx)
!   iStartT = 0 + nint(iStartZ * dTanT)
!   iStartX = -floor(iDimX/2) + nint(iStartZ * dTanX)
!   iStartY = -floor(iDimX/2) + nint(iStartZ * dTanY)
!   To translate a beam index in a global index and vice versa, we also need
!   the offset in the T/X/Y dimensions at each z-position due to the skewness
!   of the beam. This offset is given by the functions iBeamOffsetT/X/Y(iBz),
!   which determine the offsets in T/X/Y as a function of the BEAM index iBz
!   (so not as a function of the GLOBAL index iGz). The relations between the
!   global index and the beam index of a point in the 4D grid are then:
!   iBz = iGz - iStartZ
!   iBt/x/y = iGt/x/y - iStartT/X/Y - iBeamOffsetT/X/Y(iBz)
!   and vice versa:
!   iGz = iBz + iStartZ
!   iGt/x/y = iBt/x/y + iStartT/X/Y + iBeamOffsetT/X/Y(iBz)
!
!   In the Space structure, the following of the already discussed domain-
!   related parameters are stored:
!   - dTanT/X/Y       - describing the beam skewness
!   - dFnyq, dDt, dDx - the discretization frequency and step sizes
!   - iDimT/X/Y/Z     - the discrete lengths
!   - iStartT/X/Y/Z   - the discrete starting positions
!
!   Besides these parameters, in the Space structure we also have
!   - iSpaceIdentifier      - gives the type of the stored quantity; see
!                             ParnacIdentifierDef for the possible values.
!   - bYsymm                - whether the stored space exhibits symmetry
!                             in the Y-dimension. See footnote 1)
!   - iBeamIndexStartX/Y/Z  - have only significance for the convolution
!                             slices, i.e. if the stored quantity is of
!                             type iSI_XYZSLICE. See footnote 2)
!
!   Apart from these descriptive parameters, the Space structure contains the
!   Grid structure, which actually stores the values of the quantity on all
!   positions in the 4D grid. To understand how it is stored, we need to
!   consider the following:
!   - The entire grid is distributed over the local memory on all processors.
!   Some operations require the full temporal axis on a number of positions to
!   be available on a processor; this is called the T-local distribution. Some
!   operations require the full spatial domain on a number of time instants to
!   be available on a processor; this is called the XYZ-local distribution.
!   - The quantity may be either in the temporal domain, in the temporal
!   Fourier domain or in the spatiotemporal Fourier domain. In the first
!   domain, it is a real variable (stored as dp), in the latter two domains,
!   it is a complex variable (stored as dpc).
!   - In the temporal Fourier domain, we may or may not require the data to be
!   transformed with a zero-padding region in the temporal
!   dimension.
!
!   In the program, we employ the following three forms of storage for the
!   variable in the space:
!   - Distribution 0. The data is in the T-local distribution, and it is in the
!     temporal domain. The data format is double precision (dp). The total size
!     of the array stored on each processor is (iDimT*iDimX*iDimY*iDimZ)/iProcN
!     (iProcN is the total number of processors).
!   - Distribution 1. The data is in the T-local distribution, and the INTENT
!     (see footnote 3) of the distribution is that the data is in the temporal
!     Fourier domain and that it includes a zero-padding region in the temporal
!     dimension. The data format is double complex (dpc). The total size of the
!     array stored on each processor is [(iDimT+1)*iDimX*iDimY*iDimZ]/iProcN.
!     Although the temporal dimension is doubled due to the zero-padding
!     region, it is again almost halved because of the removal of the negative
!     frequencies employed in the FFTW real-to-complex transforms.
!   - XYZ-local distribution, either temporal Fourier or spatiotemporal Fourier
!     domain. The data format is double precision complex (dpc). NB: in this
!     context, reference to a 'time instant' would signify a 'frequency'.
!   In the Grid structure, we define three pointers to allocatable arrays for
!   the three different distributions: parD0, pacD1 and pacD2. Moreover, an
!   identifier iDistr denotes the distribution that the quantity is in. The
!   subroutines in ParnacDataRedist are employed to change the data from one
!   distribution to another.
!
!   For each of the distributions, the entire grid is stored as a
!   one-dimensional array, even though the grid is four-dimensional. This gives
!   us more flexibility and control over the way it is stored. However, as we
!   have to take care of the striding and of the array management by ourselves,
!   this makes the code less transparent and it may easier result in error. In
!   the Grid structure, for each distribution there is a number of descriptive
!   parameters:
!   - iD0/1/2LocSize   - the number of elements on this processor
!   - iD0/1/2GlobN     - the total number of XYZ-positions (for T-local) or
!                        T-instants (for XYZ-local)
!   - iD0/1/2LocN      - the number of XYZ-positions (for T-local) or
!                        T-instants (for XYZ-local) on this proc
!   - iD0/1/2T/X/Y/ZS  - the strides in the various dimensions
!   - iD0/1/2IS        - the stride in the XYZ-positions (for T-local)
!                        or T-instants (for XYZ-local)
!   - iD0/1/2T/X/Y/ZL  - the lengths of the various dimensions
!   - aiD0/1/2Loc      - the beam indices of the XYZ-positions (for T-local)
!                        or T-instants (for XYZ-local) stored on this processor
!
!   Other parameters in the Grid structure are:
!   - iProcN           - Total number of processors
!   - iProcID          - ID of this processor
!   and of course the iDistr parameter denoting the distribution of the data,
!   and the pointers parD0, pacD1 and pacD2. In the Grid structure, we also
!   store the cTransforms structure that contains the information for the
!   various FFT transforms performed on this quantity.

!   Footnotes:
!   1) Symmetry in Y. In the program, we have a special option for quantities
!   exhibiting a symmetry in the Y-dimension. If this is the case, the flag
!   bYsymm is set to .true., only the positive half of the Y-dimension is
!   stored, iStartY is set to 0 and iDimY is halved.
!   2) BeamIndexStartX/Y/Z. These parameters are only of significance if the
!   Space type is iSI_CONVOLUTIONSLICE. In this case, the Space is used to
!   store an XYZ-slice of one single frequency that is convolved in space
!   with another slice. For the convolution, the beam index may be negative,
!   and the negative values are stored at the end of each dimension, i.e.
!   for a range in the beam index of -3...5 we get the order
!   [0 1 2 3 4 5 -3 -2 -1]. For this case, BeamIndexStart would be -3. If we
!   introduce the array index iAi to denote the number in such an order
!   (starting at an index 0), then the relation between the array index and the
!   beam index iBi is
!   iBi = mod(iAi - iBeamIndexStart, iDim) + iBeamIndexStart.
!   3) Distribution 1 may also contain the variable in the temporal dimension
!   including the zero-padding region. In this case, the real data is stored
!   in a complex variable in the order [real(1) real(2) real(3) real(4) ...] =
!   [(real,imaginary)(1) (real,imaginary)(2) ...]. - This is used when changing
!   from Dist. 0 to Dist. 1, just before the application of the FFT. It is also
!   employed in the evaluation of the Z derivative in BeamDerivativeZ.
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types
    USE ParnacTransformDef

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   Grid          type  see the DESCRIPTION header of the module
!   Space         type  see the DESCRIPTION header of the module
!   RRMSNormtype  type  contains the Relative Root Mean Square (RRMS) norm that
!                       gives a measure of the iterative behavior of the scheme
!
! =============================================================================
!
    IMPLICIT NONE

    ! A type that contains the variable and all information on how it is
    ! stored. The information on the quantity it represents is stored in the
    ! type Space.
    type Grid

        integer(i8b) ::                        iDistr; ! Identifier for the current distribution: is -1 if there is no data
        integer(i8b) ::                        iProcN; ! Total number of processors
        integer(i8b) ::                        iProcID; ! ID number of the current processor

        ! Information on Distribution type 0
        integer(i8b) ::                        iD0LocSize; ! The amount of elements we store locally in distribution 0
        integer(i8b) ::                        iD0GlobN; ! The total number of XYZ-positions
        integer(i8b) ::                        iD0LocN; ! The number of XYZ-positions we store locally
        integer(i8b) ::                        iD0TS, iD0IS                ! The stride of this distribution with respect to t and the index i
        integer(i8b) ::                        iD0TL                                ! The number of elements in the t-axis
        integer(i8b), allocatable ::        aiD0Loc(:, :)        ! cGrid.D0Loc(1+i,:) contains the (x,y,z) coordinates of
        ! each XYZ-position i which is stored locally on this processor.

        ! Information on Distribution type 1
        integer(i8b) ::                        iD1LocSize; ! The amount of elements we store locally in distribution 1
        integer(i8b) ::                        iD1GlobN; ! The total number of XYZ-positions
        integer(i8b) ::                        iD1LocN; ! The number of XYZ-positions we store locally
        integer(i8b) ::                        iD1TS, iD1IS                ! The stride of this distribution with respect to t and the index i
        integer(i8b) ::                        iD1TL                                ! The number of elements in the t-axis
        integer(i8b), allocatable ::        aiD1Loc(:, :)        ! cGrid.D1Loc(1+i,:) contains the (x,y,z) coordinates of
        ! each XYZ-position i which is stored locally on this processor.

        ! Information on Distribution type 2
        integer(i8b) ::                        iD2LocSize; ! The amount of elements we store locally in distribution 2
        integer(i8b) ::                        iD2GlobN; ! The total number of T-instants
        integer(i8b) ::                        iD2LocN; ! The number of T-instants we store locally
        integer(i8b) ::                        iD2XS, iD2YS, iD2ZS, iD2IS        ! The stride of this distribution with respect to X, Y, Z and the iIndex
        integer(i8b) ::                        iD2XL, iD2YL, iD2ZL        ! The number of elements in the x,y, z-axis
        integer(i8b), allocatable ::        aiD2Loc(:)        ! cGrid.D2Loc(1+i) contains the (t) coordinate of
        ! each (x,y,z)-block i which is stored locally on this processor.

        ! Pointers to the data, stored in a one-dimensional array
        real(dp), pointer ::                parD0(:)
        complex(dpc), pointer ::        pacD1(:)
        complex(dpc), pointer ::        pacD2(:)

        ! The transform plans for FFTW
        type(Transforms) ::                cTransforms; 
    end type Grid

    ! A type that contains a quantity in a four-dimensional domain and all
    ! information regarding size, starting position, discretization, skewness..
    ! The data itself is stored in the type Grid.
    type Space

        ! Descriptive variables:
        integer(i8b) ::                iSpaceIdentifier                        ! Identification number for the type of space, see ParnacIdentifierDef
        logical ::                        bYSymm                                                ! Is it symmetric in Y?
        ! Then we can store it more efficiently, and use a more efficient Convolution

        !Angles of the space, denoted as the tangents = 1/tan(theta)
        real(dp) ::                        dTanT, dTanX, dTanY

        !Sampling frequency and step sizes
        real(dp) ::                 dFnyq, dDt, dDx

        !Size of the beam in points
        integer(i8b) ::         iDimT, iDimX, iDimY, iDimZ

        !Start of the global index in points, relative to the skew X and Y axes
        integer(i8b) ::                iStartT, iStartX, iStartY, iStartZ

        !Start of the beam index in points, see header remarks for explanation
        integer(i8b) ::                iBeamIndexStartX, iBeamIndexStartY, iBeamIndexStartZ

        ! Information on the grid
        type(Grid) ::                        cGrid; 
    end type Space

    !This type is defined to keep track of the RRMS norm in the different
    ! iterations. The norm is based on the last z slice in the beam under
    ! calculation and is obtained from the current field correction and
    ! the previous one.
    type RRMSNormtype

        real(dp), allocatable :: arRRMSNorm(:, :)                !Only for proc 0: the norms for all beams and iterations
        real(dp), allocatable :: arPreviouscorrection(:)!This array stores the field correction of the previous iteration on this proc
        integer(i8b) ::                         iStart                                        !Starting beam index of the last z-plane on this processor; is -1 if nothing is here

    end type RRMSNormtype

END MODULE ParnacDataDef
