MODULE ParnacDerivative

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
!   The module ParnacDerivative contains subroutines for the evaluation of 
!   finite difference approximations of arbitrary order in an arbitrary
!   dimension. -- The approach employed in this module (with the DerivLookup 
!   structure) is rather heavy compared to the usage in the Parnac Program, but 
!   it works fine.
!
! *****************************************************************************
!
!   MODULES USED
!
	USE Types
	USE ParnacGeneral
	USE ParnacParamDef
	USE ParnacDataDef
    USE ParnacIdentifierDef

	USE FDWeight, ONLY : &
        FDWEIGHTS
    USE ParnacDataSpaceInit, ONLY : &
        InitSpace, InitGrid, DestructSpace, &
        GridDistr2CreateEmpty
	USE ParnacDataMap, ONLY : &
        MapVt, MapVtInv, &
        MapVx, MapVxInv, &
        MapVy, MapVyInv
    USE ParnacDataRedist, ONLY : &
        ReorderDistr0ToDistr1,ReorderDistr1ToDistr0, &
        ReorderDistr1ToDistr2,ReorderDistr2ToDistr1, &
        Distr2ObtainXYZBlock, Distr2PutXYZBlock

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   iFDMinOrder   i8b   Minimum order of the finite difference scheme, to which
!                       it reduces as the border is coming nearer and nearer
!   DerivLookup   type  Contains the stencils for varying stencil sizes and 
!                       all possible centered or skew finite difference 
!                       approximations. The structure contains the following 
!                       variables:
!     - iOrder    i8b   maximum order of the FD approximation ( = size of the 
!                       stencil - 1), up to which the weights and positions
!                       have been determined
!     - arWeights  dp   arWeights(iOrder, iTarget, iInd) is the weight of 
!                       the position iInd within a FD stencil of order iOrder 
!                       and with the target point iTarget (= the point where 
!                       the FD approximation is obtained). For example, if the
!                       first derivative d/dt is used (determined by 
!                       bSecondOrder =.false. in DerivLookupInit), then we have 
!                       for iOrder = 2 the following arWeights(2,iTarget,iInd):
!
!                       iTarget \  iInd
!                                 \    1      2      3
!                           1        -1.5    2.0   -0.5  - right sided
!                           2         1.0   -2.0    1.0  - centered
!                           3         0.5   -2.0    1.5  - left sided
!
!                       If iOrder>iFDMinOrder, then the more skewed the stencil
!                       is, the more the order is reduced down to iFDMinOrder.
!                       This is done to prevent numerical errors from blowing
!                       up, as with the increasing order the stencil weights 
!                       become extremely large for skew stencils.
!     - aiPoints  i8b   aiPoints(iOrder, iTarget, 1:3) contains information on
!                       the location of the stencil relative to the target
!                       point. For each iOrder and iTarget, the first index 
!                       holds the starting point for the stencil relative to 
!                       the target point, the second index contains the end 
!                       point relative to the target point, and the third index
!                       contains the length of the stencil. For the example 
!                       above, we have aiPoints(2,iTarget,1:3) = 
!
!                       iTarget \ start(1) end(2) length(3)
!                           1       0        2       3      - right sided
!                           2      -1        1       3      - centered
!                           3      -2        0       3      - left sided
!
IMPLICIT NONE

	integer(i8b), parameter::			iFDMinOrder	= 6; 

	type DerivLookup
		integer(i8b)::				iOrder
		real(dp), allocatable ::		arWeights(:,:,:);
		integer(i8b), allocatable ::	aiPoints(:,:,:);
	end type

! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   DerivLookupInit      sub   Initialize the DerivLookup structure that 
!                              contains weight stencils used for the 
!                              evaluation of the Finite Difference
!   DerivLookupDestroy   sub   Destroys the DerivLookup structure
!   DerivativeReal       sub   Evaluates the finite difference for a real 
!                              variable
!   DerivativeComplex    sub   Evaluates the finite difference for a complex
!                              variable
!   BeamDerivativeZ      sub   Computes the finite difference in the z-
!                              dimension for a space containing a (skew) beam
!
! =============================================================================
!
CONTAINS

SUBROUTINE DerivLookupInit(iFDOrder, iMinOrder, cDerivLookup, bSecondOrder)

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
!   The subroutine DerivLookupInit fills a DerivLookup structure with all
!   information needed to compute a finite difference approximation with a
!   stencil of varying order and of varying skewness. 
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   iFDOrder    i   i8b   The maximum order of the FD scheme we use
!   iMinOrder   i   i8b   The order of the FD scheme near the borders
!   cDerivLookup io type(DerivLookup)   DerivLookup structure to be initialized
!   bSecondOrder i  lgt   False if the weights of the first order derivative 
!                         d/dt are required. True if the weights of the second
!                         order derivative d^2/dt^2 are required.
!
	integer(i8b), intent(in) ::     iFDOrder, iMinOrder
	type(DerivLookup), intent(inout) ::				cDerivLookup;
	logical(lgt), intent(in) ::					bSecondOrder;
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iInd       i8b   index pointing to the target point
!   iW         i8b   temporary variable 
!   iLOrder    i8b   order of the FD approximation
!   iOrder     i8b   real order of the FD approximation (always <= iLOrder)
!   iLeft      i8b   leftmost index
!   iRight     i8b   rightmost index
!   iMid       i8b   index of the target point relative to left and right
!   arPoints    dp   points at which the weights are determined
!   arWeights   dp   weights output by the FDWeights subroutine
!   iCurFDOrder i8b  order of the derivative +1 (either 2 for d/dt or 3 for 
!                    d^2/dt^2)
! 
	integer(i8b)::                  iOrder, iInd, iW, iLOrder, iLeft, iRight, iMid
	real(dp)::                      arPoints(iFDOrder+1)
	real(dp)::                      arWeights(iFDOrder+1,3)
	integer(i8b)::					iCurFDOrder;

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
!   DerivLookupDestroy
!   FDWeights
!   
! =============================================================================
	
	call PrintToLog("DerivLookupInit", 4)

	! First we will destroy the old data
	call DerivLookupDestroy(cDerivLookup);
	cDerivLookup.iOrder	= iFDOrder;

	! Initialize the data
	allocate(cDerivLookup.arWeights(iFDOrder, iFDOrder+1, iFDOrder+1));
	allocate(cDerivLookup.aiPoints(iFDOrder, iFDOrder+1, 3));
	cDerivLookup.arWeights=0.0_dp
	cDerivLookup.aiPoints=0

	if (bSecondOrder) then
		iCurFDOrder	= 3;
	else
		iCurFDOrder	= 2;
	end if

	! We will compute these lists for two separate cases:
	! 1- The case iOrder <= iMinOrder 
    !       always use iOrder, not depending on the skewness of the stencil
	! 2- The case iMinOrder < iOrder <= iFDOrder 
    !       use iOrder for the centered stencil, but for the skew stencils
    !       reduce the order stepwise down to iMinOrder for the completely
    !       left- and rightsided stencils
    
	! 1- The case iLOrder <= iMinOrder
	do iLOrder = 1, iMinOrder
	   do iInd = 1, iLOrder+1
	      arPoints(1:iLOrder+1)     = (/ (real(iW-iInd, dp), iW=1,iLOrder+1) /)
	      call FDWeights(real(0, dp), arPoints, int(iLOrder), int(iCurFDOrder,i4b), arWeights(1:iLOrder+1,:))	
	      cDerivLookup.arWeights(iLOrder, iInd, 1:iLOrder+1)= arWeights(1:iLOrder+1,iCurFDOrder)
	      cDerivLookup.aiPoints(iLOrder, iInd, 1) = 1 - iInd
	      cDerivLookup.aiPoints(iLOrder, iInd, 2) = iLOrder+1 - iInd
	      cDerivLookup.aiPoints(iLOrder, iInd, 3) = cDerivLookup.aiPoints(iLOrder, iInd, 2) - cDerivLookup.aiPoints(iLOrder, iInd, 1) + 1;
	   end do
	end do

	! 2- The case iMinOrder < iOrder <= iFDOrder
	do iLOrder = iMinOrder+1, iFDOrder
	   do iInd = 1, iLOrder+1
	      ! First decide which order to use here
	      if (iInd-1 < (iLOrder+1-iInd)) then
			 iOrder = (iInd-1)*2+iFDMinOrder;
			 if (iOrder > iLOrder) then
				iOrder = iLOrder;
			 end if

			 iLeft  = 1;
			 iMid   = iLeft + (iInd - 1)
			 iRight = iLeft+iOrder;
         else
			 iOrder = (iLOrder+1-iInd)*2+iFDMinOrder;
			 if (iOrder > iLOrder) then
				iOrder = iLOrder;
			 end if
			 
			 iRight = iLOrder+1;
			 iMid   = iRight - (iLOrder+1-iInd)
			 iLeft  = iRight-iOrder;
	      end if

	      ! Create the points
	      arPoints(1:iOrder+1)     = (/ (real(iW-iMid, dp), iW=iLeft,iRight) /)
	      call FDWeights(real(0, dp), arPoints, int(iOrder), int(iCurFDOrder,i4b), arWeights(1:iOrder+1,:))	
	      cDerivLookup.arWeights(iLOrder, iInd, 1:iOrder+1)= arWeights(1:iOrder+1,iCurFDOrder)
	      cDerivLookup.aiPoints(iLOrder, iInd, 1) = iLeft - iInd
	      cDerivLookup.aiPoints(iLOrder, iInd, 2) = iRight - iInd
	      cDerivLookup.aiPoints(iLOrder, iInd, 3) = cDerivLookup.aiPoints(iLOrder, iInd, 2) - cDerivLookup.aiPoints(iLOrder, iInd, 1) + 1;
	   end do
	end do

END SUBROUTINE DerivLookupInit

SUBROUTINE DerivLookupDestroy(cDerivLookup)

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
!   The subroutine DerivLookupDestroy destroys the DerivLookup structure.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cDerivLookup   io type(DerivLookup)   The DerivLookup structure to be 
!                                         destroyed
!
	type(DerivLookup), intent(inout)::		cDerivLookup
	
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
	
	! Deallocate each array inside
	if (allocated(cDerivLookup.arWeights)) then
	   deallocate(cDerivLookup.arWeights);
	end if
	if (allocated(cDerivLookup.aiPoints)) then
	   deallocate(cDerivLookup.aiPoints);
	end if
END SUBROUTINE DerivLookupDestroy

SUBROUTINE DerivativeReal(arIn, arOut, iLength, rDelta, arWeights, aiPoints)

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
!   The subroutine DerivativeReal computes the Finite Difference approximation
!   on an array of real data. The length of the array should be larger than
!   the order of the FD approximation - if not, an FD scheme of reduced order
!   is used.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   arIn    i   dp    Array on which to compute the FD approximation
!   arOut   o   dp    Result
!   iLength i   i8b   Length of the arrays
!   rDelta  i   dp    Step size
!   arWeights   i   dp   Array containing the weights of the stencils
!   aiPoints    i   i8b  Array containing the positions of the stencils
!
	real(dp), intent(in)::		arIn(:)
	real(dp), intent(out)::		arOut(:);
	integer(i8b), intent(in)::	iLength;
	real(dp), intent(in)::		rDelta
	real(dp), intent(in)::		arWeights(:, :, :)
	integer(i8b), intent(in)::	aiPoints(:, :, :);
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL          i8b   Position within the array
!   iOrder      i8b   Order of the FD scheme
!   iOrderhalf  i8b   Half of it
!   iIndex      i8b   Position within the FD stencil 
! 
	integer(i8b)::			iL, iOrder, iOrderhalf, iIndex;

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
			
	call PrintToLog("DerivativeReal",7)

	iOrder 	= size(arWeights, 1);

	arOut = 0;

	if ((iLength >= 1)) then
		
		if (iLength - 1 < iOrder) then
		!if Length of input is smaller than maximum order in arWeights, 
		! then use reduced order weights
			do iL = 1, iLength
				arOut(iL) = &
				sum(arIn(iL+aiPoints(iLength-1, iL,1):iL+aiPoints(iLength-1, iL,2)) * &
					arWeights(iLength-1, iL, 1:aiPoints(iLength-1,iL,3)))/rDelta;
			end do
		else
		!otherwise, use maximum order, but with reduced order at boundaries 
		!(as already included in DerivLookupInit)
			iOrderhalf = floor(real(iOrder)/2.0)
			iIndex	= 1; !iIndex keeps track of the position inside the FD stencil

			! The left part
			do iL = 1, iOrderhalf
				arOut(iL) = &
				sum(arIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
				iIndex 	= iIndex + 1;
			end do
			
			! The middle part
			do iL = iOrderhalf + 1, iLength-iOrderhalf
				arOut(iL) = &
				sum(arIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
			end do
			
			iIndex	= iIndex + 1;
			if (modulo(iOrder, 2_i8b) == 1) then
				iIndex	= iIndex + 1;
			end if

			! The right part
			do iL = iLength-iOrderhalf+1, iLength
				arOut(iL) = &
					sum(arIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
				iIndex 	= iIndex + 1;
			end do
		end if

	end if	
END SUBROUTINE DerivativeReal

SUBROUTINE DerivativeComplex(acIn, acOut, iLength, rDelta, arWeights, aiPoints)

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
!   The subroutine DerivativeComplex computes the Finite Difference approximation
!   on an array of complex data. The length of the array should be larger than
!   the order of the FD approximation - if not, an FD scheme of reduced order
!   is used.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   arIn    i   dpc   Array on which to compute the FD approximation
!   arOut   o   dpc   Result
!   iLength i   i8b   Length of the arrays
!   rDelta  i   dp    Step size
!   arWeights   i   dp   Array containing the weights of the stencils
!   aiPoints    i   i8b  Array containing the positions of the stencils
!
	complex(dpc), intent(in)::	acIn(:)
	complex(dpc), intent(inout)::	 acOut(:);
	integer(i8b), intent(in)::	iLength;
	real(dp), intent(in) ::			rDelta
	real(dp), intent(in)::	arWeights(:, :, :)
	integer(i8b), intent(in)::	aiPoints(:, :, :);
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   iL          i8b   Position within the array
!   iOrder      i8b   Order of the FD scheme
!   iOrderhalf  i8b   Half of it
!   iIndex      i8b   Position within the FD stencil 
! 
	integer(i8b)::			iL, iOrder, iOrderhalf, iIndex;

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
	
	call PrintToLog("DerivativeComplex",7)

	iOrder 	= size(arWeights, 1);
	acOut(1:iLength)	= 0;

	if (.not.(iLength == 1)) then
		
		if (iLength - 1 < iOrder) then
			do iL = 1, iLength
				acOut(iL) = &
					sum(acIn(iL+aiPoints(iLength-1, iL,1):iL+aiPoints(iLength-1, iL,2)) * &
					arWeights(iLength-1, iL, 1:aiPoints(iLength-1,iL,3)))/rDelta;
			end do
		else
			iOrderhalf = floor(real(iOrder)/2.0)
			iIndex	= 1;
			
			! The left part
			do iL = 1, iOrderhalf
				acOut(iL) = &
					sum(acIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
				iIndex 	= iIndex + 1;
			end do
			
			! The middle part
			do iL = iOrderhalf + 1, iLength-iOrderhalf
				acOut(iL) = &
					sum(acIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
			end do
			
			iIndex	= iIndex + 1;
			if (modulo(iLength, 2_i8b) == 1) then
				iIndex	= iIndex + 1;
			end if
			
			! The right part
			do iL = iLength-iOrderhalf+1, iLength
				acOut(iL) = &
					sum(acIn(iL+aiPoints(iOrder, iIndex,1):iL+aiPoints(iOrder, iIndex,2)) * &
					arWeights(iOrder, iIndex, 1:aiPoints(iOrder,iIndex,3)))/rDelta;
				iIndex 	= iIndex + 1;
			end do
		end if
	end if	

END SUBROUTINE DerivativeComplex

SUBROUTINE BeamDerivativeZ(cSpace)

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
!   The subroutine BeamDerivativeZ computes the finite difference in the z-
!   dimension for a space containing a (skew) beam.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpace   io   type(Space)   Space structure that contains the beam
!
	type(Space), intent(inout) ::		cSpace
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   cSlice   type(space)   XYZ slice at a single time instant
!   iLt       i8b       Index in t - since each complex variable in Dist.2
!                       actually contains two time instants, the index
!                       is actually running two times too slow
!   iLx       i8b       Index in x
!   iLy       i8b       Index in y
!   j         i8b       Temporary variable
!   iIndex    i8b       Start index for each TXY entry within cSlice
!   iEnd      i8b       End index for each TXY entry within cSlice
!   cDerivLookup  type(DerivLookup)  Structure containing weights and positions
!                                    Of the FD approximation stencils
!   acBuffer1   dpc     Buffer containing the input to DerivativeComplex
!   acBuffer2   dpc     Buffer containing the output from DerivativeComplex
!   iErr        i8b     Error number
! 
	type(space), target :: cSlice;
	integer(i8b)::				iLt, iLx,iLy,j, iIndex, iEnd
	type(DerivLookup)::			cDerivLookup
	complex(dpc), allocatable::				acBuffer1(:),acBuffer2(:);
	integer(i8b):: 	iErr;
	
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
!   DerivLookupInit
!   ReorderDistr0toDistr1
!   MapVt
!   ReorderDistr1toDistr2
!   InitSpace
!   InitGrid
!   GridDistr2CreateEmpty
!   alloctest
!   Distr2ObtainXYZBlock
!   MapVx
!   MapVy
!   DerivativeComplex
!   MapVxInv
!   MapVyInv
!   Distr2PutXYZBlock
!   DestructSpace
!   ReorderDistr2toDistr1
!   MapVtInv
!   ReorderDistr1toDistr0
!   DerivLookupDestroy
!   SwStop
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Initialize the DerivLookup structure
!   Reorder the Grid from Distr. 0 to Distr. 1
!   Map the Grid in the time dimension
!   Reorder the Grid from Distr. 1 to Distr. 2
!   Initialize an XYZ slice space
!   For each time instant duo (i.e. the two time instants stored in one
!     complex variable in Distr. 2)
!       Copy the XYZ data at the current time instant duo from the space to 
!         the XYZ slice
!       Map the XYZ slice in x
!       Map the XYZ slice in y
!       For each index in x
!           For each index in y
!               Extract the variable at the complete z axis for the current
!                 TXY combination from the XYZ slice.
!               Obtain the Finite Difference approximation at this TXY comb.
!               Put the result back to the XYZ slice.
!       Inverse map the XYZ slice in y
!       Inverse map the XYZ slice in x
!       Copy the result in the XYZ slice to the current time instant duo 
!         in the space
!   Destroy the XYZ slice space
!   Reorder the Grid from Distr. 2 to Distr. 1
!   Inverse map the Grid in the time dimension
!   Reorder the Grid from Distr. 1 to Distr. 0
!   Destroy the DerivLookup structure
!
! =============================================================================

	character(LEN = 2048) ::                acTemp;
	
	call PrintToLog("BeamDerivativeZ",2)

!KH in this subroutine, we perform too much work, since the wraparound regions in 
!KH the T dimension (that are implicitly allocated in Dist. 1) are differentiated as well...
!KH the remarks suggest an improved implementation, but the improvement is relatively small

	! Initialize
	call SwStartAndCount(cswDeriv)
	call DerivLookupInit(cModelParams%FDXOrder, iFDMinOrder, cDerivLookup, .false.)

	! First we have to go to distribution 2, whilst performing mappings inbetween
	call ReorderDistr0ToDistr1(cSpace%cGrid);
	call MapVt(cSpace)
!KH ! To reduce computation time, we could employ, as in the function 
!KH ! InhomContrastOperator_Ali:
!KH	call MapVt_Small(cSpace) ! instead of MapVt(cSpace) -- not yet implemented, 
!KH ! but similar as MapVx_small
!KH call SpreadFrequencies(cSpace%cGrid)
    call ReorderDistr1ToDistr2(cSpace%cGrid);

	!Allocate slice
	call PrintToLog("Allocate Slice Space", 3)
	call InitSpace(cSlice, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
						cSpace%cGrid%iProcN-1, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ, &
						0_i8b,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ,0_i8b,0_i8b,0_i8b, &
						cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	!iDimT is such that each processor gets one complex, i.e. two real slices..
    call InitGrid(cSlice, cSpace%cGrid%iProcN, cSpace%cGrid%iProcID, (/ .false., .false., .false., .false./));
	call GridDistr2CreateEmpty(cSlice%cGrid);

	call PrintToLog("Allocated", 3)

    allocate(acBuffer1(1:cSlice%iDimZ),acBuffer2(1:cSlice%iDimZ),STAT=iErr)
    call alloctest(iErr,'Buffers in BeamDerivativeZ')

	! Now we process each slice
!KH	do iLt = 0,(pcGrid%iD2locN+1)/2-1 !-- with computation time reduction trick
	do iLt = 0, cSpace%cGrid%iD2LocN-1
		write (acTemp, '("Start iLt ", I5 ," out of ", I5)') iLt, cSpace%cGrid%iD2LocN-1
		call PrintToLog(acTemp,3);
!        print *, "Punto di Interesse Odierno 6"
		call Distr2ObtainXYZBlock(cSpace%cGrid, cSlice%cGrid, iLt)
!		print *, "Punto di Interesse Odierno 7"
		call MapVx(cSlice)
		call MapVy(cSlice)
		
		do iLx = 0, cSlice%iDimX-1
			do iLy = 0, cSlice%iDimY-1
		 		! Each complex value stores two values of t, so we have to differentiate them both
				! As differentiation is a simple scalar multiplication and addition,
				! this happens automatically (real -> real, imag -> imag)

	            ! Compute the index in acSlice
	            iIndex	= 1 + (iLx * cSlice%cGrid%iD2XS + iLy * cSlice%cGrid%iD2YS);
	            iEnd	= 1 + (iLx * cSlice%cGrid%iD2XS + iLy * cSlice%cGrid%iD2YS +  (cSlice%iDimZ-1) * cSlice%cGrid%iD2ZS);
	            
	            ! Compute the derivative
	            acBuffer1=cSlice%cGrid%pacD2(iIndex:iEnd:cSpace%cGrid%iD2ZS)
	            call DerivativeComplex(acBuffer1, &
				            acBuffer2, &
				            cSlice%iDimZ, &
				            cSlice%dDx, &
				            cDerivLookup%arWeights, &
				            cDerivLookup%aiPoints)
	            cSlice%cGrid%pacD2(iIndex:iEnd:cSpace.cGrid.iD2ZS)	= acBuffer2;

			end do
		end do

		call MapVxInv(cSlice)
		call MapVyInv(cSlice)
		call Distr2PutXYZBlock(cSlice%cGrid, cSpace%cGrid, iLt);

	end do

	!Deallocate buffers and slice
	call PrintToLog("Free Buffers and Slice", 3)
    deallocate(acBuffer1,acBuffer2)
	call DestructSpace(cSlice);

	!Reorder to distribution 0
	call ReorderDistr2ToDistr1(cSpace.cGrid);
!KH    call CollectFrequencies(cSpace%cGrid) !-with computation time reduction trick
!KH    call MapVtInv_Small(cSpace) !-- instead of MapVtInv(cSpace) -- not yet implemented,
!KH    ! but similar as MapVxInv_Small
	call MapVtInv(cSpace)
	call ReorderDistr1ToDistr0(cSpace.cGrid);

	! DeInitialize
	call DerivLookupDestroy(cDerivLookup);
	call SWStop(cswDeriv)
END SUBROUTINE BeamDerivativeZ

END MODULE ParnacDerivative
