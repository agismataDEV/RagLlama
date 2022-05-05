MODULE ParnacConvolution

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
!   The module ParnacConvolution contains subroutines that compute the 
!   convolution of a contrast or primary source space with the Green's function
!   to give a field (correction) in another space. The Green's function may be
!   the function of the free space (i.e. lossless), or a lossy Green's 
!   function.
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

	USE SpecialFun, ONLY : &
        CiSi4, ExpInt_Prime        
	USE ParnacDataMap, ONLY : &
        MapVt, MapVtInv, &
        MapVx, MapVxInv, &
        MapVy, MapVyInv
	USE ParnacDataRedist, ONLY : &
        ReorderDistr0ToDistr1, ReorderDistr1ToDistr0, &
        ReorderDistr1ToDistr2, ReorderDistr2ToDistr1, &
        Distr2ObtainXYZBlock, Distr2PutXYZBlock, &
        Distr2ObtainYMirroredXYZBlock, Distr2AddToXYZBlock
	USE ParnacDataSpaceInit, ONLY : &
        InitSpace, InitGrid, DestructSpace, &
        InitConvolutionSpaces, &
        GridDistr0CreateEmpty, GridDistr0Deallocate, &
                               GridDistr1Deallocate, &
        GridDistr2CreateEmpty, GridDistr2Deallocate, &
        iBeamOffsetT, iBeamOffsetX, iBeamOffsetY
	USE ParnacTransformFilter, ONLY : &
        TransformT, TransformTInv, &
        TransformXY, TransformXYInv, &
        TransformXYZ, TransformXYZInv
            
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
!   ConvolutionGreen     sub   Computes the convolution between a (primary or
!                              contrast) source and the Green's function - sets
!                              up the spaces and calls ComputeConvolution 
!   ComputeConvolution   sub   This subroutine performs the inner loop in the
!                              evaluation of the convolution and is called by 
!                              ConvolutionGreen.
!   GreenXYZBlockFreqT_lossless   sub   Computes the part of the lossless (free
!                                       space) Green's function as requested.
!   GreenXYZBlockFreqT_lossy      sub   Computes the part of the lossy Green's
!                                       function as requested.
!
! =============================================================================
!
	CONTAINS

SUBROUTINE ConvolutionGreen(cSpaceS, cSpaceD, inplacemode)
		
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
!   The subroutine ConvolutionGreen computes the convolution between a primary 
!   or contrast source and the Green's function of the medium (either lossless
!   or lossy). This subroutine sets up the spaces, performs the temporal 
!   Fourier transform and redistributes the source data array to Dist 2, then
!   it calls ComputeConvolution to perform the inner loop of the convolution
!   and then it returns the field space to distribution 0 and wraps up.
!
!   An important assumption that is used in the program and thta is essential 
!   to the data structure and the current implementation: the temporal 
!   dimensions iDimT and angles ThetaT and therefore dTanT are equal for the 
!   source and the destination space. (of course, the start iStartT does not 
!   have to be equal).
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpaceS   io   type(space)  source space in which the primary or contrast
!                               source is located.
!   cSpaceD   io   type(space)  destination space in which the field correction
!                               or estimate will be returned after evaluation.
!   inplacemode   i   lgt       switch whether to employ a memory-efficient in-
!                               place convolution or not -- if true, then this
!                               means that cSpaceS and cSpaceD should point to
!                               the same space!!
!
		type(Space), target, intent(inout)::	cSpaceS
		type(Space), target, intent(inout)::	cSpaceD
		logical(LGT), intent(in) :: inplacemode
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!
!   cSpaceStemp   type(Space)   Temporary source space in the case that the
!                               source space does not have the correct sizes.
!   pcSpaceS      type(Space)   pointer to the Source space
!   pcSpaceD      type(Space)   pointer to the Destination space
!   pcGridS       type(Grid)    pointer to the grid in the Source space
!   pcGridD       type(Grid)    pointer to the grid in the Destination space
!   iLi              i8b        loop counter over all frequencies in Dist 2
!   iLx              i8b        loop counter in the x-dimension
!   iLy              i8b        loop counter in the y-dimension
!   iIndS            i8b        index to positions in the data array of a
!                               Source space
!   iIndS2           i8b        index to positions in the data array of a
!                               Source space
!   acTemp           char       temporary char array for output log messages
!
		type(Space), target ::		cSpaceStemp
		type(Space), pointer::		pcSpaceS, pcSpaceD
		type(Grid), pointer::		pcGridS, pcGridD
		integer						iErr
		integer(i8b) ::				iLi, iLx, iLy, iIndS, iIndS2
		character(len=1024)::		acTemp;
!		logical(lgt) ::				lStudy
	
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
!   SwStartAndCount
!   ReorderDistr0ToDistr1
!   MapVt
!   TransformT
!   ReorderDistr1ToDistr2
!   InitSpace
!   InitGrid
!   GridDistr2Deallocate
!   GridDistr0DeAllocate
!   GridDistr0DeAllocate
!   GridDistr2CreateEmpty
!   ComputeConvolution
!   DestructSpace
!   ReorderDistr2ToDistr1
!   TransformTInv
!   ReorderDistr1ToDistr0
!   SWstop
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Perform a number of checks on the source and dest spaces
!   Redistribute the source space from Distribution 0 to Distribution 1
!   Map the source space in the time dimension
!   Transform the source space in the time dimension
!   Redistribute the source space from Distr. 1 to Distr. 2
!   If the source space is a primary source or a field slice, 
!     and the dimension in z > 1
!       create a temporary source space with one layer in z
!       copy the source definition or field slice to this temporary source 
!         space
!       make sure that the convolution is performed on the temporary source
!         space instead of the original source space
!   Compute the convolution between the source space and the Green's function
!   If necessary, destruct the temporary source space
!   Redistribute the field space from Distr. 2 to Distr. 1
!   Inverse transform the field space in the time dimension
!   Inverse map the field space in the time dimension
!   Redistribute the field space from Distr. 1 to Distr. 0
!
! =============================================================================
		
		call PrintToLog("Starting ConvolutionGreen", 2)
		call SwStartAndCount(cswLinConv);

		pcSpaceS => cSpaceS
		pcSpaceD => cSpaceD
		pcGridS	=> cSpaceS%cGrid
		!in case of inplace, the source grid will be 
		!used to store the destination grid values
		if(inplacemode) then
			call PrintToLog("inplacemode mode",1)
			pcGridD	=> cSpaceS%cGrid
		else
			pcGridD	=> cSpaceD%cGrid
		end if

!!		if(pcSpaceS%iSpaceIdentifier==iSI_CONTRASTSOURCE.and.pcSpaceS%iStartZ==1.and.pcSpaceD%iStartZ>1) then
!		if(pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
!			lStudy=.true.
!		else
!			lStudy=.false.
!		end if

		! Check source space whether it is a source space
		if(pcSpaceS%iSpaceIdentifier/=iSI_PRIMARYSOURCE .and. &
			pcSpaceS%iSpaceIdentifier/=iSI_CONTRASTSOURCE .and. &
			pcSpaceS%iSpaceIdentifier/=iSI_FIELDSLICE ) then
!			print *, pcSpaceS%iSpaceIdentifier
			write (acTemp, '("Error occurred during the execution of ConvolutionGreen, aborting program")');
			call PrintToLog(acTemp, -1);
			write (acTemp, '("the space provided as input is not a PRIMARY or CONTRAST SOURCE type.")');
			call PrintToLog(acTemp, -1);
			stop
		end if

		! Check whether the time dimensions and angle are the same
		if((pcSpaceD%iDimT /= pcSpaceS%iDimT) .or. (pcSpaceD%dTanT/=pcSpaceS%dTanT)) then
			write (acTemp, '("Error occurred during the execution of ConvolutionGreen, aborting program")');
			call PrintToLog(acTemp, -1);
			write (acTemp, '("Source and Dest spaces do not have the same time dimension and/or angle")');
			call PrintToLog(acTemp, -1);
			stop
		end if
		
!if(lStudy) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSpaceS0",&
!					"cSpace", &
!					pcSpaceS, &
!					(/ 0_i8b, 0_i8b, 0_i8b, pcSpaceS%iDimZ-2 /), &
!					(/ pcSpaceS%iDimT, pcSpaceS%iDimX, pcSpaceS%iDimY, 2_i8b /), &
!					.true.);
!end if

		! 1a) Put the data from in cGridS into distribution 2
		if (pcGridS%iDistr==0) then
			call PrintToLog("Convert the source from distribution 0 to 1", 3)
			call ReorderDistr0ToDistr1(pcGridS);
			call MapVt(pcSpaceS)
			call PrintToLog("Compute the Fourier Transform with respect to the T-axis", 3)
			!print *, size( pcGridS%pacD1(:))
			!print *, "Punto di Interesse 1"
			call TransformT(pcGridS)
			!print *, size( pcGridS%pacD1(:))
			!print *, "Punto di Interesse 2"
		end if
		if (pcGridS%iDistr==1) then
			call PrintToLog("Convert the source from distribution 1 to 2", 3)
			call ReorderDistr1ToDistr2(pcGridS);
		end if

		! 1b) If the source is a primary source or a field slice, test whether we have only 1 layer in z.
		!		This may not be so since in the T-local distributions 0 and 1 we need
		!		a nice number of points in XYZ.
		!		If not, make another space with only one layer in z and copy the source layer 
        !       or field slice into it.
		if(((pcSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE) .or. &
            (pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE)) .and. &
           (pcSpaceS%iDimZ>1)) then

            if (pcSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE) then
			!Create temporary space containing the source at 
			! z= 0
			    call InitSpace(cSpaceStemp, iSI_PRIMARYSOURCE, .false., &
						    pcSpaceS%iDimT,pcSpaceS%iDimX, pcSpaceS%iDimY, 1_i8b, &
						    pcSpaceS%iStartT,0_i8b,0_i8b,0_i8b,0_i8b,0_i8b,0_i8b, &
						    pcSpaceS%dFnyq, pcSpaceS%dTanX, pcSpaceS%dTanY, pcSpaceS%dTanT)
            elseif (pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
			!Create temporary space containing the slice at 
			! z= iStartZ
			    call InitSpace(cSpaceStemp, iSI_FIELDSLICE, pcSpaceS%bYSymm, &
						    pcSpaceS%iDimT,pcSpaceS%iDimX, pcSpaceS%iDimY, 1_i8b, &
						    pcSpaceS%iStartT,0_i8b,0_i8b,pcSpaceS%iStartZ,0_i8b,0_i8b,0_i8b, &
						    pcSpaceS%dFnyq, pcSpaceS%dTanX, pcSpaceS%dTanY, pcSpaceS%dTanT)
            end if

 	        call InitGrid(cSpaceStemp, pcGridS%iProcN, pcGridS%iProcID, (/ .false., .false., .false., .false./));

            ! Some tests that should not be necessary but are here because of a strange error
            ! on the Huygens computer
			if(associated(cSpaceStemp%cGrid%parD0)) then
				call PrintToLog("******** WARNING: Space parD0 is already associated",-1)
				nullify(cSpaceStemp%cGrid%parD0)
			end if
			if(associated(cSpaceStemp%cGrid%pacD1)) then
				call PrintToLog("******** WARNING: Space pacD1 is already associated",-1)
				nullify(cSpaceStemp%cGrid%pacD1)
			end if
			if(associated(cSpaceStemp%cGrid%pacD2)) then
				call PrintToLog("******** WARNING: Space pacD2 is already associated",-1)
				nullify(cSpaceStemp%cGrid%pacD2)
			end if

            ! Copy the source/field layer into the temporary source space
			call GridDistr2CreateEmpty(cSpaceStemp%cGrid);
			do iLi=0,pcGridS%iD2locN-1
				do iLx=0,pcGridS%iD2XL-1
					do iLy=0,pcGridS%iD2YL-1
						iIndS = 1 + iLi*pcGridS%iD2IS + iLx*pcGridS%iD2XS + iLy*pcGridS%iD2YS
						iIndS2 = 1 + iLi*cSpaceStemp%cGrid%iD2IS + &
									 iLx*cSpaceStemp%cGrid%iD2XS + iLy*cSpaceStemp%cGrid%iD2YS
						cSpaceStemp%cGrid%pacD2(iIndS2)=pcGridS%pacD2(iIndS)
					end do
				end do
			end do

            ! Associate the source space and grid pointers to the temp source space
			pcSpaceS => cSpaceStemp
			pcGridS	=> cSpaceStemp%cGrid

			!!deallocate original source space!!
			call GridDistr2Deallocate(cSpaceS%cGrid)			
		end if

		! 2) If the grid in Dest is not already initialized 
		!	  (which would already be so if e.g. it is equal to the Source space),
		!	  initialize Grid in Dest
		if (pcGridD%iDistr/=2) then
			if(pcGridD%iDistr==0) then
				call GridDistr0DeAllocate(pcGridD)
			else if(pcGridD%iDistr==1) then
				call GridDistr1DeAllocate(pcGridD)
			end if
			call GridDistr2CreateEmpty(pcGridD)
		end if
        
		! 3) Compute the convolution with the Green's function
		!		In this subroutine, the real work is performed
		call ComputeConvolution(pcSpaceS,pcGridS,pcSpaceD,pcGridD)

		! If there was a temporary source space, destruct it now
		if((cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .and. cSpaceS%iDimZ>1).or. &
			(cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE .and. cSpaceS%iDimZ>1)) then
			call DestructSpace(cSpaceStemp)
			nullify(pcSpaceS,pcGridS)
			call PrintToLog("Temp source space destructed", 3)
		end if

		!in case of inplacemode, the source grid has been
		!filled with the destination grid values
		!now reference destination grid to source grid
		if(inplacemode) then
			pcSpaceD%cGrid%pacD2=>pcSpaceS%cGrid%pacD2
			pcSpaceD%cGrid%iDistr=2
			pcGridD	=> cSpaceD%cGrid
		end if
		
		! 6) Redistribute cGridD from distribution 2 to 0
		call PrintToLog("Redistribute the result to distribution 1", 3)
		call ReorderDistr2ToDistr1(pcGridD);
		call PrintToLog("Compute the inverse Fourier Transform with respect to the T-axis", 3)
		call TransformTInv(pcGridD)
		call MapVtInv(pcSpaceD)
		call PrintToLog("Redistribute the result to distribution 0", 3)
		call ReorderDistr1ToDistr0(pcGridD);
		
		call PrintToLog("Finished ConvolutionGreen", 2)
		call SWStop(cswLinConv);

END SUBROUTINE ConvolutionGreen

SUBROUTINE ConvolutionGreenConj(cSpaceS, cSpaceD, inplacemode)
		
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
!   The subroutine ConvolutionGreen computes the convolution between a primary 
!   or contrast source and the Green's function of the medium (either lossless
!   or lossy). This subroutine sets up the spaces, performs the temporal 
!   Fourier transform and redistributes the source data array to Dist 2, then
!   it calls ComputeConvolution to perform the inner loop of the convolution
!   and then it returns the field space to distribution 0 and wraps up.
!
!   An important assumption that is used in the program and thta is essential 
!   to the data structure and the current implementation: the temporal 
!   dimensions iDimT and angles ThetaT and therefore dTanT are equal for the 
!   source and the destination space. (of course, the start iStartT does not 
!   have to be equal).
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpaceS   io   type(space)  source space in which the primary or contrast
!                               source is located.
!   cSpaceD   io   type(space)  destination space in which the field correction
!                               or estimate will be returned after evaluation.
!   inplacemode   i   lgt       switch whether to employ a memory-efficient in-
!                               place convolution or not -- if true, then this
!                               means that cSpaceS and cSpaceD should point to
!                               the same space!!
!
		type(Space), target, intent(inout)::	cSpaceS
		type(Space), target, intent(inout)::	cSpaceD
		logical(LGT), intent(in) :: inplacemode
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!
!   cSpaceStemp   type(Space)   Temporary source space in the case that the
!                               source space does not have the correct sizes.
!   pcSpaceS      type(Space)   pointer to the Source space
!   pcSpaceD      type(Space)   pointer to the Destination space
!   pcGridS       type(Grid)    pointer to the grid in the Source space
!   pcGridD       type(Grid)    pointer to the grid in the Destination space
!   iLi              i8b        loop counter over all frequencies in Dist 2
!   iLx              i8b        loop counter in the x-dimension
!   iLy              i8b        loop counter in the y-dimension
!   iIndS            i8b        index to positions in the data array of a
!                               Source space
!   iIndS2           i8b        index to positions in the data array of a
!                               Source space
!   acTemp           char       temporary char array for output log messages
!
		type(Space), target ::		cSpaceStemp
		type(Space), pointer::		pcSpaceS, pcSpaceD
		type(Grid), pointer::		pcGridS, pcGridD
		integer(i8b) ::				iLi, iLx, iLy, iIndS, iIndS2
		character(len=1024)::		acTemp;
!		logical(lgt) ::				lStudy
	
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
!   SwStartAndCount
!   ReorderDistr0ToDistr1
!   MapVt
!   TransformT
!   ReorderDistr1ToDistr2
!   InitSpace
!   InitGrid
!   GridDistr2Deallocate
!   GridDistr0DeAllocate
!   GridDistr0DeAllocate
!   GridDistr2CreateEmpty
!   ComputeConvolution
!   DestructSpace
!   ReorderDistr2ToDistr1
!   TransformTInv
!   ReorderDistr1ToDistr0
!   SWstop
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Perform a number of checks on the source and dest spaces
!   Redistribute the source space from Distribution 0 to Distribution 1
!   Map the source space in the time dimension
!   Transform the source space in the time dimension
!   Redistribute the source space from Distr. 1 to Distr. 2
!   If the source space is a primary source or a field slice, 
!     and the dimension in z > 1
!       create a temporary source space with one layer in z
!       copy the source definition or field slice to this temporary source 
!         space
!       make sure that the convolution is performed on the temporary source
!         space instead of the original source space
!   Compute the convolution between the source space and the Green's function
!   If necessary, destruct the temporary source space
!   Redistribute the field space from Distr. 2 to Distr. 1
!   Inverse transform the field space in the time dimension
!   Inverse map the field space in the time dimension
!   Redistribute the field space from Distr. 1 to Distr. 0
!
! =============================================================================
		
		call PrintToLog("Starting ConvolutionGreen", 2)
		call SwStartAndCount(cswLinConv);

		pcSpaceS => cSpaceS
		pcSpaceD => cSpaceD
		pcGridS	=> cSpaceS%cGrid
		!in case of inplace, the source grid will be 
		!used to store the destination grid values
		if(inplacemode) then
			call PrintToLog("inplacemode mode",1)
			pcGridD	=> cSpaceS%cGrid
		else
			pcGridD	=> cSpaceD%cGrid
		end if

!!		if(pcSpaceS%iSpaceIdentifier==iSI_CONTRASTSOURCE.and.pcSpaceS%iStartZ==1.and.pcSpaceD%iStartZ>1) then
!		if(pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
!			lStudy=.true.
!		else
!			lStudy=.false.
!		end if

		! Check source space whether it is a source space
		if(pcSpaceS%iSpaceIdentifier/=iSI_PRIMARYSOURCE .and. &
			pcSpaceS%iSpaceIdentifier/=iSI_CONTRASTSOURCE .and. &
			pcSpaceS%iSpaceIdentifier/=iSI_FIELDSLICE ) then
!			print *, pcSpaceS%iSpaceIdentifier
			write (acTemp, '("Error occurred during the execution of ConvolutionGreen, aborting program")');
			call PrintToLog(acTemp, -1);
			write (acTemp, '("the space provided as input is not a PRIMARY or CONTRAST SOURCE type.")');
			call PrintToLog(acTemp, -1);
			stop
		end if

		! Check whether the time dimensions and angle are the same
		if((pcSpaceD%iDimT /= pcSpaceS%iDimT) .or. (pcSpaceD%dTanT/=pcSpaceS%dTanT)) then
			write (acTemp, '("Error occurred during the execution of ConvolutionGreen, aborting program")');
			call PrintToLog(acTemp, -1);
			write (acTemp, '("Source and Dest spaces do not have the same time dimension and/or angle")');
			call PrintToLog(acTemp, -1);
			stop
		end if
		
!if(lStudy) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSpaceS0",&
!					"cSpace", &
!					pcSpaceS, &
!					(/ 0_i8b, 0_i8b, 0_i8b, pcSpaceS%iDimZ-2 /), &
!					(/ pcSpaceS%iDimT, pcSpaceS%iDimX, pcSpaceS%iDimY, 2_i8b /), &
!					.true.);
!end if

		! 1a) Put the data from in cGridS into distribution 2
		if (pcGridS%iDistr==0) then
			call PrintToLog("Convert the source from distribution 0 to 1", 3)
			call ReorderDistr0ToDistr1(pcGridS);
			call MapVt(pcSpaceS)
			call PrintToLog("Compute the Fourier Transform with respect to the T-axis", 3)
			!print *, size( pcGridS%pacD1(:))
			!print *, "Punto di Interesse 1"
			call TransformT(pcGridS)
			!print *, size( pcGridS%pacD1(:))
			!print *, "Punto di Interesse 2"
		end if
		if (pcGridS%iDistr==1) then
			call PrintToLog("Convert the source from distribution 1 to 2", 3)
			call ReorderDistr1ToDistr2(pcGridS);
		end if

		! 1b) If the source is a primary source or a field slice, test whether we have only 1 layer in z.
		!		This may not be so since in the T-local distributions 0 and 1 we need
		!			 a nice number of points in XYZ.
		!		If not, make another space with only one layer in z and copy the source layer 
        !       or field slice into it.
		if(((pcSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE) .or. &
            (pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE)) .and. &
           (pcSpaceS%iDimZ>1)) then

            if (pcSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE) then
			!Create temporary space containing the source at 
			! z= 0
			    call InitSpace(cSpaceStemp, iSI_PRIMARYSOURCE, .false., &
						    pcSpaceS%iDimT,pcSpaceS%iDimX, pcSpaceS%iDimY, 1_i8b, &
						    pcSpaceS%iStartT,0_i8b,0_i8b,0_i8b,0_i8b,0_i8b,0_i8b, &
						    pcSpaceS%dFnyq, pcSpaceS%dTanX, pcSpaceS%dTanY, pcSpaceS%dTanT)
            elseif (pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
			!Create temporary space containing the slice at 
			! z= iStartZ
			    call InitSpace(cSpaceStemp, iSI_FIELDSLICE, pcSpaceS%bYSymm, &
						    pcSpaceS%iDimT,pcSpaceS%iDimX, pcSpaceS%iDimY, 1_i8b, &
						    pcSpaceS%iStartT,0_i8b,0_i8b,pcSpaceS%iStartZ,0_i8b,0_i8b,0_i8b, &
						    pcSpaceS%dFnyq, pcSpaceS%dTanX, pcSpaceS%dTanY, pcSpaceS%dTanT)
            end if

 	        call InitGrid(cSpaceStemp, pcGridS%iProcN, pcGridS%iProcID, (/ .false., .false., .false., .false./));

            ! Some tests that should not be necessary but are here because of a strange error
            ! on the Huygens computer
			if(associated(cSpaceStemp%cGrid%parD0)) then
				call PrintToLog("******** WARNING: Space parD0 is already associated",-1)
				nullify(cSpaceStemp%cGrid%parD0)
			end if
			if(associated(cSpaceStemp%cGrid%pacD1)) then
				call PrintToLog("******** WARNING: Space pacD1 is already associated",-1)
				nullify(cSpaceStemp%cGrid%pacD1)
			end if
			if(associated(cSpaceStemp%cGrid%pacD2)) then
				call PrintToLog("******** WARNING: Space pacD2 is already associated",-1)
				nullify(cSpaceStemp%cGrid%pacD2)
			end if

            ! Copy the source/field layer into the temporary source space
			call GridDistr2CreateEmpty(cSpaceStemp%cGrid);
			do iLi=0,pcGridS%iD2locN-1
				do iLx=0,pcGridS%iD2XL-1
					do iLy=0,pcGridS%iD2YL-1
						iIndS = 1 + iLi*pcGridS%iD2IS + iLx*pcGridS%iD2XS + iLy*pcGridS%iD2YS
						iIndS2 = 1 + iLi*cSpaceStemp%cGrid%iD2IS + &
									 iLx*cSpaceStemp%cGrid%iD2XS + iLy*cSpaceStemp%cGrid%iD2YS
						cSpaceStemp%cGrid%pacD2(iIndS2)=pcGridS%pacD2(iIndS)
					end do
				end do
			end do

            ! Associate the source space and grid pointers to the temp source space
			pcSpaceS => cSpaceStemp
			pcGridS	=> cSpaceStemp%cGrid

			!!deallocate original source space!!
			call GridDistr2Deallocate(cSpaceS%cGrid)			
		end if

		! 2) If the grid in Dest is not already initialized 
		!	  (which would already be so if e.g. it is equal to the Source space),
		!	  initialize Grid in Dest
		if (pcGridD%iDistr/=2) then
			if(pcGridD%iDistr==0) then
				call GridDistr0DeAllocate(pcGridD)
			else if(pcGridD%iDistr==1) then
				call GridDistr1DeAllocate(pcGridD)
			end if
			call GridDistr2CreateEmpty(pcGridD)
		end if
        
		! 3) Compute the convolution with the Green's function
		!		In this subroutine, the real work is performed
!		print *, "Point of Interest"
		call ComputeConvolutionConj(pcSpaceS,pcGridS,pcSpaceD,pcGridD)

		! If there was a temporary source space, destruct it now
		if((cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .and. cSpaceS%iDimZ>1).or. &
			(cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE .and. cSpaceS%iDimZ>1)) then
			call DestructSpace(cSpaceStemp)
			nullify(pcSpaceS,pcGridS)
			call PrintToLog("Temp source space destructed", 3)
		end if

		!in case of inplacemode, the source grid has been
		!filled with the destination grid values
		!now reference destination grid to source grid
		if(inplacemode) then
			pcSpaceD%cGrid%pacD2=>pcSpaceS%cGrid%pacD2
			pcSpaceD%cGrid%iDistr=2
			pcGridD	=> cSpaceD%cGrid
		end if
		
		! 6) Redistribute cGridD from distribution 2 to 0
		call PrintToLog("Redistribute the result to distribution 1", 3)
		call ReorderDistr2ToDistr1(pcGridD);
		call PrintToLog("Compute the inverse Fourier Transform with respect to the T-axis", 3)
		call TransformTInv(pcGridD)
		call MapVtInv(pcSpaceD)
		call PrintToLog("Redistribute the result to distribution 0", 3)
		call ReorderDistr1ToDistr0(pcGridD);
		
		call PrintToLog("Finished ConvolutionGreen", 2)
		call SWStop(cswLinConv);

END SUBROUTINE ConvolutionGreenConj

SUBROUTINE ComputeConvolution(cSpaceS, cGridS, cSpaceD, cGridD)

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
!   The subroutine ComputeConvolution performs the inner loop in the 
!   convolution of a primary or contrast source space with the lossless or 
!   lossy Green's function of the medium. The Grids are assumed to be in Distr.
!   2. Each frequency is separately addressed. The routine is capable of
!   performing an inplace evaluation (if cGridS and cGridD point to the same
!   data array), and it has special memory-efficient implementations for 
!   a source plane and for symmetry in the y-dimension.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpaceS   i   type(Space)   Source space - used only for the parameters
!   cGridS    i   type(Grid)    Source grid  - contains the source data array
!   cSpaceD   i   type(Space)   Destination space -used only for the parameters
!   cGridD    io  type(Grid)    Destination grid - upon exiting the subroutine,
!                               contains the field data array 
!
		type(Space), intent(in) :: cSpaceS
		type(Grid), intent(in) :: cGridS
		type(Space), intent(in) :: cSpaceD
		type(Grid), intent(inout) :: cGridD 
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   cSliceS   type(Space)   XYZ-slice of a single frequency of the source, 
!                           including wraparound regions in all dimensions
!   cSliceD   type(Space)   XYZ-slice of a single frequency of the field or the
!                           Green's function, including wraparound regions in 
!                           all dimensions
!   iDiffT        i8b   Difference in time steps between the destination space
!                       and the source space
!   iDimT         i8b   Time dimension of the destination space and the source
!                       space
!   iNumYBlocks   i8b   Number of blocks in the y-dimension: 2 if y-symmetry is
!                       employed, 1 if it is not.
!   iYBlocks      i8b   Loop counter over all blocks in the y-dimension
!   iOmega        i8b   Loop counter over all frequencies on this processor
!   iLx           i8b   Loop counter over all indices in the x-dimension
!   iLy           i8b   Loop counter over all indices in the y-dimension
!   iLz           i8b   Loop counter over all indices in the z-dimension
!   iIndS         i8b   Index into the data array of the source space
!   iIndG         i8b   Index into the data array of the destination space
!   iIndSUpdate   i8b   Increase of iIndS when looping over all z-positions
!                       Normally, this is equal to iD2ZS, but for the primary
!                       source plane and the field plane, iIndSUpdate = 0
!   dOmega         dp   Step size in the angular frequency axis
!   dCorrection    dp   Correction factor for the forward FFT's
!   dKcutoff       dp   Spatial angular cutoff frequency
!   acTemp       char   Temporary char array for output log messages
!
		type(Space) ::		cSliceS, cSliceD
		integer(i8b) ::		iDiffT, iDimT, iNumYBlocks, iYBlocks, &
				 iOmega, iLx, iLy, iLz, iIndS, iIndG, iIndSUpdate
		real(dp) ::			dOmega, dCorrection, dKcutoff
		character(len=1024)::		acTemp;
!		logical(lgt) ::				lStudy
	
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
!   InitConvolutionSpaces
!   InitGrid
!   GridDistr2CreateEmpty
!   Distr2ObtainXYZBlock
!   MapVx
!   MapVy
!   GreenXYZBlockFreqT_lossless
!   GreenXYZBlockFreqT_lossy
!   TransformXY
!   TransformXYZ
!   TransformXYInv
!   TransformXYZInv
!   MapVxInv
!   MapVyInv
!   Distr2PutXYZBlock
!   Distr2ObtainYMirroredXYZBlock
!   Distr2AddToXYZBlock
!   DestructSpace
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Initialize source and destination XYZ slices that include the wraparound
!     regions in all dimensions
!   If symmetry in Y is employed, set the number of blocks in Y to 2, otherwise
!     set to 1
!   For each frequency on this processor do
!       Copy the XYZ data at the current frequency from the source space to 
!         the source XYZ slice
!       For all blocks in Y do
!           Map the source XYZ slice in x and y
!           Fill the destination XYZ slice with the lossless or the lossy 
!             Green's function
!           Map the destination XYZ slice in x and y
!           If the source space contains a primary source or a field plane
!               Transform the source and destination XYZ slices in x and y
!               Apply Fourier correction factor for each Forward transform
!           Else the source space contains a volume contrast source
!               Transform the source and destination XYZ slices in x, y and z
!               Apply Fourier correction factor for each Forward transform
!           Multiply the source and destination XYZ slices elementwise, and
!             put the result into the destination XYZ slice
!           If the source space contains a primary source or a field plane
!               Inverse transform the destination XYZ slice in x and y
!           Else the source space contains a volume contrast source
!               Inverse transform the destination XYZ slice in x, y and z
!           Inverse map the destination XYZ slice in x and y
!           If symmetry in Y is not employed
!               Copy the destination XYZ slice to the current frequency in the
!                 destination space
!           Else
!               If it is the first Y-block
!                   First obtain the XYZ data at the current frequency from the
!                     source space to the source XYZ slice in an order mirrored 
!                     in y and displaced in y to the negative y-axis
!                   Then copy the destination XYZ slice to the current 
!                     frequency in the destination space (order of these two 
!                     operations is important for inplace evaluation)
!               If it is the second Y-block
!                   Add the destination XYZ slice to the current frequency in
!                   the destination space
!   De-initialize source and destination XYZ slices           
!
! =============================================================================

		call PrintToLog("ComputeConvolution", 3)

!!		if(pcSpaceS%iSpaceIdentifier==iSI_CONTRASTSOURCE.and.pcSpaceS%iStartZ==1.and.pcSpaceD%iStartZ>1) then
!		if(pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
!			lStudy=.true.
!		else
!			lStudy=.false.
!		end if

		! 3) Initialize convolution slice spaces
		call PrintToLog("Allocate Slice spaces S and D", 3)
		call InitConvolutionSpaces(cSliceS, cSliceD, cSpaceS, cSpaceD)
	    call InitGrid(cSliceS, cGridS%iProcN, cGridS%iProcID, (/ .false., .false., .false., .false./));
	    call InitGrid(cSliceD, cGridD%iProcN, cGridD%iProcID, (/ .false., .false., .false., .false./));
		call GridDistr2CreateEmpty(cSliceS%cGrid);
		call GridDistr2CreateEmpty(cSliceD%cGrid);

		! 4) Now we perform the convolution on each xyz-block we have stored locally and store the result
		call PrintToLog("Start the inner loop of the convolution", 3)
		! For the T axis we assume that pcSpaceD%iDimT = pcSpaceS%iDimT and pcSpaceD%dTanT=pcSpaceS%dTanT
		! These assumptions are already important for the mapping and transform in T in dist. 1..
		iDiffT = cSpaceD%iStartT - cSpaceS%iStartT
		
		iDimT = cSpaceS%iDimT
		dKcutoff = two_pi * ( 1.0_dp + 0.5_dp/cSpaceS%iDimX ) * cSpaceS%dFnyq
		
		!Y symmetry: if true, then we calculate two blocks of a contrast source
		!Not for ISI_PRIMARYSOURCE: not necessary because when evaluating the primary-source-to-field
        ! convolution we are not close to maximum memory anyway
		if ((cSpaceS%iSpaceIdentifier==ISI_CONTRASTSOURCE.or.cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE)&
			.and.(cSpaceS%bYSymm .EQV. .true.)) then
			iNumYBlocks=2
		else
			iNumYBlocks=1
		end if

!if(lStudy) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSpaceS1",&
!					"cSpace", &
!					cSpaceS, &
!					(/ 0_i8b, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ cSpaceS%iDimT+1, cSpaceS%iDimX, cSpaceS%iDimY, cSpaceS%iDimZ /), &
!					.false.);
!end if
		!Main loop in the convolution: for each omega, calculate the convolution in the 
		! spectral domain
		
		do iOmega = 0, cGridS%iD2LocN-1
			write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iOmega, cGridS.iD2LocN-1
			call PrintToLog(acTemp,3);

			dOmega = two_pi * cGridS%aiD2Loc(1+iOmega) * cSpaceS%dFnyq / cSpaceS%iDimT
			
			!KH change for periodical T to 
			!KH dOmega = two_pi * cGridS%aiD2Loc(1+iOmega) * 2.0_dp*cSpaceS%dFnyq/cSpaceS%iDimT

			call PrintToLog("Obtain Source slice", 4)
			
			call Distr2ObtainXYZBlock(cGridS, cSliceS%cGrid, iOmega)   !!! L.D. 11-11-2011
			do iYblocks=1,iNumYBlocks !stays 1 if bYSymm==.false.

				!Obtain and map Source slice
				call PrintToLog("Map Source slice", 4)
				call MapVx(cSliceS)
				call MapVy(cSliceS)

				!Obtain and map Green's function slice
				if(cMediumParams%a_att==0.0_dp) then
					call PrintToLog("Obtain and map lossless Green's function slice", 4)
	 				call GreenXYZBlockFreqT_lossless(cSliceD, cSliceS, dOmega, dKcutoff, iDiffT, iDimT)
				else
					call PrintToLog("Obtain and map lossy Green's function slice", 4)
	 				call GreenXYZBlockFreqT_lossy(cSliceD, cSliceS, dOmega, dKcutoff, iDiffT, iDimT)
				end if

!if(lStudy.and.pcGridD%aiD2Loc(1+iOmega)==36.and.iYblocks==1) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSliceD0",&
!					"cSlice", &
!					cSliceD, &
!					(/ cSliceD%cGrid%iProcID, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ 1_i8b, cSliceD%iDimX, cSliceD%iDimY, cSliceD%iDimZ /), &
!					.false.);
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSliceS0",&
!					"cSlice", &
!					cSliceS, &
!					(/ cSliceD%cGrid%iProcID, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ 1_i8b, cSliceS%iDimX, cSliceS%iDimY, cSliceS%iDimZ /), &
!					.false.);
!
!end if

				call MapVx(cSliceD);
				call MapVy(cSliceD);

				if(cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .or. cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
				!iDimZ=1, don't transform in z

					call PrintToLog("Transform Source and Green's function slice in XY", 4)
					call TransformXY(cSliceS%cGrid)
					call TransformXY(cSliceD%cGrid)

					!All factors for the forward transforms of G and S and the inverse transform combined 
					!dCorrection	= cSliceD%dDx**4 * cSliceD%dDt * cSliceD%dDf * cSliceD%dDkx * cSliceD%dDky;
                    ! or equivalently:
					dCorrection	= cSliceD%dDx**2/( real(2*cSpaceS%iDimT * cSliceD%iDimX,dp) * cSliceD%iDimY); 
					!the real() command is to prevent integer overflow error from occurring
					iIndSupdate = 0

				else
				!iSI_CONTRASTSOURCE, iDimZ can be anything...

					call PrintToLog("Transform Source and Green's function slice in XYZ", 4)
!					print *, "Punto di Interesse Odierno 3 Bis New"
					call TransformXYZ(cSliceS%cGrid,.false.)
!					print *, "Punto di Interesse Odierno 4"
					call TransformXYZ(cSliceD%cGrid,.true.)

					!All factors for the forward transforms of G and S and the inverse transform combined 
					!dCorrection	= cSliceD%dDx**6 * cSliceD%dDt * cSliceD%dDf * cSliceD%dDkx * cSliceD%dDky * cSliceD%dDkz;
                    ! or equivalently:
					dCorrection	= cSliceD%dDx**3/( real(2*cSpaceS%iDimT * cSliceD%iDimX,dp) * cSliceD%iDimY * cSliceD%iDimZ); !should give the same result
					!the real command is to prevent integer overflow error from occurring
					iIndSupdate = cSliceS%cgrid%iD2ZS
		
				end if

				! Take the inproduct of cGridS and of cGridG. Put them in distribution 2 of cGridG.
				! since we also have to multiply the result of the Fourier transforms with a constant factor, 
				! we'll do  it all in one go and multiply each value with the correction factor dCorrection. Here 
				call PrintToLog("Compute the inproduct for the two xyz-blocks", 4)
				if (cSpaceS%iSpaceIdentifier==iSI_CONTRASTSOURCE) then
					! for contrastsource only..:
					cSliceD%cGrid%pacD2	= cSliceD%cGrid%pacD2 * cSliceS%cGrid%pacD2 * dCorrection;
				else
					do iLx = 0, cSliceD%cgrid%iD2XL-1
						do iLy = 0, cSliceD%cgrid%iD2YL-1
							iIndS	=	1 + iLx * cSliceD%cgrid%iD2XS + iLy * cSliceD%cgrid%iD2YS;
							iIndG	=	1 + iLx * cSliceD%cgrid%iD2XS + iLy * cSliceD%cgrid%iD2YS;
							do iLz = 0, cSliceD%cgrid%iD2ZL-1
								cSliceD%cGrid%pacD2(iIndG)	= cSliceD%cGrid%pacD2(iIndG) * &
															  cSliceS%cGrid%pacD2(iIndS) * dCorrection;
								iIndS	= iIndS + iIndSupdate
								iIndG	= iIndG + cSliceD%cgrid%iD2ZS;
								!if iSI_PRIMARYSOURCE, iIndsSupdate=0 and we therefore stick to z=0...
							end do
						end do
					end do
				endif

				! 5) Take the inverse transform of cGridP with respect to the XY axes
				if(cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .or. cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
				!iDimZ=1, don't transform in z

					call PrintToLog("Inverse transform Field slice in XY", 4)
					call TransformXYInv(cSliceD%cGrid)

				else

					call PrintToLog("Inverse transform Field slice in XYZ", 4)
					call TransformXYZInv(cSliceD%cGrid, .false.)

				end if

 				call MapVxInv(cSliceD)
 				call MapVyInv(cSliceD)

				!Storing of the block: depends on whether we have Y symmetry...
				if (cSpaceS%bYSymm .EQV. .false.) then

					! Put the slice into the destination grid
					call Distr2PutXYZBlock(cSliceD%cGrid, cGridD, iOmega)

				else

					if (iYBlocks==1) then

						!Y symmetry: fill cSliceS with a mirrored version of pcSpaceS
						! BEFORE pcSpaceD is written (it could be the same space...)
						call PrintToLog("Obtain Mirrored Source slice", 4)
						call Distr2ObtainYMirroredXYZBlock(cGridS, cSliceS%cGrid,  iOmega)
						! set iStartY parameter of cSliceS to get the negative starting point
						cSliceS%iStartY=-(cSpaceS%iDimY-1)
						! Put the slice into the destination grid
						call Distr2PutXYZBlock(cSliceD%cGrid, cGridD, iOmega)

					else

						!add result to cSpaceD
						call Distr2AddToXYZBlock(cSliceD%cGrid, cGridD, iOmega)
						! reset iStartY parameter of cSliceS to normal
						cSliceS%iStartY=0

					end if

				end if

			end do 
			
		end do

		! Free the slices and if necessary the temporary and original source space 
		call DestructSpace(cSliceS)
		call DestructSpace(cSliceD)
		call PrintToLog("Slice spaces S and D destructed", 3)

END SUBROUTINE ComputeConvolution

SUBROUTINE ComputeConvolutionConj(cSpaceS, cGridS, cSpaceD, cGridD)

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
!   The subroutine ComputeConvolution performs the inner loop in the 
!   convolution of a primary or contrast source space with the lossless or 
!   lossy Green's function of the medium. The Grids are assumed to be in Distr.
!   2. Each frequency is separately addressed. The routine is capable of
!   performing an inplace evaluation (if cGridS and cGridD point to the same
!   data array), and it has special memory-efficient implementations for 
!   a source plane and for symmetry in the y-dimension.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cSpaceS   i   type(Space)   Source space - used only for the parameters
!   cGridS    i   type(Grid)    Source grid  - contains the source data array
!   cSpaceD   i   type(Space)   Destination space -used only for the parameters
!   cGridD    io  type(Grid)    Destination grid - upon exiting the subroutine,
!                               contains the field data array 
!
		type(Space), intent(in) :: cSpaceS
		type(Grid), intent(in) :: cGridS
		type(Space), intent(in) :: cSpaceD
		type(Grid), intent(inout) :: cGridD 
	
! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   cSliceS   type(Space)   XYZ-slice of a single frequency of the source, 
!                           including wraparound regions in all dimensions
!   cSliceD   type(Space)   XYZ-slice of a single frequency of the field or the
!                           Green's function, including wraparound regions in 
!                           all dimensions
!   iDiffT        i8b   Difference in time steps between the destination space
!                       and the source space
!   iDimT         i8b   Time dimension of the destination space and the source
!                       space
!   iNumYBlocks   i8b   Number of blocks in the y-dimension: 2 if y-symmetry is
!                       employed, 1 if it is not.
!   iYBlocks      i8b   Loop counter over all blocks in the y-dimension
!   iOmega        i8b   Loop counter over all frequencies on this processor
!   iLx           i8b   Loop counter over all indices in the x-dimension
!   iLy           i8b   Loop counter over all indices in the y-dimension
!   iLz           i8b   Loop counter over all indices in the z-dimension
!   iIndS         i8b   Index into the data array of the source space
!   iIndG         i8b   Index into the data array of the destination space
!   iIndSUpdate   i8b   Increase of iIndS when looping over all z-positions
!                       Normally, this is equal to iD2ZS, but for the primary
!                       source plane and the field plane, iIndSUpdate = 0
!   dOmega         dp   Step size in the angular frequency axis
!   dCorrection    dp   Correction factor for the forward FFT's
!   dKcutoff       dp   Spatial angular cutoff frequency
!   acTemp       char   Temporary char array for output log messages
!
		type(Space) ::		cSliceS, cSliceD
		integer(i8b) ::		iDiffT, iDimT, iNumYBlocks, iYBlocks, &
				 iOmega, iLx, iLy, iLz, iIndS, iIndG, iIndSUpdate
		real(dp) ::			dOmega, dCorrection, dKcutoff
		character(len=1024)::		acTemp;
!		logical(lgt) ::				lStudy
	
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
!   InitConvolutionSpaces
!   InitGrid
!   GridDistr2CreateEmpty
!   Distr2ObtainXYZBlock
!   MapVx
!   MapVy
!   GreenXYZBlockFreqT_lossless
!   GreenXYZBlockFreqT_lossy
!   TransformXY
!   TransformXYZ
!   TransformXYInv
!   TransformXYZInv
!   MapVxInv
!   MapVyInv
!   Distr2PutXYZBlock
!   Distr2ObtainYMirroredXYZBlock
!   Distr2AddToXYZBlock
!   DestructSpace
!   
! *****************************************************************************
!
!   PSEUDOCODE
!
!   Initialize source and destination XYZ slices that include the wraparound
!     regions in all dimensions
!   If symmetry in Y is employed, set the number of blocks in Y to 2, otherwise
!     set to 1
!   For each frequency on this processor do
!       Copy the XYZ data at the current frequency from the source space to 
!         the source XYZ slice
!       For all blocks in Y do
!           Map the source XYZ slice in x and y
!           Fill the destination XYZ slice with the lossless or the lossy 
!             Green's function
!           Map the destination XYZ slice in x and y
!           If the source space contains a primary source or a field plane
!               Transform the source and destination XYZ slices in x and y
!               Apply Fourier correction factor for each Forward transform
!           Else the source space contains a volume contrast source
!               Transform the source and destination XYZ slices in x, y and z
!               Apply Fourier correction factor for each Forward transform
!           Multiply the source and destination XYZ slices elementwise, and
!             put the result into the destination XYZ slice
!           If the source space contains a primary source or a field plane
!               Inverse transform the destination XYZ slice in x and y
!           Else the source space contains a volume contrast source
!               Inverse transform the destination XYZ slice in x, y and z
!           Inverse map the destination XYZ slice in x and y
!           If symmetry in Y is not employed
!               Copy the destination XYZ slice to the current frequency in the
!                 destination space
!           Else
!               If it is the first Y-block
!                   First obtain the XYZ data at the current frequency from the
!                     source space to the source XYZ slice in an order mirrored 
!                     in y and displaced in y to the negative y-axis
!                   Then copy the destination XYZ slice to the current 
!                     frequency in the destination space (order of these two 
!                     operations is important for inplace evaluation)
!               If it is the second Y-block
!                   Add the destination XYZ slice to the current frequency in
!                   the destination space
!   De-initialize source and destination XYZ slices           
!
! =============================================================================

		call PrintToLog("ComputeConvolution", 3)

!!		if(pcSpaceS%iSpaceIdentifier==iSI_CONTRASTSOURCE.and.pcSpaceS%iStartZ==1.and.pcSpaceD%iStartZ>1) then
!		if(pcSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
!			lStudy=.true.
!		else
!			lStudy=.false.
!		end if

		! 3) Initialize convolution slice spaces
		call PrintToLog("Allocate Slice spaces S and D", 3)
		call InitConvolutionSpaces(cSliceS, cSliceD, cSpaceS, cSpaceD)
	    call InitGrid(cSliceS, cGridS%iProcN, cGridS%iProcID, (/ .false., .false., .false., .false./));
	    call InitGrid(cSliceD, cGridD%iProcN, cGridD%iProcID, (/ .false., .false., .false., .false./));
		call GridDistr2CreateEmpty(cSliceS%cGrid);
		call GridDistr2CreateEmpty(cSliceD%cGrid);
!		print *, " Test Function S"
!			print *, size(cGridS%pacD2)
!			print *, " Test Function D"
!			print *, size(cGridD%pacD2)

! call PrintToLog("Variables in cSliceS",-1)
! call SpaceInfo(cSliceS)
! call PrintToLog("Variables in cSliceD",-1)
! call SpaceInfo(cSliceD)

		! 4) Now we perform the convolution on each xyz-block we have stored locally and store the result
		call PrintToLog("Start the inner loop of the convolution", 3)
		! For the T axis we assume that pcSpaceD%iDimT = pcSpaceS%iDimT and pcSpaceD%dTanT=pcSpaceS%dTanT
		! These assumptions are already important for the mapping and transform in T in dist. 1..
		iDiffT = cSpaceD%iStartT - cSpaceS%iStartT
		iDimT = cSpaceS%iDimT
		dKcutoff = two_pi * ( 1.0_dp + 0.5_dp/cSpaceS%iDimX ) * cSpaceS%dFnyq

		!Y symmetry: if true, then we calculate two blocks of a contrast source
		!Not for ISI_PRIMARYSOURCE: not necessary because when evaluating the primary-source-to-field
        ! convolution we are not close to maximum memory anyway
		if ((cSpaceS%iSpaceIdentifier==ISI_CONTRASTSOURCE.or.cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE)&
			.and.(cSpaceS%bYSymm .EQV. .true.)) then
			iNumYBlocks=2
		else
			iNumYBlocks=1
		end if

!if(lStudy) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSpaceS1",&
!					"cSpace", &
!					cSpaceS, &
!					(/ 0_i8b, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ cSpaceS%iDimT+1, cSpaceS%iDimX, cSpaceS%iDimY, cSpaceS%iDimZ /), &
!					.false.);
!end if

		!Main loop in the convolution: for each omega, calculate the convolution in the 
		! spectral domain
		do iOmega = 0, cGridS%iD2LocN-1
			write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iOmega, cGridS.iD2LocN-1
			call PrintToLog(acTemp,3);

			dOmega = two_pi * cGridS%aiD2Loc(1+iOmega) * cSpaceS%dFnyq / cSpaceS%iDimT
			!KH change for periodical T to 
			!KH dOmega = two_pi * cGridS%aiD2Loc(1+iOmega) * 2.0_dp*cSpaceS%dFnyq/cSpaceS%iDimT

			call PrintToLog("Obtain Source slice", 4)
!			print *, " Punto di Interesse Odierno 4 New"
			
			call Distr2ObtainXYZBlock(cGridS, cSliceS%cGrid, iOmega)   !!! L.D. 11-11-2011
!            print *, " Punto di Interesse Odierno 5"
			do iYblocks=1,iNumYBlocks !stays 1 if bYSymm==.false.

				!Obtain and map Source slice
				call PrintToLog("Map Source slice", 4)
				call MapVx(cSliceS)
				call MapVy(cSliceS)

				!Obtain and map Green's function slice
				if(cMediumParams%a_att==0.0_dp) then
					call PrintToLog("Obtain and map lossless Green's function slice", 4)
	 				call GreenXYZBlockFreqT_lossless(cSliceD, cSliceS, dOmega, dKcutoff, iDiffT, iDimT)
				else
					call PrintToLog("Obtain and map lossy Green's function slice", 4)
	 				call GreenXYZBlockFreqT_lossy(cSliceD, cSliceS, dOmega, dKcutoff, iDiffT, iDimT)
				end if

!if(lStudy.and.pcGridD%aiD2Loc(1+iOmega)==36.and.iYblocks==1) then
!
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSliceD0",&
!					"cSlice", &
!					cSliceD, &
!					(/ cSliceD%cGrid%iProcID, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ 1_i8b, cSliceD%iDimX, cSliceD%iDimY, cSliceD%iDimZ /), &
!					.false.);
!call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "cSliceS0",&
!					"cSlice", &
!					cSliceS, &
!					(/ cSliceD%cGrid%iProcID, 0_i8b, 0_i8b, 0_i8b /), &
!					(/ 1_i8b, cSliceS%iDimX, cSliceS%iDimY, cSliceS%iDimZ /), &
!					.false.);
!
!end if

				call MapVx(cSliceD);
				call MapVy(cSliceD);

				if(cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .or. cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
				!iDimZ=1, don't transform in z

					call PrintToLog("Transform Source and Green's function slice in XY", 4)
					call TransformXY(cSliceS%cGrid)
					call TransformXY(cSliceD%cGrid)

					!All factors for the forward transforms of G and S and the inverse transform combined 
					!dCorrection	= cSliceD%dDx**4 * cSliceD%dDt * cSliceD%dDf * cSliceD%dDkx * cSliceD%dDky;
                    ! or equivalently:
					dCorrection	= cSliceD%dDx**2/( real(2*cSpaceS%iDimT * cSliceD%iDimX,dp) * cSliceD%iDimY); 
					!the real() command is to prevent integer overflow error from occurring
					iIndSupdate = 0

				else
				!iSI_CONTRASTSOURCE, iDimZ can be anything...

					call PrintToLog("Transform Source and Green's function slice in XYZ", 4)
!					print *, "Punto di Interesse Odierno 3 Bis New"
					call TransformXYZ(cSliceS%cGrid,.false.)
!					print *, "Punto di Interesse Odierno 4"
					call TransformXYZ(cSliceD%cGrid,.true.)

					!All factors for the forward transforms of G and S and the inverse transform combined 
					!dCorrection	= cSliceD%dDx**6 * cSliceD%dDt * cSliceD%dDf * cSliceD%dDkx * cSliceD%dDky * cSliceD%dDkz;
                    ! or equivalently:
					dCorrection	= cSliceD%dDx**3/( real(2*cSpaceS%iDimT * cSliceD%iDimX,dp) * cSliceD%iDimY * cSliceD%iDimZ); !should give the same result
					!the real command is to prevent integer overflow error from occurring
					iIndSupdate = cSliceS%cgrid%iD2ZS
		
				end if

				! Take the inproduct of cGridS and of cGridG. Put them in distribution 2 of cGridG.
				! since we also have to multiply the result of the Fourier transforms with a constant factor, 
				! we'll do  it all in one go and multiply each value with the correction factor dCorrection. Here 
				call PrintToLog("Compute the inproduct for the two xyz-blocks", 4)
				do iLx = 0, cSliceD%cgrid%iD2XL-1
					do iLy = 0, cSliceD%cgrid%iD2YL-1
						iIndS	=	1 + iLx * cSliceD%cgrid%iD2XS + iLy * cSliceD%cgrid%iD2YS;
						iIndG	=	1 + iLx * cSliceD%cgrid%iD2XS + iLy * cSliceD%cgrid%iD2YS;
						do iLz = 0, cSliceD%cgrid%iD2ZL-1
							cSliceD%cGrid%pacD2(iIndG)	= CONJG(cSliceD%cGrid%pacD2(iIndG)) * &
														  cSliceS%cGrid%pacD2(iIndS) * dCorrection;
							!cSliceD%cGrid%pacD2(iIndG)	= (cSliceD%cGrid%pacD2(iIndG)) * &
							!							  cSliceS%cGrid%pacD2(iIndS) * dCorrection;
							
							iIndS	= iIndS + iIndSupdate
							iIndG	= iIndG + cSliceD%cgrid%iD2ZS;
							!if iSI_PRIMARYSOURCE, iIndsSupdate=0 and we therefore stick to z=0...
						end do
					end do
				end do

				! or simply(for contrastsource only..):
				!cSliceD%cGrid%pacD2	= cSliceD%cGrid%pacD2 * cSliceS%cGrid%pacD2 * dCorrection;

				! 5) Take the inverse transform of cGridP with respect to the XY axes

				if(cSpaceS%iSpaceIdentifier==iSI_PRIMARYSOURCE .or. cSpaceS%iSpaceIdentifier==iSI_FIELDSLICE) then
				!iDimZ=1, don't transform in z

					call PrintToLog("Inverse transform Field slice in XY", 4)
					call TransformXYInv(cSliceD%cGrid)

				else

					call PrintToLog("Inverse transform Field slice in XYZ", 4)
					call TransformXYZInv(cSliceD%cGrid, .false.)

				end if

 				call MapVxInv(cSliceD)
 				call MapVyInv(cSliceD)

				!Storing of the block: depends on whether we have Y symmetry...
				if (cSpaceS%bYSymm .EQV. .false.) then

					! Put the slice into the destination grid
					call Distr2PutXYZBlock(cSliceD%cGrid, cGridD, iOmega)

				else

					if (iYBlocks==1) then

						!Y symmetry: fill cSliceS with a mirrored version of pcSpaceS
						! BEFORE pcSpaceD is written (it could be the same space...)
						call PrintToLog("Obtain Mirrored Source slice", 4)
						call Distr2ObtainYMirroredXYZBlock(cGridS, cSliceS%cGrid,  iOmega)
						! set iStartY parameter of cSliceS to get the negative starting point
						cSliceS%iStartY=-(cSpaceS%iDimY-1)
						! Put the slice into the destination grid
						call Distr2PutXYZBlock(cSliceD%cGrid, cGridD, iOmega)

					else

						!add result to cSpaceD
						call Distr2AddToXYZBlock(cSliceD%cGrid, cGridD, iOmega)
						! reset iStartY parameter of cSliceS to normal
						cSliceS%iStartY=0

					end if

				end if

			end do 

		end do

		! Free the slices and if necessary the temporary and original source space 
		call DestructSpace(cSliceS)
		call DestructSpace(cSliceD)
		call PrintToLog("Slice spaces S and D destructed", 3)

END SUBROUTINE ComputeConvolutionConj

SUBROUTINE GreenXYZBlockFreqT_lossless(cGreen, cSource, dOmega, dKcutoff, iDiffT, iDimT)
	
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
!   The subroutine GreenXYZBlockFreqT_lossless computes the value of the 
!   spatiotemporally filtered and windowed Lossless Green's function (i.e. free
!   space Green's function) on a beam space in x/y/z for a specific value of 
!   temporal angular frequency omega. It takes into account the skewness of the
!   source and field beams in T/X/Y and the translations in T/X/Y/Z from the 
!   field beam with respect to the source beam. The axes contain the values on 
!   positive as well as negative positions, which are ordered in the FFT 
!   fashion - start with zero, positive positions, negative positions.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGreen  io   type(Space)   Space of type CONVOLUTIONSLICE that will contain
!                              the Green's function over the entire space for a
!                              single frequency; its location, size and 
!                              skewness parameters are taken to be those of the
!                              field space
!   cSource  i   type(Space)   Space of type CONVOLUTIONSLICE containing the 
!                              source location, size and skewness parameters.
!   dOmega     i     dp        temporal angular frequency for which the Green's
!                              function is to be determined (relative to freq0)
!   dKcutoff   i     dp        spatial angular cutoff frequency (re k0 = w0/c0)
!   iDiffT     i     i8b       translational difference between the time axis 
!                              of the source space and the field space
!   iDimT      i     i8b       size of the temporal axis
!
	type(Space), target, intent(inout)::	cGreen;
	type(Space), target, intent(in)::	cSource;
	real(dp), intent(in) ::			dOmega, dKcutoff
	integer(i8b), intent(in) :: 	iDiffT, iDimT
	
! *****************************************************************************
!
!   LOCAL PARAMETERS      
!
!   pcGrid   type(Grid)   pointer to the grid in the Green's function space
!   iLX   i8b      array index in x
!   iLY   i8b      array index in y
!   iLZ   i8b      array index in z
!   iIndZ   i8b    beam index in z
!   viX   i8b      Array giving the global index in x as a function of the 
!                  array index in x
!   viY   i8b      Array giving the global index in y as a function of the 
!                  array index in y
!   viZ   i8b      Array giving the global index in z as a function of the 
!                  array index in z
!   i   i8b        temporary counter
!   dStartT   dp   start of the window in the temporal dimension
!   dEndT     dp   end of the window in the temporal dimension
!   dRad      dp   radius, or |x| in the Green's function
!   dShiftFactor   dpc  Factor that accounts for the time difference between
!                       source and destination space
!   error   i8b   Error number
!   acTemp       char   Temporary char array for output log messages
!
!   pcGrid   type(Grid)  pointer to the grid in the Green's function space
!   iLX        i8b       array index in x
!   iLY        i8b       array index in y
!   iLZ        i8b       array index in z
!   iIndZ      i8b       beam index in z
!   viX        i8b       Array giving the global index in x as a function of 
!                        the array index in x
!   viY        i8b       Array giving the global index in y as a function of
!                        the array index in y
!   viZ        i8b       Array giving the global index in z as a function of
!                        the array index in z
!   i          i8b       temporary counter
!   dStartT     dp       start of the window in the temporal dimension
!   dEndT       dp       end of the window in the temporal dimension
!   dRad        dp       radius, or |x| in the Green's function
!   dShiftFactor   dpc   Factor that accounts for the time difference between
!                        source and destination space
!   dSi1        dp       Sine integral temporary variable 1
!   dSi2        dp       Sine integral temporary variable 2
!   dCi1        dp       Cosine integral temporary variable 1
!   dCi2        dp       Cosine integral temporary variable 2
!   error      i8b       Error number
!   acTemp     char      Temporary char array for output log messages
!
	type(Grid), pointer ::			pcGrid;
	integer(i8b), allocatable ::	viX(:), viY(:), viZ(:);
	integer(i8b) :: 		i, iLi, iLX, iLY, iLZ, iIndZ;
	real(dp) ::				dStartT, dEndT, dRad
	complex(dpc) ::			dShiftFactor;
	real(dp) ::				dSi1, dSi2, dCi1, dCi2;
	integer(i8b) :: 		error
	character(1024) ::		acTemp
	
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
!   SwStartAndCount
!   PrintToLog
!   iBeamOffSetT
!   iBeamOffSetX
!   iBeamOffSetY
!   cisi4
!   SWstop
!   
! =============================================================================
    
   	call SwStartAndCount(cswGreen)

	call PrintToLog("GreenXYZBlockFreqT_lossless", 5);

	! For ease of notation we introduce this pointer
	pcGrid=>cGreen%cGrid;
	
	! Initialize the X, Y and Z axes.
	allocate(viX(cGreen%iDimX));
	allocate(viY(cGreen%iDimY));
	allocate(viZ(cGreen%iDimZ));

	viZ		= (/ (i, i=0,cGreen%iBeamIndexStartZ+cGreen%iDimZ-1), (i, i=cGreen%iBeamIndexStartZ,-1) /) &
				+ cGreen%iStartZ-cSource%iStartZ;

	pcGrid%pacD2	= 0;
	do iLz=0,cGreen%iDimZ-1
		!Take care, iLz is the array index, iIndz is the beam index, viZ is the global index...
		iIndZ = mod(iLz - cGreen%iBeamIndexStartZ, cGreen%iDimZ) + cGreen%iBeamIndexStartZ;

		! Calculate the time interval we're interested in
!KH Change here for periodical T
		dEndT		= (iDiffT + iBeamOffSetT(cGreen,iIndZ) + iDimT + 0.5) * cGreen%dDt;
		dStartT		= (iDiffT + iBeamOffSetT(cGreen,iIndZ) - iDimT + 0.5) * cGreen%dDt;

		!KH The iBeamOffset needs to be added here, because the 
		!KH Green's function is defined in the (still unmapped) beam coordinates of the Green space.
		!KH Mapping needs to be done outside the GreenXYZBlockFreqT subroutine
		viX		= (/ (i, i=0,cGreen%iBeamIndexStartX+cGreen%iDimX-1), (i, i=cGreen%iBeamIndexStartX,-1) /) &
						+ iBeamOffSetX(cGreen,iIndZ) + cGreen%iStartX - cSource%iStartX;
		viY		= (/ (i, i=0,cGreen%iBeamIndexStartY+cGreen%iDimY-1), (i, i=cGreen%iBeamIndexStartY,-1) /) &
						+ iBeamOffSetY(cGreen,iIndZ) + cGreen%iStartY - cSource%iStartY;

		! if we include the zero axis in time, we have to include a filter factor
		!KH Change here for periodical T, always use filter term
		if (dStartT <= 0 .and. dEndT > 0) then
			do iLX=0,cGreen%iDimX-1
				do iLY=0,cGreen%iDimY-1
					dRad =  cGreen%dDx * sqrt(real(viX(1+iLX)**2 + viY(1+iLY)**2 + viZ(1+iLZ)**2,dp))
					call cisi4((dKcutoff-dOmega)*dRad, dCi1, dSi1)
					call cisi4((dKcutoff+dOmega)*dRad, dCi2, dSi2)
					if (dRad == 0) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								=  (dKcutoff/(2.0_dp*pi**2)&
									+ dOmega/(4.0_dp*pi**2)*log((dKcutoff-dOmega)/(dKcutoff+dOmega)) &
									- im*dOmega/(4.0_dp*pi))
					else if ((dRad>=dStartT) .and. (dRad<=dEndT)) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
							=(1.0_dp/(dRad*4.0_dp*pi)* &
								( cos(dOmega*dRad)*(dSi1+dSi2)&
								 +sin(dOmega*dRad)*(dCi1-dCi2-im*pi))/pi);
					else
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
							= (1.0_dp/(dRad*4.0_dp*pi)* &
								( cos(dOmega*dRad)*(dSi1+dSi2-pi) &
								 +sin(dOmega*dRad)*(dCi1-dCi2))/pi);
					end if
					
				end do
			end do
		else
		!Otherwise, use the unfiltered chopped Green's function
			do iLX=0,cGreen%iDimX-1
				do iLY=0,cGreen%iDimY-1
					dRad =  cGreen%dDx * sqrt(real(viX(1+iLX)**2 + viY(1+iLY)**2 + viZ(1+iLZ)**2,dp))
					if ((dRad >= dStartT).and.(dRad <= dEndT)) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 	(1.0_dp/(dRad*4.0_dp*pi)*(cos(dOmega * dRad)-im*sin(dOmega * dRad)));
					else 
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 0.0_dp;
					end if
				end do
			end do	
		end if
	
	end do

	!account for overall t shift factor
	dShiftFactor=exp(im*dOmega*iDiffT*cGreen%dDt)
	pcGrid.pacD2 = pcGrid.pacD2 * dShiftFactor

	deallocate(viX);
	deallocate(viY);
	deallocate(viZ);

	call SWStop(cswGreen)

END SUBROUTINE GreenXYZBlockFreqT_lossless

SUBROUTINE GreenXYZBlockFreqT_lossy(cGreen, cSource, dOmega, dKcutoff, iDiffT, iDimT)
	
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
!   The subroutine GreenXYZBlockFreqT_lossy computes the value of the 
!   spatiotemporally filtered and windowed Lossy Green's function (attenuation
!   of the frequency power law type a*f^b, powers 1<b<=2, with the accompanying
!   dispersion) on a beam space in x/y/z for a specific value of temporal  
!   angular frequency omega. It takes into account the skewness of the source 
!   and field beams in T/X/Y and the translations in T/X/Y/Z from the field 
!   beam with respect to the source beam. The axes contain the values on 
!   positive as well as negative positions, which are ordered in the FFT 
!   fashion - start with zero, positive positions, negative positions.
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   cGreen  io   type(Space)   Space of type CONVOLUTIONSLICE that will contain
!                              the Green's function over the entire space for a
!                              single frequency; its location, size and 
!                              skewness parameters are taken to be those of the
!                              field space
!   cSource  i   type(Space)   Space of type CONVOLUTIONSLICE containing the 
!                              source location, size and skewness parameters.
!   dOmega     i     dp        temporal angular frequency for which the Green's
!                              function is to be determined (relative to freq0)
!   dKcutoff   i     dp        spatial angular cutoff frequency (re k0 = w0/c0)
!   iDiffT     i     i8b       translational difference between the time axis 
!                              of the source space and the field space
!   iDimT      i     i8b       size of the temporal axis
!
	type(Space), target, intent(inout)::	cGreen;
	type(Space), target, intent(in)::	cSource;
	real(dp), intent(in) ::			dOmega, dKcutoff
	integer(i8b), intent(in) :: 	iDiffT, iDimT
	
! *****************************************************************************
!
!   LOCAL PARAMETERS      
!
!   pcGrid   type(Grid)  pointer to the grid in the Green's function space
!   iLX        i8b       array index in x
!   iLY        i8b       array index in y
!   iLZ        i8b       array index in z
!   iIndZ      i8b       beam index in z
!   viX        i8b       Array giving the global index in x as a function of 
!                        the array index in x
!   viY        i8b       Array giving the global index in y as a function of
!                        the array index in y
!   viZ        i8b       Array giving the global index in z as a function of
!                        the array index in z
!   i          i8b       temporary counter
!   dStartT     dp       start of the window in the temporal dimension
!   dEndT       dp       end of the window in the temporal dimension
!   dRad        dp       radius, or |x| in the Green's function
!   dShiftFactor   dpc   Factor that accounts for the time difference between
!                        source and destination space
!   dKangular  dpc       complex wavenumber including losses
!   E1mm       dpc       Primed exp integral temorary variable 1
!   E1mp       dpc       Primed exp integral temorary variable 2
!   E1pm       dpc       Primed exp integral temorary variable 3
!   E1pp       dpc       Primed exp integral temorary variable 4
!   error      i8b       Error number
!   acTemp     char      Temporary char array for output log messages
!
	type(Grid), pointer ::			pcGrid;
	integer(i8b), allocatable ::	viX(:), viY(:), viZ(:);
	integer(i8b) :: 		i, iLi, iLX, iLY, iLZ, iIndZ;
	real(dp) ::				dStartT, dEndT, dRad
	complex(dpc) ::			dShiftFactor;
	complex(dpc) :: dKangular,E1mm,E1mp,E1pm,E1pp
	integer(i8b) :: 		error
	character(1024) ::		acTemp

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
!   SWstartandcount
!   PrintToLog
!   iBeamOffSetT
!   iBeamOffSetX
!   iBeamOffSetY
!   Expint_prime
!   SWstop
!   
! =============================================================================

	call SwStartAndCount(cswGreen)

	call PrintToLog("GreenXYZBlockFreqT_lossy", 5);

	! For ease of notation we introduce this pointer
	pcGrid=>cGreen%cGrid;
	
	! Initialize the X, Y and Z axes. The order we use here is described in the description of this function
	allocate(viX(cGreen%iDimX));
	allocate(viY(cGreen%iDimY));
	allocate(viZ(cGreen%iDimZ));

	!dKangular is the complex wavenumber (normalized wrt c0/f0) that includes losses
	dKangular = cMediumParams%c0/cModelParams%freq0 * (dOmega*cModelParams%freq0/cMediumParams%c0	&
		+ cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*dOmega*cModelParams%freq0		&
		* (abs(dOmega*cModelParams%freq0)**(cMediumParams%b_att-1) -				&
		(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1))						&
		- im*cMediumParams%alpha0_att*abs(dOmega*cModelParams%freq0)**cMediumParams%b_att)

	viZ		= (/ (i, i=0,cGreen%iBeamIndexStartZ+cGreen%iDimZ-1), (i, i=cGreen%iBeamIndexStartZ,-1) /) &
				+ cGreen%iStartZ-cSource%iStartZ;

	pcGrid%pacD2	= 0;
	do iLz=0,cGreen%iDimZ-1
		!Take care, iLz is the array index, iIndz is the beam index, viZ is the global index...
		iIndZ = mod(iLz - cGreen%iBeamIndexStartZ, cGreen%iDimZ) + cGreen%iBeamIndexStartZ;

		! Calculate the time interval we're interested in
!KH Change here for periodical T
		dEndT		= (iDiffT + iBeamOffSetT(cGreen,iIndZ) + iDimT + 0.5) * cGreen%dDt;
		dStartT		= (iDiffT + iBeamOffSetT(cGreen,iIndZ) - iDimT + 0.5) * cGreen%dDt;

		!KH The iBeamOffset needs to be added here, because the 
		!KH Green's function is defined in the (still unmapped) beam coordinates of the Green space.
		!KH Mapping needs to be done outside the GreenXYZBlockFreqT subroutine
		viX		= (/ (i, i=0,cGreen%iBeamIndexStartX+cGreen%iDimX-1), (i, i=cGreen%iBeamIndexStartX,-1) /) &
						+ iBeamOffSetX(cGreen,iIndZ) + cGreen%iStartX - cSource%iStartX;
		viY		= (/ (i, i=0,cGreen%iBeamIndexStartY+cGreen%iDimY-1), (i, i=cGreen%iBeamIndexStartY,-1) /) &
						+ iBeamOffSetY(cGreen,iIndZ) + cGreen%iStartY - cSource%iStartY;

		! if we include the zero axis in time, we have to include a filter factor
!KH Change here for periodical T, then always use filter term
		if (dStartT <= 0 .and. dEndT > 0) then
			do iLX=0,cGreen%iDimX-1
				do iLY=0,cGreen%iDimY-1
					dRad =  cGreen%dDx * sqrt(real(viX(1+iLX)**2 + viY(1+iLY)**2 + viZ(1+iLZ)**2,dp))
					call Expint_prime(-im*(dKcutoff-dKangular)*dRad,E1mm,error)
					call Expint_prime( im*(dKcutoff-dKangular)*dRad,E1pm,error)
					call Expint_prime(-im*(dKcutoff+dKangular)*dRad,E1mp,error)
					call Expint_prime( im*(dKcutoff+dKangular)*dRad,E1pp,error)
					if (dRad == 0) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								=  (dKcutoff/(2.0_dp*pi**2)&
									+ dKangular/(4.0_dp*pi**2)*log((dKcutoff-dKangular)/(dKcutoff+dKangular)) &
									- im*dKangular/(4.0_dp*pi))
					else if (dRad<=dEndT) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 1.0_dp/(4.0_dp*pi*dRad)*(exp(-im*dKangular*dRad)	&
										+(-exp( im*dKcutoff*dRad)*(E1mm+E1mp)&
										  +exp(-im*dKcutoff*dRad)*(E1Pm+E1pp))/(im*two_pi))	
					else
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 1.0_dp/(8.0_dp*pi**2*im*dRad)*(				&
										-exp( im*dKcutoff*dRad)*(E1mm+E1mp)		&
										+exp(-im*dKcutoff*dRad)*(E1Pm+E1pp))	
					end if
				end do
			end do
		else
		!Otherwise, use the unfiltered chopped Green's function
			do iLX=0,cGreen%iDimX-1
				do iLY=0,cGreen%iDimY-1
					dRad =  cGreen%dDx * sqrt(real(viX(1+iLX)**2 + viY(1+iLY)**2 + viZ(1+iLZ)**2,dp))
					if ((dRad >= dStartT).and.(dRad <= dEndT)) then
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 1.0_dp/(dRad*4.0_dp*pi)*(exp(-im * dKangular * dRad));
					else 
						pcGrid.pacD2(1+iLX*pcGrid%iD2XS+iLY*pcGrid%iD2YS+iLZ*pcGrid%iD2ZS)  &
								= 0.0_dp;
					end if
				end do
			end do	
		end if
	end do
	
	!account for overall t shift factor
	dShiftFactor=exp(im*dOmega*iDiffT*cGreen%dDt)
	pcGrid.pacD2 = pcGrid.pacD2 * dShiftFactor
	
	deallocate(viX);
	deallocate(viY);
	deallocate(viZ);

	call SWStop(cswGreen)
	
END SUBROUTINE GreenXYZBlockFreqT_lossy

END MODULE ParnacConvolution

