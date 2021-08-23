MODULE ParnacBubbleContrast

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
  USE ParnacDataRedist, ONLY : &
       ReorderDistr0ToDistr1,ReorderDistr1ToDistr0, &
	   ReorderDistr1ToDistr2,ReorderDistr2ToDistr1
  USE ParnacDataSpaceInit, ONLY : &
       iBeamOffsetX, iBeamOffsetY, InitSaveSlices
  USE ParnacTransformFilter, ONLY : &
       TransformT, TransformTInv, &
       TransformT_sml, TransformTInv_sml, &
       dTaperingWindow, INTERP1DFREQ, GREENS1DFREQ
  USE ParnacParamInit, ONLY :  BubbleInit
  USE ParnacOutput, ONLY : ExportSlice
  USE SpecialFun, ONLY : &
       dSinc,INTERP1D, INTERP3D, LINSPACE, FIND_CLOSEST_DIVISOR
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

    USE ParnacParamInit, ONLY : &
         BubbleInit
    USE RP
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

    type(Space), target, intent(inout)  ::	    cSpace    

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

    integer 	         		::  	iErr, iStat(MPI_STATUS_SIZE), iREQUEST
    integer(i8b)         		:: 		iBeamIndex(3),iBeamIndexN(3), XYZindex , iTargetIndex(3)
    integer(i8b)				::		ones(8,3), XYZ_i1(8), p(10), XYZ_i2(size(p,1)**3), NeighbGridPoints(size(p,1)**3,3), iBeamIndexNGP(size(p,1)**3,3)
    real(dp)             		:: 		dGlobXYZ(3), GlobalPos(8,3) , xyz_factor(3), diff(8,3), Domain_Range(2,3), BubbleEnhancedP(1E3)
    real(dp)             		:: 		dRadius(1), dLambdaMM, dK_c
    integer(i8b)         		:: 		iDimX, iDimY, iDimZ, iDimT
    integer(i8b)         		:: 		iBubble,BubbleProc , iIndex, i, j, l, i_neighbgp
    integer						::		M, N, K, ibDimT
    integer						::		i_loc, log_isnan
    integer						::		Bubble_disps(cSpace%cGrid%iProcN),Bubble_locdisps(cSpace%cGrid%iProcN),Bubble_condisps(cSpace%cGrid%iProcN), Bubble_count(cSpace%cGrid%iProcN)

    character(len = 1024)		::  	actemp, sBubbleDir 

    integer(i8b),allocatable	::  	arBufferXYZ(:,:), BubbleID(:), tempbuffer(:) , NeighbGPIndex(:)
    real(dp),allocatable		::  	dBubbleContrast(:), dBubbleLocationN(:,:) 
    real(dp),allocatable		::		temploc(:), tempcon(:), temp_globloc(:), LocPres(:,:)

    ! Parameters for Rayleigh - Plesset Solver
    integer(i8b),parameter		::  	n_pad = 8, OverSampling_Grid = 1

    real(dp)		 			::  	RealPressure(cSpace%iDimT)				, RealTime(cSpace%iDimT)			,  R(cSpace%iDimT,3) 			
    real(dp)             		::  	RealPressurePad(cSpace%iDimT *n_pad)  	, RealTimePad(cSpace%iDimT *n_pad)	,  R_Pad(cSpace%iDimT *n_pad,3)
    real(dp)             		::  	RealPressure_i1(8,cSpace%iDimT)			, RealPressure_Own(8,cSpace%iDimT)	,  Interp_PressureBL(cSpace%iDimT)
    real(dp)             		::  	V_dd_norm(OverSampling_Grid*cSpace%iDimT)	,  V_dd_Pad(cSpace%iDimT *n_pad), RealTimePadOut(cSpace%iDimT *n_pad)
    LOGICAL(lgt)            	::  	exists, exists_loc, result ,last_iter

    !Filtering and windowing parameters 
    integer(i8b)				::		iDimW, ibDimW, iStart
    real(dp)					::		dTaperSupportWindow(cSpace%iDimT), dTaperSupportWindowN(cSpace%iDimT*n_pad), dTaperMaxFreqWindow(cSpace%iDimT), dLeftBand, dRightBand
    real(dp) 					::		dDOmega, dFFTFactor, dMultFactor(cSpace%iDimT/2+1)

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
    !   CreateBubbleCluster
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
    !10			  Save only the bubbles that have neighbouring points outside of the processor
    !11	  For each bubble,
    !11   		  Find in which of the processors' grid is the bubble located and continue with this one,
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
    !22	  And correct the value fo the scattered pressure field.
    !
    ! =============================================================================
    call BubbleInit()

    sBubbleDir = 'Bubbles/'
    INQUIRE(FILE=trim(trim(sOutputDir)//trim(sBubbleDir)),EXIST =exists_loc)
    if (.NOT. exists_loc) result = MAKEDIRQQ(trim(trim(sOutputDir)//trim(sBubbleDir)))

    write(acTemp,"('BubbleContrastOperator')");    call PrintToLog(acTemp,2)
    !**************************************** INITIALIZE VALUES **********************************************
    iDimT   = cSpace%iDimT   !t dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
	
    iDimT   = OverSampling_Grid*iDimT
    iDimW   = iDimT/2 + 1
    dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor λ in [mm]
    dK_c = pi / cSpace%dDx

    ibDimT   = cSpace%iDimT ! This should not be changed because ibDimT is Integer

    !Allocate for memory management
    ALLOCATE(arBufferXYZ(cSpace.cGrid.iD0LocN,3))
    ALLOCATE(temploc(BubbleParams%BubbleN*3), BubbleID(BubbleParams%BubbleN),temp_globloc(3+cSpace%iDimT))
    ALLOCATE(dBubbleLocationN(BubbleParams%BubbleN,3) )

    iMemAllocated=iMemAllocated + (PRODUCT(SHAPE(temploc)) + PRODUCT(SHAPE(temp_globloc))  + PRODUCT(SHAPE(dBubbleLocationN)) )*dpS +  ( PRODUCT(SHAPE(arBufferXYZ)) + PRODUCT(SHAPE(BubbleID)) ) * i4bS

    arBufferXYZ = cSpace%cGrid%aiD0Loc ! src: ParnacDataDef,  [X*Y*Z,3] Array with all grid points coordinates

    !Initialize for this variables , close to 0, double precission
    dBubbleLocationN = 1.0D-30
    dBubbleContrast	 = 1.0D-30
	BubbleID 		 = 1.0D-30
    temploc 		 = 1.0D-30       
	temp_globloc 	 = 1.0D-30
    i_loc			 = 0
    !**************************************** CREATE BUBBLE CLUSTER / LOAD BUBBLE PARAMETERS **********************************************

    !!This is only for iBeam =1!! I should change this
    INQUIRE(FILE=trim(trim(sOutputDir)//trim(sBubbleDir)//'Bubble_Location'//int2str(0)),EXIST =exists_loc)
    ! If the file for the first iteration exists , then just load the file from the location
    ! See later for more information about the process of loading
    if (.NOT. exists_loc) then
       write(acTemp,"('Randomly Position the microbubbles in the domain') ");    call PrintToLog(acTemp,2)

       ! Generate the cluster and save result on slices defined in input file
       call CreateBubbleCluster_Random(cSpace,dBubbleLocationN, Domain_Range)
       call SaveExtraBubbleSlices(cSpace, dBubbleLocationN, Domain_Range)
	   
	   ! call RANDOM_SEED() ! Fully random , even at each run of INCS
	   ! call RANDOM_NUMBER(BubbleEnhancedP)
	   ! BubbleEnhancedP = INT(CEILING(BubbleEnhancedP*(size(BubbleEnhancedP,1)-1)+1))
       ! call MPI_BCAST(BubbleEnhancedP,size(BubbleEnhancedP,1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iErr)
       write(acTemp,"('No. of Bubbles consisting the cluster: ', I<INT(log10(real(BubbleParams%BubbleN,dp)))+1>, ' in a volume of ', F5.3, ' [mL].')") BubbleParams%BubbleN, PRODUCT(maxval(dBubbleLocationN,1)-minval(dBubbleLocationN,1))*dLambdaMM**3*1.0D-3;    call PrintToLog(acTemp,2);
    else
       ! Load each Bubble's Position 
       open (11, file=trim(trim(sOutputDir)//trim(sBubbleDir)//trim('Bubble_Location')//int2str(0)), status='OLD')
       do iBubble= 1 ,BubbleParams%BubbleN
          read(11, *) dBubbleLocationN(iBubble,1),dBubbleLocationN(iBubble,2),dBubbleLocationN(iBubble,3)
       end do
       close(11)

       ! Load i_loc which is the number of the scatterers in each cpu domain
       open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_iloc')//int2str(cSpace%cGrid%iProcID)), status='OLD')
       read (8,*) i_loc
       close(8)
       ALLOCATE(tempcon(i_loc*iDimT))

       ! Load the radius of the scatterers
       if (BubbleParams%Distribution == 'polydisperse') then
          open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_Radius')//int2str(0)), status='OLD')
          read (8,*) BubbleParams%R0
          close(8)
       endif
    endif

    !********************************************* START OF TAPPERING *********************************************************
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !    use tapering of the field at start and end to prevent wraparound leakage
       !    build tapering window, the first two periods and the last two periods are tapered
       !    in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       cSpace%cGrid%parD0= cSpace%cGrid%parD0 * (/SPREAD(dTaperSupportWindow,2,cSpace%cGrid%iD0LocN)/)
    end if
    !************************************* START OF 3D INTERPOLATION FOR BUBBLES ***************************************************

    ! Create the global positions for every point in order to find the proper pressure from the interpolation in respect to Bubble Location
    ! Then use message passing interface in order for the processors to communicate with each other the information for the scatterers
    ! that are positioned in corner cases (close to each cpu domain boundaries)

    ! This is done to reduce the memory usage
    call ProcBubbleCount(cSpace,dBubbleLocationN,tempcon,i_loc)
    if (.NOT. ALLOCATED(tempcon)) ALLOCATE(tempcon(i_loc*cSpace%iDimT))
    ALLOCATE(LocPres(8,(3+iDimT)*i_loc))

    ALLOCATE(NeighbGPIndex(size(p,1)**3*BubbleParams%BubbleN))
    iMemAllocated=iMemAllocated + (PRODUCT(SHAPE(LocPres)) + PRODUCT(SHAPE(tempcon)) ) * dpS + (PRODUCT(SHAPE(NeighbGPIndex)) ) * i4bS
    NeighbGPIndex = 0
    LocPres=1D-30
    i_loc=0 
	i_neighbgp =0

    ! This is to generate all the 8 neighbouring points for the trilinear interpolation
    ones = 1 ; ones(:,1) = (/-1,1,-1,1,-1,1,-1,1/); ones(:,2) = (/-1,-1,-1,-1,1,1,1,1/); ones(:,3) = (/-1,-1,1,1,-1,-1,1,1/)
    ! This is to produce results on sufficient gridpoints to have accurate results. By increasing the size of p (even numbers) , the accuracy increases
    p = (/ (i, i=-size(p,1)/2,size(p,1)/2-1) /) ; 
    NeighbGridPoints= RESHAPE( (/ ( ( ( p(i),p(j),p(l) , i=1,size(p,1)) , j =1,size(p,1)), l =1,size(p,1)) /) ,(/ size(p,1)**3,3 /),order=(/2,1/) )

    write(acTemp,"('Calculating Trilinear Interpolation') ");    call PrintToLog(acTemp,2)
    do iBubble =  1,BubbleParams%BubbleN
       ! Determine Beam Index of the point source
       ! In ParnacDataDef it is specified that cSpace%iStartT/X/Y/Z should be included 
       ! The difference here is that the dBubbleLocationN parameter is in the normalized space
       ! So if we change this in real space then 
       ! dBubbleLocation = (dBubbleLocationN + (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaNN
       ! So there is no difference if the starting position is not added in all the equations
       iBeamIndex(3) = nint(dBubbleLocationN(iBubble,3)/cSpace%dDx) 
       iBeamIndex(1) = nint(dBubbleLocationN(iBubble,1)/cSpace%dDx) - iBeamOffsetX(cSpace,iBeamIndex(3))
       iBeamIndex(2) = nint(dBubbleLocationN(iBubble,2)/cSpace%dDx) - iBeamOffsetY(cSpace,iBeamIndex(3))

       !Find the Index in the arBufferXYZ of the BubbleLocation
       XYZIndex = (iBeamIndex(3)*iDimY+iBeamIndex(2))*iDimX+iBeamIndex(1)-cSpace%cGrid%iProcID * cSpace%cGrid%iD0LocN
       if (XYZIndex>=0 .AND. XYZIndex<cSpace%cGrid%iD0LocN) i_loc = i_loc+1  ! Update i_loc if Bubble inside processor

       ! Determine current location in global grid
       ! Compare the global coordinates with the scatterer location
       dGlobXYZ(1) = (iBeamIndex(1) + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
       dGlobXYZ(2) = (iBeamIndex(2) + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
       dGlobXYZ(3) = (iBeamIndex(3) ) * cSpace.dDx 

       dRadius(1) = sqrt(sum((dGlobXYZ-dBubbleLocationN(iBubble,:))**2)) 

       ! If the scatterer is not in the grid point (dRadius(1)>1e-10),
       ! then calculate pressure based on the 8 neighbouring points, trilinear interpolation
       if (dRadius(1)>1.0D-10) then		  
          !each if the point is found in a processor different that XYZIndex location
          ! Xyz_factor is to find if the bubble real location is below or above the location of the grid point closer to the bubble
          xyz_factor = floor(dBubbleLocationN(iBubble,:) - dGlobXYZ(:))+1.0D0/2.0D0	   
		  
          do i =1,8
             iBeamIndexN = iBeamIndex + int(xyz_factor+ones(i,:)*1.0D0/2.0D0)
             XYZ_i1(i) = iBeamIndexN(3)*iDimY*iDimX+iBeamIndexN(2)*iDimX+iBeamIndexN(1) - cSpace%cGrid%iProcID * cSpace%cGrid%iD0LocN

             if (XYZ_i1(i)>=0 .AND. XYZ_i1(i)<cSpace%cGrid%iD0LocN) then

                RealPressure_i1(i,:)= cSpace%cGrid%parD0(1 +XYZ_i1(i)*iDimT: (XYZ_i1(i)+1)*iDimT)
                GlobalPos(i,1) = (arBufferXYZ(XYZ_i1(i)+1,1) + iBeamOffsetX(cSpace,arBufferXYZ(XYZ_i1(i)+1,3))) * cSpace.dDx
                GlobalPos(i,2) = (arBufferXYZ(XYZ_i1(i)+1,2) + iBeamOffsetY(cSpace,arBufferXYZ(XYZ_i1(i)+1,3))) * cSpace.dDx
                GlobalPos(i,3) = (arBufferXYZ(XYZ_i1(i)+1,3) ) * cSpace.dDx 
                temp_globloc(1:3) = GlobalPos(i,:)
                temp_globloc(4:iDimT+3) = RealPressure_i1(i,:)	

                ! Bubbles' neighbours that are positioned in a different processor, are transferred via MPI_Send to the correct cpu
                if (XYZIndex<0 .OR. XYZIndex>=cSpace%cGrid%iD0LocN) then
                   BubbleProc = (iBeamIndex(1)+iBeamIndex(2)*iDimX+iBeamIndex(3)*iDimX*iDimY)/cSpace%cGrid%iD0LocN					
                   call MPI_Send(temp_globloc, iDimT+3, MPI_DOUBLE_PRECISION, BubbleProc, (iBubble-1)*8+i,MPI_COMM_WORLD,iErr)	
                else  
                   LocPres(i,1+(3+iDimT)*(i_loc-1):(3+iDimT)*i_loc) = temp_globloc
                endif
             endif

          enddo          
		  iBeamIndexNGP(:,1) = iBeamIndex(1) + NeighbGridPoints(:,1)+INT(xyz_factor(1)+1.0D0/2.0D0)
          iBeamIndexNGP(:,2) = iBeamIndex(2) + NeighbGridPoints(:,2)+INT(xyz_factor(2)+1.0D0/2.0D0)
          iBeamIndexNGP(:,3) = iBeamIndex(3) + NeighbGridPoints(:,3)+INT(xyz_factor(3)+1.0D0/2.0D0)
          XYZ_i2 = iBeamIndexNGP(:,3)*iDimY*iDimX+iBeamIndexNGP(:,2)*iDimX+iBeamIndexNGP(:,1)
          NeighbGPIndex(i_neighbgp+1:i_neighbgp+size(p,1)**3) = XYZ_i2
          i_neighbgp = i_neighbgp + size(p,1)**3

       else if(.NOT. ANY(NeighbGPIndex(1:i_neighbgp) == XYZIndex) ) then
          i_neighbgp = i_neighbgp+1
          NeighbGPIndex(i_neighbgp) = XYZIndex! Save neighbourloc points in order to calculate pressure only there
       endif
    enddo
	
    DEALLOCATE(arBufferXYZ)

    call I_unista(NeighbGPIndex(1:i_neighbgp), i_neighbgp)
    iMemAllocated=iMemAllocated - ( PRODUCT(SHAPE(NeighbGPIndex)) + PRODUCT(SHAPE(arBufferXYZ)) ) * i4bs
    ! NeighGPIndex has more elements than needed. Here we transfer only the one needed to save memory
    ALLOCATE(tempbuffer(i_neighbgp)); tempbuffer = NeighbGPIndex(1:i_neighbgp) ; DEALLOCATE(NeighbGPIndex) ; 
    ALLOCATE(NeighbGPIndex(i_neighbgp)); NeighbGPIndex = tempbuffer; DEALLOCATE(tempbuffer);
    iMemAllocated=iMemAllocated + PRODUCT(SHAPE(NeighbGPIndex)) * i4bs

    call MPI_BARRIER(MPI_COMM_WORLD, iErr);
    !************************************* 3D INTEPR , RP SOLVER, POINT CONTRAST **************************************************
    call PrintToLog("Calculating Bubble Contrast",2)

    RealTime = [(i , i=cSpace%iStartT,iDimT-1+cSpace%iStartT)]*cSpace%dDt/cModelParams%freq0     ! Array with the real time unnormalized 
    RealTimePad = 1D-30
    call linspace(RealTimePad, RealTime(1), RealTime(size(RealTime,1)), iDimT*n_pad)
	dTaperSupportWindowN = dTaperingWindow(iDimT*n_pad,cSpace%dDt/n_pad,dLeftBand,dRightBand)
	
    ALLOCATE(BubbleParams%P_driv(iDimT*n_pad),BubbleParams%T_driv(iDimT*n_pad))

    i_loc = 1
    do iBubble =  1,BubbleParams%BubbleN
       ! Determine Beam Index of the point source
       iBeamIndex(3) = nint(dBubbleLocationN(iBubble,3)/cSpace%dDx)
       iBeamIndex(1) = nint(dBubbleLocationN(iBubble,1)/cSpace%dDx) - iBeamOffsetX(cSpace,iBeamIndex(3))
       iBeamIndex(2) = nint(dBubbleLocationN(iBubble,2)/cSpace%dDx) - iBeamOffsetY(cSpace,iBeamIndex(3))
		
       ! Find the Index in the arBufferXYZ of the BubbleLocation
       XYZIndex = iBeamIndex(3)*iDimY*iDimX+iBeamIndex(2)*iDimX+iBeamIndex(1)-cSpace%cGrid%iProcID * cSpace%cGrid%iD0LocN
       ! Check if the cpu domain contains this scatterer
       if (XYZIndex>=0 .AND. XYZIndex<cSpace%cGrid%iD0LocN) then
          ! Find the scatterer position and do inside if all the necessary actions because of spatial decomposition
          ! Place the Bubble in respect to i_loc for later use for storing its values to a file
          RealPressure = cSpace%cGrid%parD0(1+XYZIndex*iDimT: (XYZIndex+1)*iDimT);   ! Array buffer with length size of iDimT for every XYZ position

          ! Determine current location in global grid
          dGlobXYZ(1) = (iBeamIndex(1)  + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
          dGlobXYZ(2) = (iBeamIndex(2)  + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
          dGlobXYZ(3) = (iBeamIndex(3)  ) * cSpace.dDx
          ! Compare the global coordinates with the scatterer location
          dRadius(1) = sqrt(sum((dGlobXYZ-dBubbleLocationN(iBubble,:))**2))
          if (dRadius(1)>1.0D-10) then

             do i = 1,8
                ! Save to temp_globloc which is a temporary buffer the position and pressure of the bubbles that are not positioned close to
                ! each processor's grid.
                RealPressure_i1(i,:) = 1D-30;  GlobalPos(i,:) = 1D-30; temp_globloc = 1D-30;
                temp_globloc = LocPres(i,1+(3+iDimT)*(i_loc-1):(3+iDimT)*i_loc)

                !If some bubble is positioned close to the boundaries, receive the sent message
                if (maxval(LocPres(i,4+(3+iDimT)*(i_loc-1):(3+iDimT)*i_loc)) <1D-20 .OR. &
					maxval(abs(dBubbleLocationN(iBubble,:) - LocPres(i,1+(3+iDimT)*(i_loc-1):3+(3+iDimT)*(i_loc-1)))) > cSpace%dDx) then
                   xyz_factor = floor(dBubbleLocationN(iBubble,:) - dGlobXYZ(:))+1.0D0/2.0D0 ;
                   iBeamIndexN = iBeamIndex + int(xyz_factor+ones(i,:)*1.0D0/2.0D0)
                   BubbleProc = (iBeamIndexN(3)*iDimY*iDimX+iBeamIndexN(2)*iDimX+iBeamIndexN(1))/cSpace%cGrid%iD0LocN
                   call MPI_Recv(temp_globloc, iDimT+3, MPI_DOUBLE_PRECISION, BubbleProc, (iBubble-1)*8+i,MPI_COMM_WORLD,iStat,iErr) ! IMprove with MPI_Irecv	
                endif

                ! Assign the values of the temp_globloc to the parameters that will be used for the interpolation
                GlobalPos(i,:) = temp_globloc(1:3)
                RealPressure_i1(i,:) =temp_globloc(4:iDimT+3)

                ! if Location is transferred in a correct way. If the difference with the bubble location is bigger than a space step,
                ! then it means that the transfer did not succeed.
                log_isnan = sum(isnan(GlobalPos(i,:))*-1) + sum( isnan(RealPressure_i1(i,:))*-1)
                if ( maxval(abs(dBubbleLocationN(iBubble,:) - GlobalPos(i,:))) > cSpace%dDx .OR. maxval(RealPressure_i1(i,:))<1D-20 .OR. log_isnan>0) then
                   write(acTemp,"('Error in Trilinear Interpolation in Location , Bubble No. ', I7, ' Neighbouring ',I5, ' .')") iBubble, i
                   call PrintToLog(acTemp,-1) 
                   write(*,*) "Exists : " , exists , XYZIndex ,cSpace%cGrid%iProcID
                   stop
                endif

                RealPressure_Own(i,:) = 1.0D-30;
                diff(i,:) = dBubbleLocationN(iBubble,:) - GlobalPos(i,:)
                if (exists_loc) RealPressure_Own(i,:) = tempcon(1+(i_loc-1)*iDimT:i_loc*iDimT)*(cModelParams%freq0/cMediumParams%c0) *(dK_c/pi)**3

                iTargetIndex(3) = nint(GlobalPos(i,3)/cSpace%dDx) 
                iTargetIndex(1) = nint(GlobalPos(i,1)/cSpace%dDx) - iBeamOffsetX(cSpace,iBeamIndex(3))
                iTargetIndex(2) = nint(GlobalPos(i,2)/cSpace%dDx) - iBeamOffsetY(cSpace,iBeamIndex(3))

                !Needs some more improvement, maybe not fully remove the exact "own" scattering ~ 5-10% error because more points should be included
                call GREENS1DFREQ( RealPressure_Own(i,:) , RealPressure_Own(i,:), cSpace ,iTargetIndex,diff(i,:))
                RealPressure_i1(i,:) = RealPressure_i1(i,:) - RealPressure_Own(i,:)
             enddo
             !Do the 3D interpolation in order to acquire the pressure , RealPressureFinal    
             call INTERP3D(GlobalPos,RealPressure_i1,dBubbleLocationN(iBubble,:),Interp_PressureBL)
             RealPressure = 1D-30
             RealPressure  = Interp_PressureBL
          else
             RealPressure_Own(1,:) = 0.0D0;
             diff(1,:) = dBubbleLocationN(iBubble,:) - dGlobXYZ
             if (exists_loc) RealPressure_Own(1,:) = tempcon(1+(i_loc-1)*iDimT:i_loc*iDimT)*(cModelParams%freq0/cMediumParams%c0) *(dK_c/pi)**3
             ! Needs some more improvement, maybe not fully remove the exact "own" scattering ~ 5-10% error
             call GREENS1DFREQ( RealPressure_Own(1,:) , RealPressure_Own(1,:), cSpace ,iBeamIndex,diff(1,:))
             RealPressure = RealPressure - RealPressure_Own(1,:)
          endif

          !0 for every pressure point for zero-padding later and for initial condition. This should be done for every bubble
          RealPressurePad = 0.0D0 ;  R_Pad = 0.0D0
		  
          ! You can use INTERP1DFREQ (frequency interpolation) instead, but this is more accurate and efficient
          call INTERP1D( RealTime, RealPressure, RealTimePad, RealPressurePad)  
          ! call INTERP1DFREQ(RealPressure,RealPressurePad, cSpace,0)
		  
		  RealPressurePad = RealPressurePad * dTaperSupportWindowN
          ! Open File
		  ! open (7, file=trim(trim(sOutputDir)//trim(sBubbleDir)//trim('output_file_pressure')//int2str(cSpace%cGrid%iProcID)//'.txt'),status = 'UNKNOWN', action='write',position='append')				
          ! write(7,*) iBubble
          ! write(7,*) RealTime, RealPressure,RealTimePad,RealPressurePad
          ! close(7) 
			
          ! Marmottant solver ( ODE Solver) , Result is R_Pad with R(:,1) is the radius , R(:,2) is the velocity (R_dot) 
          ! and R(:,3) is the acceleration (R_dotdot)
          call RP_SOLVER(RealPressurePad, RealTimePad,RealTime(1), RealTime(size(RealTime,1)), iDimT*n_pad,iBubble, R_Pad,RealTimePadOut)
          
		  R_Pad(:,1)= R_Pad(:,1) * dTaperSupportWindowN
		  R_Pad(:,2)= R_Pad(:,2) * dTaperSupportWindowN
		  R_Pad(:,3)= R_Pad(:,3) * dTaperSupportWindowN
		  ! Find volume acceleration d^2V/dt^2 [m^3/s^2] 
          ! This way the temporal derivative is calculated analytically so a spectral difference method is not needed.

          V_dd_pad = 4.0D0*pi*R_Pad(:,1)*(R_Pad(:,1)*R_Pad(:,3)+2.0D0*R_Pad(:,2)**2)

          ! Return to the initial dimensions!
          ! call INTERP1DFREQ(V_dd_Pad, V_dd_norm, cSpace,0)
          call INTERP1D( RealTimePadOut, V_dd_Pad, RealTime, V_dd_norm)

          ! V_dd_norm = 0.0D0;
          ! BubbleParams%R0(iBubble) = 1D-4

          ! Normalize by freq based on the order of the spectral derivative
          call LinScatterer(RealPressure , RealPressure, cSpace, 4.0/3.0*pi*BubbleParams%R0(iBubble)**3/cMediumParams%c0**2*cModelParams%freq0**2, 2 )
          ! call NonLinScatterer(RealPressure, RealPressure, cSpace, 1.0D-23*cModelParams%freq0**2, 2, 2 )
          tempcon(1+(i_loc-1)*iDimT:i_loc*iDimT) = cMediumParams.rho0*V_dd_norm + RealPressure;

          ! Multiply the volume acceleration with density so the units change to kg/s^2 , which do not include distance inside
          temploc(3*i_loc-2:3*i_loc) = [dBubbleLocationN(iBubble,1),dBubbleLocationN(iBubble,2),dBubbleLocationN(iBubble,3)] 
          BubbleID(i_loc) = iBubble
          i_loc = i_loc+1

       endif
    enddo

    DEALLOCATE(BubbleParams%P_driv,BubbleParams%T_driv)
    ! Deallocate to save memory
    if (ALLOCATED(LocPres)) then
       iMemAllocated=iMemAllocated - PRODUCT(SHAPE(LocPres)) * dpS
       DEALLOCATE(LocPres)
    endif
    iMemAllocated=iMemAllocated - PRODUCT(SHAPE(temp_globloc)) * dpS
    DEALLOCATE(temp_globloc)

    ! Create Bubble_count which contains the number of the bubbles inside each processors domain
    ! This is done to enhance speed . Other way is to save and load but it is totally inefficient
    Bubble_count = 0;
    call MPI_Allgather(i_loc-1,1,MPI_INTEGER,Bubble_count,1,MPI_INTEGER,MPI_COMM_WORLD,iErr)
    if ( sum(Bubble_count) .NE. BubbleParams%BubbleN) then
       write(acTemp,"('Error! Less bubbles that total are considered .')") iBubble, i
       call PrintToLog(acTemp,-1) 
       stop
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,iErr)

    !************************************* START OF TRANSFERRING DATA BETWEEN PROCESSORS ************************************************
    write(acTemp,"('Communicate Bubble Location and Pressure')");    call PrintToLog(acTemp,2)
    Bubble_disps=0; Bubble_locdisps=0; Bubble_condisps=0
    do i = 1,cSpace%cGrid%iProcN
       if (Bubble_count(i)>0)	then
          Bubble_disps(i) 		= sum(Bubble_count(1:i))		-	Bubble_count(i)			! Displacement variable in order to gather the value from every process
          Bubble_locdisps(i) 	= sum(3*Bubble_count(1:i))		-	3*Bubble_count(i)		! The same for location , now the variable is multiplied by 3 (x,y,z)
          Bubble_condisps(i) 	= sum(ibDimT*Bubble_count(1:i))	-	ibDimT*Bubble_count(i)	! The same for contrast source term , multiplied by length of T
       endif
    enddo

    ! allocate temporary buffer in order to change the order of scatterers based on the ID of the proc
    ALLOCATE(tempbuffer(i_loc-1))
    tempbuffer=BubbleID(1:i_loc-1)
    call MPI_Allgatherv(tempbuffer,i_loc-1,MPI_INTEGER8,BubbleID,Bubble_count,Bubble_disps,MPI_INTEGER8,MPI_COMM_WORLD,iErr)
    DEALLOCATE(tempbuffer)

    ! Same proceedure but in order to reduce time, do this instead
    dBubbleLocationN = dBubbleLocationN( BubbleID,:)	
    BubbleParams%R0 = BubbleParams%R0( BubbleID)

    ALLOCATE(dBubbleContrast(ibDimT*BubbleParams%BubbleN))
    iMemAllocated=iMemAllocated + PRODUCT(SHAPE(dBubbleContrast)) * dpS

    ! Same idea with previous. This is done in order each processor to know the contrast and location of all the scatterers
    dBubbleContrast=1D-30 ! Change below to MPI_DOUBLE_PRECISION if error
    call MPI_Allgatherv(tempcon,ibDimT*(i_loc-1),MPI_DOUBLE_PRECISION,dBubbleContrast,ibDimT*Bubble_count,Bubble_condisps,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,iErr)	
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); ! Wait for all the processors that the bubbles are located in to store the values

    !Wait until each processor has checked in , in order to finalize the number of bubbles in each processor and the proceed to next step
    iMemAllocated=iMemAllocated - ( PRODUCT(SHAPE(temploc)) + PRODUCT(SHAPE(tempcon)) )* dpS
    DEALLOCATE(temploc,tempcon)
    !****************************************************** CALCULATE PRESSURE FIELD **********************************************

    write(acTemp,"('Calculate the pressure field due to medium nonlinearity ')");    call PrintToLog(acTemp,2)

    ! Iniatilize pressure field with 0
    ! Iterate through each bubble and find the pressure field in the whole domain
    ! Add each field to the previous to find the total pressure distribution in the grid
    cSpace%cGrid%parD0 =  1.0D-30
    ! call NonlinContrastOperator_Ali(cSpace); ! Add Nonlinearities ! If beta=0 then it will result on 0 
	
    !First, transform to W-domain (use small transform, no wraparound regions)
    write(acTemp,"('Calculate pressure field due to the Bubble Cloud')");    call PrintToLog(acTemp,2)
    call CalculateCloudPressure(cSpace,i_neighbgp, NeighbGPIndex, dBubbleLocationN, dBubbleContrast)
	
    !*************************************************** END OF CALCULATE PRESSURE FIELD *******************************************    
    ! if (cSpace%cGrid%iProcID <=1) then

       ! Save each Bubble's Position 
       open (11, file=trim(trim(sOutputDir)//trim(BubbleParams%sBubbleDir)//'Bubble_Location'//int2str(cSpace%cGrid%iProcID)), status='REPLACE')
       do iBubble= 1 ,BubbleParams%BubbleN
          write(11, *) dBubbleLocationN(iBubble,1),dBubbleLocationN(iBubble,2),dBubbleLocationN(iBubble,3)
       end do
       close(11) 

       ! Save each Bubble's Contrast term for the correction of the "Own's Scattering"
       open (8,  file=trim(trim(sOutputDir)//trim(BubbleParams%sBubbleDir)//'Bubble_Contrast'//int2str(cSpace%cGrid%iProcID)), status='REPLACE')
       write(8,*) dBubbleContrast
       close(8) 

       ! For the polydisperse case , save each bubble's radius
       if (BubbleParams%Distribution == 'polydisperse') then
          open (8, file=trim(trim(sOutputDir)//'/Bubbles/'//trim('Bubble_Radius')//int2str(cSpace%cGrid%iProcID)), status='REPLACE')
          write(8,*) BubbleParams%R0
          close(8)
       endif
    ! endif
    ! This is for reducing memory
    open (8,  file=trim(trim(sOutputDir)//trim(BubbleParams%sBubbleDir)//'Bubble_iloc'//int2str(cSpace%cGrid%iProcID)), status='REPLACE')
    write(8,*) i_loc-1
    close(8) 

    iMemAllocated=iMemAllocated - ( PRODUCT(SHAPE(dBubbleContrast)) + PRODUCT(SHAPE(dBubbleLocationN)) ) * dpS - PRODUCT(SHAPE(BubbleID)) * i4bs
    DEALLOCATE(dBubbleContrast)
    DEALLOCATE(BubbleID)
    DEALLOCATE(dBubbleLocationN)

    call PrintToLog("End of Pressure Calculation",2) 
    call MPI_BARRIER(MPI_COMM_WORLD, iErr);

  END SUBROUTINE BubbleContrastOperator

  SUBROUTINE CreateBubbleCluster(cSpace,dBubbleLocationN, Domain_Range , cubes)

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
    !   The subroutine CreateBubbleCluster creates a cluster of microbubbles based
    !   on the minimum distance between the microbubbles and the domain that encloses
    !   all the bubbles ( maximum distance). This is done by dividing the cube into 
    !   smaller cubes and randomly position each bubble inside based on the two
    !	aforementioned dimensions
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace             io   type(space)  Space from which the data array needs
    !                                        to be stored.
    !   MinDistance		i    i8b	Minimum distance between microbubbles
    !   dBubbleLocationN	r    dp		Location of each microbubble (X,Y,Z)
    !   Domain_Range		r    dp         Dimensions of Domain (2,3) that encloses
    !						all the microbubbles
    !
    ! *****************************************************************************
    real(dp), intent(inout)					::		dBubbleLocationN(:,:)
    type(Space),  intent(in)				::		cSpace
    real(dp),intent(out)					::		Domain_Range(2,3)


    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   GridDist	r      Random number for the positioning of microbubbles
    !   i               i8b    Loop counter for each layer in Y + Z axis
    !   j               i8b    Loop counter for each layer in X axis
    !   Nx		i8b    Number of microbubbles in the X dimension of Domain
    !   Ny		i8b    Number of microbubbles in the Y dimension of Domain
    !   Nz		i8b    Number of microbubbles in the Z dimension of Domain
    !
    ! *****************************************************************************
    ! Cubes
    integer(i4b)         		::  	iErr
    real(dp),allocatable		::      GridDist(:,:)
    real(dp)					::      MaxDist(3), dLambdaMM, MinDist
    integer(i8b)				::		i ,j
    integer(i8b)				::		Nx, Ny, Nz
    logical 					::		cubes
    character(len=1024)			::      acTemp

    ! Lines

    integer(i8b)			::		arBubble_i(BubbleParams%BubbleN), arBubble_j(BubbleParams%BubbleN)
    integer(i8b)			:: 		iBubble, counter, diff,iBubbleInside
    integer(i8b)			::		LocXYZ(BubbleParams%BubbleN) , LocX(BubbleParams%BubbleN), LocY(BubbleParams%BubbleN), LocZ(BubbleParams%BubbleN) 
    real(dp)				::		randpos(BubbleParams%BubbleN)

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
    !   BubbleInit
    !   RANDOM_SEED
    !   RANDOM_NUMBER
    !   ALLOCATE/DEALLOCATE
    !	SaveExtraBubbleSlices
    !
    ! =============================================================================
    if (cSpace%cGrid%iProcID==0) then

       dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor 
       !unnormalized values of the resized Space
       Domain_Range(:,1) = (BubbleParams%ClusterDimsRatio(:,1) * (cSpace%iDimX-1) * cSpace.dDx + cSpace%iStartX*cSpace%dDx)* dLambdaMM  
       Domain_Range(:,2) = (BubbleParams%ClusterDimsRatio(:,2) * (cSpace%iDimY-1) * cSpace.dDx + cSpace%iStartY*cSpace%dDx)* dLambdaMM
       Domain_Range(:,3) = (BubbleParams%ClusterDimsRatio(:,3) * (cSpace%iDimZ-1) * cSpace.dDx + cSpace%iStartZ*cSpace%dDx)* dLambdaMM
       MinDist = BubbleParams%MinInBetweenDist*1e+3

		if (minval(Domain_Range(1,:) - (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaMM) < 0.0D0  .OR. &
		  maxval(Domain_Range(2,:) - (/cSpace.iDimX,cSpace%iDimY,cSpace%iDimZ/) * cSpace.dDx*dLambdaMM) >= 0.0D0 )then
		  write(acTemp,"('Error in CreateBubbleCluster, Domain out of boundaries')")
          call PrintToLog(acTemp,-1)
          stop
       endif

       if (cubes) then
          MaxDist = (  ( Domain_Range(2,1)-Domain_Range(1,1) )*( Domain_Range(2,2)-Domain_Range(1,2) )*( Domain_Range(2,3)-Domain_Range(1,3) )/( size(dBubbleLocationN,1) )  )**(1.0D0/3)   ! Maximum inbetween distance of bubbles in the specified Domain Range

          Nx = ceiling((Domain_Range(2,1)-Domain_Range(1,1))/MaxDist(1))
          Ny = ceiling((Domain_Range(2,2)-Domain_Range(1,2))/MaxDist(2))
          Nz = ceiling(size(dBubbleLocationN,1)*1.0/(Nx*Ny))

          if (Nx*Ny*Nz < size(dBubbleLocationN,1)) then
             write(acTemp,"('Error in CreateBubbleCluster, less discretized points than Bubble number')")
             call PrintToLog(acTemp,-1)
             stop
          endif

          ALLOCATE(GridDist(Nx*Ny*Nz,3))

          MaxDist = (Domain_Range(2,:)-Domain_Range(1,:))/(/ Nx,Ny,Nz /)    
          if (minval(MaxDist) < MinDist) then
             write(acTemp,"('Error ! Bubble Cluster size too small for ',I5,' Bubbles.')") size(dBubbleLocationN,1)
             call PrintToLog(acTemp,-1)
             write(acTemp,"('Max Bubble number for this Cluster size',I5,' Bubbles.')") size(dBubbleLocationN,1)/((MinDist/MaxDist(1))**(1.0D0/3))
             call PrintToLog(acTemp,-1)
             stop
          endif

          call RANDOM_SEED(PUT = [(1, i = 1,Nx*Ny*Nz*3)])
          call RANDOM_NUMBER(GridDist)

          do i = 0,Ny*Nz-1
             GridDist(1 + i*Nx : Nx + i*Nx,1) = Domain_Range(1,1) + 1.0D0*(1.0-GridDist(1 + i*Nx : Nx + i*Nx,1))*MinDist + &
                  1.0D0*MaxDist(1)*GridDist(1 + i*Nx : Nx + i*Nx,1) + MaxDist(1) * ([(j, j = 0,Nx-1)])
             GridDist(1 + i*Nx : Nx + i*Nx,2) = Domain_Range(1,2) + 1.0D0*(1.0-GridDist(1 + i*Nx : Nx + i*Nx,2))*MinDist + &
                  1.0D0*MaxDist(2)*GridDist(1 + i*Nx : Nx + i*Nx,2) + MaxDist(2) * mod(i,Ny)
             GridDist(1 + i*Nx : Nx + i*Nx,3) = Domain_Range(1,3) + 1.0D0*(1.0-GridDist(1 + i*Nx : Nx + i*Nx,3))*MinDist + &
                  1.0D0*MaxDist(3)*GridDist(1 + i*Nx : Nx + i*Nx,3) + MaxDist(3) * i/Ny
          enddo

          dBubbleLocationN = GridDist(1:size(dBubbleLocationN,1),:)

          DEALLOCATE(GridDist)

       else  ! Lines random type

          Nx = INT((Domain_Range(2,1) - Domain_Range(1,1))/MinDist,8) + 1
          Ny = INT((Domain_Range(2,2) - Domain_Range(1,2))/MinDist,8) + 1
          Nz = INT((Domain_Range(2,3) - Domain_Range(1,3))/MinDist,8) + 1

          if (Nx*Ny*Nz < BubbleParams%BubbleN) then
             write(acTemp,"('Too big minimum in-between distance or too small domain')")
             call PrintToLog(acTemp,-1)
             stop
          endif

          arBubble_i = 1
          LocXYZ = 0
          call RANDOM_SEED(PUT = [(1, i = 1,10*BubbleParams%BubbleN*3)]) !! Generate the same random number , when ran
          call RANDOM_SEED()	!GEnerate different random numbers for every run
          call RANDOM_NUMBER(randpos)
          arBubble_i = INT(randpos*Nx*Ny*Nz,8)
          call I_unista(arBubble_i, iBubble)
          LocXYZ(1:iBubble) = arBubble_i(1:iBubble)

          call RANDOM_SEED(PUT = [(1, i = 1,10*BubbleParams%BubbleN*3)])
          do while (iBubble < BubbleParams%BubbleN)
             counter=0; arBubble_i =0;diff=0;
             !			call RANDOM_SEED()	!GEnerate different random numbers for every run
             call RANDOM_NUMBER(HARVEST = randpos)
             diff = BubbleParams%BubbleN-iBubble
             arBubble_i(1:diff) = INT(randpos(1:diff)*Nx*Ny*Nz,8)
             LocXYZ(iBubble+1:iBubble+diff) = arBubble_i(1:diff)
             call I_unista(LocXYZ(1:iBubble+diff), counter)
             iBubble  = counter
          enddo

          LocX = mod(LocXYZ,Nx);
          arBubble_i  = (LocXYZ - LocX)/Nx;
          LocY = mod(  arBubble_i , Ny);
          arBubble_j  = ( arBubble_i - LocY)/Ny;
          LocZ = arBubble_j;

          dBubbleLocationN(:,1) = LocX*MinDist + Domain_Range(1,1);
          dBubbleLocationN(:,2) = LocY*MinDist + Domain_Range(1,2);
          dBubbleLocationN(:,3) = LocZ*MinDist + Domain_Range(1,3);
       endif

       dBubbleLocationN(:,1) = dBubbleLocationN(:,1)/dLambdaMM - cSpace%iStartX* cSpace%dDx
       dBubbleLocationN(:,2) = dBubbleLocationN(:,2)/dLambdaMM - cSpace%iStartY* cSpace%dDx
       dBubbleLocationN(:,3) = dBubbleLocationN(:,3)/dLambdaMM - cSpace%iStartZ* cSpace%dDx
       dBubbleLocationN = dBubbleLocationN + EPSILON(1.0D0) ! This is added because when the scatterer is exactly positioned at a gridpoint , the code is stucked

       if (minval( (/ minval(dBubbleLocationN,1),&
            (cSpace%iDimX-1) * cSpace.dDx - maxval(dBubbleLocationN(:,1),1), &
            (cSpace%iDimY-1) * cSpace.dDx - maxval(dBubbleLocationN(:,2),1), &
            (cSpace%iDimZ-1) * cSpace.dDx - maxval(dBubbleLocationN(:,3),1) /)) <= 0) then
          write(acTemp,"('Error in CreateBubbleCluster, Bubble positioned outside domain boundaries')")
          call PrintToLog(acTemp,-1)
          stop
       endif

    endif
	
    call MPI_BCAST(dBubbleLocationN,BubbleParams%BubbleN*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iErr)

  END SUBROUTINE CreateBubbleCluster

  SUBROUTINE CreateBubbleCluster_Random(cSpace,dBubbleLocationN, Domain_Range)

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
    !   The subroutine CreateBubbleCluster_Random creates a cluster of microbubbles based
    !   on the minimum distance between the microbubbles and the domain that encloses
    !   all the bubbles ( maximum distance). This is done by making a grid with the max
    !	distance in each dimension divided by the minimum mutual distance. Then generate
    !	randomly N points in this grid equally to the number of bubbles. Take the unique
    !	values from these numbers and if these are less than N , iterate until all the
    !	unique random numbers equal the number of microbubbles N.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace             io   type(space)  Space from which the data array needs
    !                                        to be stored.
    !   MinDistance		i    i8b	Minimum distance between microbubbles
    !   dBubbleLocationN	r    dp		Location of each microbubble (X,Y,Z)
    !
    ! *****************************************************************************
    real(dp), intent(inout)					::		dBubbleLocationN(:,:)
    type(Space),  intent(in)				::		cSpace
    real(dp),intent(out)					::		Domain_Range(2,3)


    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   Domain_Range		r    dp         Dimensions of Domain (2,3) that encloses all the microbubbles
    !   randpos				r    dp         vector with random values in the generated grid
    !   MinAllowableDist    r    dp		    Value of minimum mutual distance 
    !	dLambdaMM			r    dp	     	value of wavelength , normalization factor
    !   Nx,Ny,Nz			i8b    			Number of gridpoints in new generated domain in x, y and z
    !   iBubble,counter		i8b    			Index counter
    !   diff				i8b   			difference between N and the duplicate random numbers
    !   i,j					r    dp         temporary vectors in order for the proceedure to be clear
    !	LocX,LocY,LocZ		i				vectors with the location of bubbles in X,Y,Z in the generated domain with mindist
    !	LocXYZ				i				vector with one values for a triplet of (X,Y,Z)
    !
    ! *****************************************************************************
    integer(i4b)         	::  	iErr
    real(dp),parameter		::		MultFactor = 1.0D0
    integer(i8b)         	::  	Bubble_count(INT(MultFactor*BubbleParams%BubbleN,8)), BubbleID(INT(MultFactor*BubbleParams%BubbleN,8))
    integer				    ::		Bcounter, iBubble, BubbleIdx, i, MaxBubbles
    integer					::		Num_counter(cSpace%cGrid%iProcN),Bubble_Outliers(cSpace%cGrid%iProcN), Geom_Series(cSpace%cGrid%iProcN)
    real(dp)				::      dLambdaMM, MinAllowableDist,randpos(INT(BubbleParams%BubbleN*MultFactor,8),3), Dist , RadFactor(INT(BubbleParams%BubbleN*MultFactor))
    character(len=1024)		::      acTemp
    !Geometric series parameters 
    real(dp)				:: 		a, r
    integer(i8b)			::		n, S, iProc
    integer(i8b), allocatable 	:: 		GridDist(:)

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
    !	SaveExtraBubbleSlices
    !
    ! =============================================================================

    dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor 
    !unnormalized values of the resized Space
    Domain_Range(:,1) = (BubbleParams%ClusterDimsRatio(:,1) * (cSpace%iDimX-1) * cSpace.dDx + cSpace%iStartX*cSpace%dDx)* dLambdaMM  
    Domain_Range(:,2) = (BubbleParams%ClusterDimsRatio(:,2) * (cSpace%iDimY-1) * cSpace.dDx + cSpace%iStartY*cSpace%dDx)* dLambdaMM
    Domain_Range(:,3) = (BubbleParams%ClusterDimsRatio(:,3) * (cSpace%iDimZ-1) * cSpace.dDx + cSpace%iStartZ*cSpace%dDx)* dLambdaMM

    RadFactor(1:BubbleParams%BubbleN) = BubbleParams%MinInBetweenDist/BubbleParams%R0*1.0D+3

    if (minval(Domain_Range(1,:) - (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaMM) < 0.0D0  .OR. &
        maxval(Domain_Range(2,:) - (/cSpace.iDimX,cSpace%iDimY,cSpace%iDimZ/) * cSpace.dDx*dLambdaMM) >= 0.0D0 )then
       write(acTemp,"('Error in CreateBubbleCluster, Domain out of boundaries')")
       call PrintToLog(acTemp,-1)
       stop
    endif
	
    ! ================================================ Geometric Series ====================================================

    n = cSpace%cGrid%iProcN
    S = INT(MultFactor*BubbleParams%BubbleN,8)
    a = S*1.0D0/n
    r =  Solve_Geometric_Series(a,n,S)
    Geom_Series = 0
    Geom_Series(2:n) = INT((/ (a*r**iProc , iProc = n-2,0,-1) /),8)
    Geom_Series(1) = S- sum(Geom_Series)

    !=================================================== RANDOM LOCATION OF SCATTERERS =====================================
    ! call RANDOM_SEED(PUT = [(1, i = 1,10*BubbleParams%BubbleN*3)])  ! Random but the same at each run of INCS
    call RANDOM_SEED() ! Fully random , even at each run of INCS
    call RANDOM_NUMBER(randpos)
    randpos(:,1) = randpos(:,1)*(Domain_Range(2,1)-Domain_Range(1,1)) + Domain_Range(1,1)
    randpos(:,2) = randpos(:,2)*(Domain_Range(2,2)-Domain_Range(1,2)) + Domain_Range(1,2)
    randpos(:,3) = randpos(:,3)*(Domain_Range(2,3)-Domain_Range(1,3)) + Domain_Range(1,3)
    call MPI_BCAST(randpos,BubbleParams%BubbleN*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iErr)
	
    ! ============================== FIND HOW MANY ARE NOT REQUIREMENT THE HYPOTHESIS OF DISTANCE ==========================
    Bcounter=1
    MaxBubbles = sum(Geom_Series(1:cSpace%cGrid%iProcID))
    do iBubble =1,Geom_Series(cSpace%cGrid%iProcID+1)

       BubbleIdx = iBubble+MaxBubbles
       if( ANY(sqrt( (randpos(BubbleIdx,1) - randpos(1:BubbleIdx-1,1) )**2 + &
					 (randpos(BubbleIdx,2) - randpos(1:BubbleIdx-1,2) )**2 + &
					 (randpos(BubbleIdx,3) - randpos(1:BubbleIdx-1,3) )**2 ) - &
					 RadFactor(1:BubbleIdx-1)*BubbleParams%R0(1:BubbleIdx-1) <0 )) then

          Bubble_count(Bcounter)=BubbleIdx
          Bcounter=Bcounter+1
       endif

    enddo
    ! ================================== TRANSFER THEM FROM EACH PROCESSOR =======================================================
    call MPI_Allgather(Bcounter-1,1,MPI_INTEGER,Num_counter,1,MPI_INTEGER,MPI_COMM_WORLD,iErr) ! Array with the number of bubble outliers
    Bubble_outliers=0;

    do i = 1,cSpace%cGrid%iProcN
       if (Num_counter(i)>0)	then
          Bubble_outliers(i) 	= sum(Num_counter(1:i))  -	Num_counter(i)		! Displacement variable in order to gather the value from every process
       endif
    enddo

    ALLOCATE(GridDist(Bcounter-1))
    GridDist=Bubble_count(1:Bcounter-1)
	Bubble_count = 0
    call MPI_Allgatherv(GridDist,Bcounter-1,MPI_INTEGER8,Bubble_count,Num_counter,Bubble_outliers,MPI_INTEGER8,MPI_COMM_WORLD,iErr)
    DEALLOCATE(GridDist)
	
    ! ======================================= GENERATE THE LOCATION MATRIX FROM THE CORRECT VALUES ===================================
    randpos( (/ Bubble_count(1:sum(Num_counter)) /),1) = Domain_Range(1,1) - 1.0D0
    MaxBubbles = MAXLOC( (/ (i, i =1,sum(Num_counter) )/),1,Bubble_count(1:sum(Num_counter))-(/ (i, i =0,sum(Num_counter)-1) /) < BubbleParams%BubbleN)

    ! if there is no scatterer too close to the others, change the variables of interest and proceed without any more actions
    if(sum(Num_counter)==0) then
       Bubble_count(MaxBubbles)=BubbleParams%BubbleN
       dBubbleLocationN(1:BubbleParams%BubbleN,:) = randpos( 1:BubbleParams%BubbleN,:)
       BubbleParams%R0(1:BubbleParams%BubbleN)= BubbleParams%R0(1:BubbleParams%BubbleN)
    ! if there a number of scatterers less than the number of interest that are too close and 
    elseif(Bubble_count(MaxBubbles)<BubbleParams%BubbleN .AND. MaxBubbles==sum(Num_counter)) then
       BubbleID(1:S-MaxBubbles) = (/ pack( (/ (i , i = 1,S ) /), randpos(:,1) >= Domain_Range(1,1) )		/)
       dBubbleLocationN(1:S-MaxBubbles,:) = randpos( (/ BubbleID(1:S-MaxBubbles) /) ,:)
       BubbleParams%R0(1:S)= BubbleParams%R0( (/ BubbleID(1:S-MaxBubbles) , INT(pack( (/ (i , i = 1,S ) /), randpos(:,1) < Domain_Range(1,1) ),8) /) )
    else
       MaxBubbles= MaxBubbles+1
       randpos(BubbleParams%BubbleN+MaxBubbles:S,1) = Domain_Range(1,1)-1
       BubbleID(1:BubbleParams%BubbleN) = (/ pack( (/ (i , i = 1,Bubble_count(MaxBubbles) ) /), randpos(1:Bubble_count(MaxBubbles),1) >= Domain_Range(1,1) ) /)
       dBubbleLocationN(1:BubbleParams%BubbleN,:) = randpos( (/ BubbleID(1:BubbleParams%BubbleN) /),:)
       BubbleParams%R0= BubbleParams%R0( (/ BubbleID(1:BubbleParams%BubbleN) /) )
    endif

    ! ======================================= FILL THE MATRIX WITH THE REST CORRECT VALUES ==========================================

    iBubble=1
    i=0
    if (Bubble_count(MaxBubbles)<BubbleParams%BubbleN .AND. cSpace%cGrid%iProcID==0) then
       BubbleIdx = S-MaxBubbles
       do while (iBubble < MaxBubbles+1)
          if (i==0) then
             call RANDOM_SEED()
             ! call RANDOM_SEED(PUT = [(1, i = 1,10*BubbleParams%BubbleN*3)])
             call RANDOM_NUMBER(randpos)
             randpos(:,1) = randpos(:,1)*(Domain_Range(2,1)-Domain_Range(1,1)) + Domain_Range(1,1)
             randpos(:,2) = randpos(:,2)*(Domain_Range(2,2)-Domain_Range(1,2)) + Domain_Range(1,2)
             randpos(:,3) = randpos(:,3)*(Domain_Range(2,3)-Domain_Range(1,3)) + Domain_Range(1,3)
             i=1
          elseif (i==S-1) then
             i=0
          else
             BubbleIdx = S-MaxBubbles+iBubble
             if(ALL(sqrt( (randpos(i,1) - dBubbleLocationN(1:BubbleIdx-1,1) )**2    + &
						  (randpos(i,2) - dBubbleLocationN(1:BubbleIdx-1,2) )**2    + &
						  (randpos(i,3) - dBubbleLocationN(1:BubbleIdx-1,3) )**2  ) - &
						  RadFactor(1:BubbleIdx-1)*BubbleParams%R0(1:BubbleIdx-1) >0)) then
						
                dBubbleLocationN(BubbleIdx,:)  = randpos(i,:)
                iBubble = iBubble+1
             endif
             i=i+1
          endif

       enddo
    endif
    ! ================================================================================================================================
    call MPI_BCAST(dBubbleLocationN,BubbleParams%BubbleN*3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iErr)
    call MPI_BCAST(BubbleParams%R0,BubbleParams%BubbleN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,iErr)
    ! dBubbleLocationN(1,:)=  (/0.576333333333333 ,   0.0823333333333333 ,    5.187/) ! At a gridpoint 
    ! dBubbleLocationN(1,:)=  (/0.576333333333333  ,   0.0823333333333333  ,    5.187/)  ! Middle of gridpoints 
    ! dBubbleLocationN(2,:)=  (/-0.411666666666667 ,   0.082333333333333  ,    6.25733333333333  /)

    dBubbleLocationN(:,1) = dBubbleLocationN(:,1)/dLambdaMM - cSpace%iStartX* cSpace%dDx
    dBubbleLocationN(:,2) = dBubbleLocationN(:,2)/dLambdaMM - cSpace%iStartY* cSpace%dDx
    dBubbleLocationN(:,3) = dBubbleLocationN(:,3)/dLambdaMM - cSpace%iStartZ* cSpace%dDx
    dBubbleLocationN = dBubbleLocationN + EPSILON(1.0D0) ! This is added because when the scatterer is exactly positioned at a gridpoint , the code is stucked

    if (minval( (/ minval(dBubbleLocationN,1),&
         (cSpace%iDimX-1) * cSpace.dDx - maxval(dBubbleLocationN(:,1),1), &
         (cSpace%iDimY-1) * cSpace.dDx - maxval(dBubbleLocationN(:,2),1), &
         (cSpace%iDimZ-1) * cSpace.dDx - maxval(dBubbleLocationN(:,3),1) /)) <= 0) then
       write(acTemp,"('Error in CreateBubbleCluster, Bubble positioned outside domain boundaries')")
       call PrintToLog(acTemp,-1)
       stop
    endif

  END SUBROUTINE CreateBubbleCluster_Random

  SUBROUTINE ProcBubbleCount(cSpace, dBubbleLocationN, tempcon, iProc_NumBubbles)

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
    ! 	each cpu used).Besides this , the values of the time signature of each scatterer
    !	is loaded . Then, only the values of the bubbles inside each subdomain are kept.
    !	The rest are discarded.
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

    real(dp), intent(in)					::		dBubbleLocationN(:,:)
    real(dp),allocatable, intent(inout)		::		tempcon(:)
    type(Space),  intent(in)				::		cSpace
    integer, intent(inout)					::		iProc_NumBubbles


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
    !	dBubbleCon		dp	   vector that is used for transferring the time signature array
    !
    ! *****************************************************************************

    integer(i4b)         	::  	iErr , BubbleID(size(dBubbleLocationN,1))
    integer(i8b)			::		iBubble ,XYZIndex,  iBeamIndex(3),iDimT
    real(dp),allocatable   	::		dBubbleCon(:)
    character(len=1024)		::      acTemp   


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
    !02	  load the saved time signatures 
    !03   Initialize the number of bubbles inside each cpu domain 
    !04	  Count for each scatterer inside the domain
    !04	  This is done by generating one value for each triplet of (X,Y,Z) 
    !04	  and checking if this value is between 0 and max number of gridpoints
    !05	  if 2nd iter or higher, transfer the time signatures of the bubbles 
    !05	  located in each cpu domain.
    !
    ! =============================================================================

    !    !Calculate Number of bubbles below half wavelength of Fnyq
    !    integer(i8b)			::		count5(BubbleParams%BubbleN/cSpace%cGrid%iProcN),count6(size(count5,1)),count9(size(count5,1)),count10(size(count5,1)), Bubble_Loc
    !    integer(i8b)			::		count51(BubbleParams%BubbleN/cSpace%cGrid%iProcN),count61(size(count51,1)),count91(size(count51,1)),count101(size(count51,1)) , length(4)
    !    real(dp)		     	:: 		diff(BubbleParams%BubbleN),dLambdaMM
    !    
    !    dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor 
    !	count5=0; count6=0;count9=0;count10=0;
    !	count51=0; count61=0;count91=0;count101=0;
    !	do iBubble =  1,BubbleParams%BubbleN/cSpace%cGrid%iProcN
    !		diff=0;
    !	   Bubble_Loc = iBubble+BubbleParams%BubbleN/cSpace%cGrid%iProcN*cSpace%cGrid%iProcID

    !       diff = sqrt((dBubbleLocationN(Bubble_Loc,1) - dBubbleLocationN(:,1))**2 +(dBubbleLocationN(Bubble_Loc,2) - dBubbleLocationN(:,2))**2+ (dBubbleLocationN(Bubble_Loc,3) - dBubbleLocationN(:,3))**2)
    !       count51(iBubble) = count(diff *dLambdaMM*1e-3 < 1482.0/2/9E+6*1) -1
    !       count61(iBubble) = count(diff *dLambdaMM*1e-3 < 1482.0/2/9E+6*2) -1
    !       count91(iBubble) = count(diff *dLambdaMM*1e-3 < 1482.0/2/9E+6*3) -1
    !       count101(iBubble) = count(diff *dLambdaMM*1e-3 < 1482.0/2/9E+6*4) -1

    !       count5(iBubble) = count(diff(Bubble_Loc+1:BubbleParams%BubbleN) *dLambdaMM*1e-3 < 1482.0/2/9E+6*1)
    !       count6(iBubble) = count(diff(Bubble_Loc+1:BubbleParams%BubbleN) *dLambdaMM*1e-3 < 1482.0/2/9E+6*2) 
    !       count9(iBubble) = count(diff(Bubble_Loc+1:BubbleParams%BubbleN) *dLambdaMM*1e-3 < 1482.0/2/9E+6*3) 
    !       count10(iBubble) = count(diff(Bubble_Loc+1:BubbleParams%BubbleN) *dLambdaMM*1e-3 < 1482.0/2/9E+6*4)
    !     enddo
    !   		
    !	length = INT( log10( real( (/ sum(count5(:)) , sum(count6(:)),sum(count9(:)),sum(count10(:)) /),dp) ),8)+3
    !	
    !	write(*,'(I<length(1)>,I<length(2)>,I<length(3)>,I<length(4)>)',advance='no') sum(count5(:)) , sum(count6(:)),sum(count9(:)),sum(count10(:))
    !	
    !	length = INT( log10( (/ sum(count51(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN, sum(count61(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN,sum(count91(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN,sum(count101(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN /) ),8)+3

    !	write(*,'(F<length(1)+8>.5, F<length(2)+8>.5, F<length(3)+8>.5, F<length(4)+8>.5)',advance='no') sum(count51(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN, sum(count61(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN,sum(count91(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN,sum(count101(:))*1.0/BubbleParams%BubbleN*cSpace%cGrid%iProcN
    !	write(*,'(4F)')

    !    call MPI_Barrier(MPI_COMM_WORLD,iErr)
    !***********************************CALCULATE HOW MANY MICROBUBBLES ARE LOCATED IN EACH PROCESSOR***********************************************	

    iDimT = cSpace%iDimT
    if (iProc_NumBubbles>0) then
       ALLOCATE(dBubbleCon(iDimT*BubbleParams%BubbleN))
       open(9,  file=trim(trim(sOutputDir)//trim(BubbleParams%sBubbleDir)//trim('Bubble_Contrast')//int2str(0)), status='OLD')
       read(9,*) dBubbleCon
       close(9)
    endif

    iProc_NumBubbles=0
    do iBubble =  1,BubbleParams%BubbleN
       ! Determine Beam Index of the point source
       ! In ParnacDataDef it is specified that cSpace%iStartT/X/Y/Z should be included 
       ! The difference here is that the dBubbleLocationN parameter is in the normalized space
       ! So if we change this in real space then 
       ! dBubbleLocation = (dBubbleLocationN + (/cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ/)*cSpace%dDx*dLambdaNN
       ! So there is no difference if the starting position is not added in all the equations
       iBeamIndex(3) = nint(dBubbleLocationN(iBubble,3)/cSpace%dDx) 
       iBeamIndex(1) = nint(dBubbleLocationN(iBubble,1)/cSpace%dDx) - iBeamOffsetX(cSpace,iBeamIndex(3))
       iBeamIndex(2) = nint(dBubbleLocationN(iBubble,2)/cSpace%dDx) - iBeamOffsetY(cSpace,iBeamIndex(3))

       !Find the Index in the cSpace%cGrid%aiD0Loc of the BubbleLocation
       XYZIndex = iBeamIndex(3)*cSpace%iDimY*cSpace%iDimX+iBeamIndex(2)*cSpace%iDimX+iBeamIndex(1)-(cSpace%cGrid%iProcID) * cSpace%cGrid%iD0LocN
       if (XYZIndex>=0 .AND. XYZINdex<cSpace%cGrid%iD0LocN) then
          iProc_NumBubbles = iProc_NumBubbles +1
          if (ALLOCATED(tempcon)) tempcon(1+(iProc_NumBubbles-1)*iDimT:iProc_NumBubbles*iDimT) = dBubbleCon(1+(iBubble-1)*iDimT:iBubble*iDimT)
       endif

    enddo

    if (ALLOCATED(dBubbleCon)) DEALLOCATE(dBubbleCon)
  END SUBROUTINE ProcBubbleCount

  SUBROUTINE CalculateCloudPressure(cSpace,i_neighbgp, NeighbGPIndex, dBubbleLocationN, dBubbleContrast)

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
    !	cloud.  This is done in two phases. The first phase is for each iteration ,except
    !   the last one. For this phase, all the neighbouring point of microbubbles are used
    !	and then the pressure is calculated in these points by all the microbubbles.
    !	In order to keep memory low , this is done by each cpu using the neighbouring points
    !	that are located inside its subdomain. The second phase is the last iteration.
    !	At this iteration ,the pressure is calculated at each point of the domain.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace              	 io   type(space)  The space in which the point has to be given		
    !   dBubbleLocationN         r   dp   a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
    !   dBubbleContrast	         r   dp   a vector with the time signatures of the microbubbles
    !   NeighbGPIndex	         i   i8b  a vector with the location of the neighbouring points of microbubbles in each cpu domain
    !   i_neighbgp	         	 i   i8b  number of neighbouring points ( Do not get confused, not all elements of NeighbGPIndex
    !									  have the values , only 1:i_neighbgp)
    !   last_iter		         l        logical value which states if it is the last iteration or not
    !
    ! *****************************************************************************
    real(dp), intent(inout) 	:: 		dBubbleLocationN(:,:)
    real(dp), intent(inout)		:: 		dBubbleContrast(:) 
    integer(i8b),intent(in)		::		i_neighbgp
    integer(i8b),intent(in)		::		NeighbGPIndex(:)
    type(Space), intent(inout) 	:: 		cSpace

    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   P_Scattered       		dpc   Scattered Pressure Matrix for MKL
    !   dSincMult     			dpc   a matrix that contains spatially filtered dirac
    !   NeighbLocInd			i8b   a vector which contains the indexes of the neighbouring points
    !   dGlobBXYZ		        r     a matrix of the gridpoints location of each cpu domain (same dim as dBubbleLocationN)
    !   M, ibDimT	       		i     integer value for MKL , Time length
    !   N, GridPointsNum    	i	  integer value for MKL , Space length 
    !   K, BLayer				i 	  integer value for MKL , Bubbles per layer
    !   iIndex, i, iBubble	    i8b   index variables
    !   B_Div			   		i8b   divisor of total numbers of variables 
    !	dLambdaMM				dp	  value of wavelength , normalization factor
    !	dK_c					dp	  Kutoff frequency 
    !
    ! *****************************************************************************


    real(dp),allocatable			::  P_scattered(:,:), dSincMult(:,:), P_scattered_Total(:,:), dGlobBXYZ(:,:) 
    integer(i8b), allocatable		::  NeighbLocInd(:), XYZIndex(:)
    integer, allocatable 			::	GridPoint_OnProc(:,:)

    integer					::	M, N, K, iDimT, BLayer, iDimX, iDimY, iDimZ, iDimXYZ, i,j, BubbleProc
    integer					::	GridPoint_Count(cSpace%cGrid%iProcN), GridPoint_Disps(cSpace%cGrid%iProcN), TimePoints(cSpace%iDimT)
    integer 	         	::  iErr, iStat(MPI_STATUS_SIZE), iREQUEST
    integer(i8b)         	::  iIndex, iLI , x, y, z, iBubble, B_Div , GridPointsNum
    real(dp)             	::  dLambdaMM, dK_c
	
	real(dp)				:: dRightBand(BubbleParams%BubbleN,3), dLeftBand(BubbleParams%BubbleN,3), dBand(BubbleParams%BubbleN,3)
    real(dp)				:: dTaperingWindowX(cSpace%iDimx,BubbleParams%BubbleN) ,dTaperingWindowY(cSpace%iDimY,BubbleParams%BubbleN) ,dTaperingWindowZ(cSpace%iDimZ,BubbleParams%BubbleN)
	character(len = 1024)	:: actemp 

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
    !01	  and set the points for calculation only to neighbouring, otherwise calculate for every point of domain
    !02	  Set a good choice for the number of Bubbles per Layer ( should check memory usage)
    !03	  Iterate for the number of Layers (Total Bubbles / No. Bubbles per Layer)
    !04	  Calculate dirac function matrix
    !05	  Do the matrix multiplication based on the space points for calculation     
    !06	  if not in last iteration , transfer the calculated values to a matrix for easy replacement to cSpace var
    !07	  else just replace the calculated values
    !
    ! =============================================================================

    iDimT    = cSpace%iDimT
    iDimX  	 = cSpace%iDimX
    iDimY  	 = cSpace%iDimY
    iDimZ  	 = cSpace%iDimZ
							 
    dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor λ in [mm]
    dK_c = pi /cSpace%dDx    

    ! At this point, the normalized location of the gridpoints of interest are calculated
    if (BubbleParams%GridPointsPressure == 'all') then  
       ! In this case, every gridpoint's location is computed
       write(acTemp,"('Calculate the contrast source term at each spatial gridpoint')"); call PrintToLog(acTemp,3)
       GridPointsNum = cSpace%cGrid%iD0LocN
       ALLOCATE(dGlobBXYZ(GridPointsNum,3))
       ! This is the way to compute them and in this order. This is important because later they will be used
       ! for the computation of the total scattered pressure and they should be placed in this order
       do i = 0, GridPointsNum-1
          dGlobBXYZ(1+i,1) = (cSpace%cGrid%aiD0Loc(i+1,1) + iBeamOffsetX(cSpace,cSpace%cGrid%aiD0Loc(i+1,3))) * cSpace.dDx
          dGlobBXYZ(1+i,2) = (cSpace%cGrid%aiD0Loc(i+1,2) + iBeamOffsetY(cSpace,cSpace%cGrid%aiD0Loc(i+1,3))) * cSpace.dDx
          dGlobBXYZ(1+i,3) = (cSpace%cGrid%aiD0Loc(i+1,3) ) * cSpace.dDx  
       enddo
       BLayer=2000; 
    else 
       ! In this case, only the gridpoints that are added in the main module above, are computed
       ! The number of the neighbouring points for every scatterer equals to size(p,1)**3
       write(acTemp,"('Calclulate the contrast source term in the gridpoints closer to the scatterers')"); call PrintToLog(acTemp,3)
	   
       GridPointsNum = i_neighbgp/cSpace%cGrid%iProcN   ! Do this for all processors 
	   
	   ! But for the last process compute the number of gridpoints left, which will not be always equal to the previous value.
       if (cSpace%cGrid%iProcID==cSpace%cGrid%iProcN-1) GridPointsNum=i_neighbgp-(cSpace%cGrid%iProcN-1)*(i_neighbgp/cSpace%cGrid%iProcN)
	   
       ALLOCATE(dGlobBXYZ(GridPointsNum,3))
	   ! Because all cpus see all the gridpoints, we have to use a way to separate those
	   ! that are located in a specific processor. Generate a matrix, where first column is the global value of each gridpoint
	   ! and the 2nd value is the processor, where it is located.
       ALLOCATE(GridPoint_OnProc(2,GridPointsNum))
	   
       do i = 0,GridPointsNum-1 
          ! Divide the number of gridpoints based on the number of cpus
		  j  = NeighbGPIndex( (i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID+i+1);
		  
          GridPoint_OnProc(1,i+1) = i
          GridPoint_OnProc(2,i+1) = j/cSpace%cGrid%iD0LocN
		  
		  ! This is the way to compute the normalized global values of each gridpoint  in space.
          dGlobBXYZ(1+i, 1)        = mod(j, cSpace%iDimX)*cSpace%dDx;;
          j         = (j - mod(j, cSpace%iDimX)) / cSpace%iDimX;
          dGlobBXYZ(1+i, 2)        = mod(j, cSpace%iDimY)*cSpace%dDx;
          j        = (j - mod(j, cSpace%iDimY)) / cSpace%iDimY;
          dGlobBXYZ(1+i, 3)        = j*cSpace%dDx;
       end do
       ! Increase the number of scatterers per layer, because there are less points than final iteration
       BLayer=1E4;
    endif
    M = iDimT ; N = GridPointsNum; K =BLayer

    ! Because of big size arrays for clusters higher than 10^5 bubbles
    ! due to memory , a big array can not be created ,so I implement the idea
    ! of layering in the matrix, So I do a loop for every K bubbles and calculate
    ! the scattered pressure. I do this until all the bubbles hve been included 
    if (BubbleParams%BubbleN <= K) K = BubbleParams%BubbleN
    B_Div = BubbleParams%BubbleN/K;

    ALLOCATE(P_scattered(M,N))
    ALLOCATE(P_scattered_Total(M+1,N))
    ALLOCATE(dSincMult(K,N))

    iMemAllocated=iMemAllocated + ( PRODUCT(SHAPE(dSincMult)) + PRODUCT(SHAPE(P_scattered)) + PRODUCT(SHAPE(P_scattered_Total)) )* dpS + PRODUCT(SHAPE(dGlobBXYZ)) * i4bs
    write(acTemp,"('Create ', I3,' Layer(s) of ', I<INT(log10(K*1.0D0))+1>, ' Bubble(s) for Dirac Array')") B_Div,K;    call PrintToLog(acTemp,3)

    ! Here we generate only the contrast source term of interest! This array will include the nonlinear contrast source term for a 
    ! number of time steps equal to (cSpace%iDimT+1)/iProcN in every processor. This means that instead of dividing the domain based on
    ! the number of cores, the time is divided. In this way we have the same number of computations per processor and every processor 
    ! takes part in the computations
    P_scattered_Total = 1.0D-30;
	
    do iBubble = 1,B_Div ! This is the loop for the number of layers of the bubbles 
       P_scattered = 1.0D-30
       dSincMult   = 1.0D-30        
       ! Calculate the spatial signature in the point source definition for every scatterer
       do i = 0, N-1  					  				
          dSincMult(:,1+i) =	 dSinc(dK_c*(dBubbleLocationN(K*(iBubble-1)+1:K*iBubble,1) - dGlobBXYZ(1+i,1) ))  &
								*dSinc(dK_c*(dBubbleLocationN(K*(iBubble-1)+1:K*iBubble,2) - dGlobBXYZ(1+i,2) ))  &
								*dSinc(dK_c*(dBubbleLocationN(K*(iBubble-1)+1:K*iBubble,3) - dGlobBXYZ(1+i,3) )) 
       enddo
       ! DO the multiplication with the temporal factor , in order to calculte the scattered pressure
       ! Each processor computes the scattered pressure in the gridpoints of interest only for the iterated layer of microbubbles
       CALL DGEMM('N','N',M,N,K,1.0D0,RESHAPE(dBubbleContrast(1+K*(iBubble-1)*M:K*iBubble*M),(/M,K/)), M,dSincMult ,K,0.0D0,P_scattered,M)	 
       P_scattered_Total(2:M+1,:)= P_scattered_Total(2:M+1,:) + P_Scattered
    enddo
	! We need this row in order to know for which value of the gridpoint, the pressure of each row(gridpoint) is computed so we can send it to the right processor
    P_scattered_Total(1,:) = NeighbGPIndex(1+(i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID:(i_neighbgp/cSpace%cGrid%iProcN)*cSpace%cGrid%iProcID+N)*1.0D0

    iMemAllocated=iMemAllocated - ( PRODUCT(SHAPE(dSincMult)) + PRODUCT(SHAPE(P_scattered)) )* dpS - PRODUCT(SHAPE(dGlobBXYZ)) * i4bs
    DEALLOCATE(P_Scattered)	
    DEALLOCATE(dSincMult) 
    DEALLOCATE(dGlobBXYZ)

    write(acTemp,"('Matrix Multiplication has finished.')");    call PrintToLog(acTemp,3)
    write(acTemp,"('GridPoints Distribution to each processor starts.')");    call PrintToLog(acTemp,3)
    if (BubbleParams%GridPointsPressure == 'all') then
       cSpace%cGrid%parD0 = cSpace%cGrid%parD0 + RESHAPE(P_scattered_Total(2:M+1,:)*(dK_c/pi)**3.0 *(cModelParams%freq0/cMediumParams%c0),(/M*N/))
    else
	   ! Here, the array for the number of elements GridPoint_Count and 
	   ! the array for the starting index of the location of the elements per processor GridPoint_Disps is generated.
       GridPoint_Count=0;
	   GridPoint_Disps = 0;
	   
       do BubbleProc = 0,cSpace%cGrid%iProcN-1
          j = i_neighbgp/cSpace%cGrid%iProcN
          i = j; 
		  if (BubbleProc == cSpace%cGrid%iProcN-1) i = i_neighbgp - j*(cSpace%cGrid%iProcN-1)
		  ! Count how many gridpoints are located in the working processor from all the processors in order to be able to communicate them later 
          GridPoint_Count(BubbleProc+1) = COUNT( NeighbGPIndex(j*BubbleProc+1:j*BubbleProc+i)/cSpace%cGrid%iD0LocN == cSpace%cGrid%iProcID)
		  
		  ! Compute the index location for the total number of gridpoints for each processor 
          if (GridPoint_Count(BubbleProc+1)>0)	then
             GridPoint_Disps(BubbleProc+1) 	= sum((M+1)*GridPoint_Count(1:BubbleProc+1))	-  (M+1)*GridPoint_Count(BubbleProc+1)		
          endif
       enddo

       N = sum(GridPoint_Count);
       ALLOCATE(P_Scattered((M+1)*N,1))
       iMemAllocated=iMemAllocated + PRODUCT(SHAPE(P_Scattered))* dpS 	
	   P_scattered = 0.0D0
       
	   do BubbleProc = 0,cSpace%cGrid%iProcN-1
		  ! Count how many points on the running processor are located at the iterated value of processor
          j = count( GridPoint_OnProc(2,:) == BubbleProc) 
          ALLOCATE(XYZIndex(j))
          XYZIndex = (/ pack( GridPoint_OnProc(1,:)+1, GridPoint_OnProc(2,:) == BubbleProc ) /) ! This has the value of the gridpoints relevant to the iterated value of processor
		  ! Send then the relevant values from each processor to the iterated value of processor.
		  ! In this way, each processor will get the pressure field that was computed only on the gridpoints that are located 
          call MPI_Gatherv(P_scattered_Total(:,XYZIndex),(M+1)*j,MPI_DOUBLE_PRECISION,P_scattered,(M+1)*GridPoint_Count,GridPoint_Disps,MPI_DOUBLE_PRECISION,BubbleProc,MPI_COMM_WORLD,iErr)	
          DEALLOCATE(XYZIndex)
       enddo
	   
       P_scattered_Total =  RESHAPE(P_scattered,(/M+1,N/))
       iMemAllocated=iMemAllocated -  PRODUCT(SHAPE(P_Scattered))* dpS 	
       DEALLOCATE(P_Scattered)	
       DEALLOCATE(GridPoint_OnProc)

       if (N>0) then
          ALLOCATE(NeighbLocInd(M*N))
          TimePoints  = (/(i , i = 1,M)/);  NeighbLocInd = (/ (INT(P_scattered_Total(1,j)- cSpace%cGrid%iProcID * cSpace%cGrid%iD0LocN)*M+TimePoints , j =1,N) /)
          cSpace%cGrid%parD0(NeighbLocInd) = cSpace%cGrid%parD0(NeighbLocInd)+ RESHAPE(P_scattered_Total(2:M+1,:)*(dK_c/pi)**3.0 *(cModelParams%freq0/cMediumParams%c0),(/M*N/))	
          DEALLOCATE(NeighbLocInd)
       endif
    endif
    iMemAllocated=iMemAllocated -  PRODUCT(SHAPE(P_scattered_Total))* dpS 	
    DEALLOCATE(P_scattered_Total)
    write(acTemp,"('Proceedure finished.')");    call PrintToLog(acTemp,3)

  END SUBROUTINE CalculateCloudPressure

  SUBROUTINE LinScatterer(InputVal , OutputVal, cSpace, Magnitude, OrderFFT )

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
    real(dp), intent(in) 		:: 		InputVal(:)
    real(dp), intent(out) 		:: 		OutputVal(:)
    integer(i8b),intent(in)		::		OrderFFT
    real(dp), intent(in) 		:: 		Magnitude
    type(Space), intent(in) 	:: 		cSpace


    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   InputValW       		dpc   Frequency components of input signal
    !   OutputValW     			dpc   Frequency components of output signal
    !   InputLen   				i8b   Length of Input Signal
    !   EvenOrOdd       		i8b   Variable, 0 if even or 1 if odd 
    !   iNumPlanTransform     	i8b   plan for FFT
    !   iNumPlanTransform_inv	i8b   plan for IFFT
    !   dTaperMaxFreqWindow     dp    Tapering of Max Frequency to remove artifact
    !   dLeftBand     			dp    Left limit of banded tapering
    !   dRightBand		   		dp    Right limit of banded tapering
    !

    complex(dpc) 				::	 		InputValW(size(InputVal)/2+1)
    complex(dpc),allocatable	::			OutputValW(:)
    complex(dpc),allocatable	::			dMultFactor(:)
    integer(i8b) 				::			InputLen, OutputLen, EvenOrOdd, iNumPlanTransform, iNumPlanTransform_inv, iDimW

    !Filtering and windowing parameters 
    real(dp) 	::			dTaperMaxFreqWindow(size(InputVal)/2+1), dLeftBand,dRightBand, dDomega
    integer		::			i,iErr

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
    InputLen  = size(InputVal)
    EvenOrOdd = mod(InputLen,2_i8b)
    iDimW = InputLen/2+1

    ! Do a FFT of InputLen-point (Initial Variable Length) with complex numbers 
    call dfftw_plan_dft_r2c_1d(iNumPlanTransform, InputLen, InputVal, InputValW, FFTW_ESTIMATE )
    call dfftw_execute_dft_r2c(iNumPlanTransform, InputVal, InputValW)
    call dfftw_destroy_plan(iNumPlanTransform)

    !Allocate to upsample or decimate . OutputValW has to have more or equal points to the initial signal 
    ALLOCATE(OutputValW(max(iDimW,OutputLen/2)))

    ! Tapering of the max frequency because of artifact created due to higher frequencies
    dLeftBand = 0.0
    dRightBand =0.1_dp  
    dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(InputLen*cSpace%dDt),dLeftBand,dRightBand) 

    ALLOCATE(dMultFactor(iDimW))
    dMultFactor = 1.0D0
    if (OrderFFT /= 0) then
       dDOmega = two_pi/real(InputLen * cSpace%dDt,dp) 
       dMultFactor(1:iDimW) = ( (/ (i,i=0,iDimW-1) /) * dDOmega * im)**OrderFFT
    endif

    ! Zero-padding after the values which will be used for interpolation
    ! This is done by initializing OutputValW with 0 ( Really small number due to precision)
    ! Replace the first half part with the values of initial variable.
    OutputValW = 0.0D0
    OutputValW(1:iDimW) = InputValW(1:iDimW) * dMultFactor(1:iDimW) * dTaperMaxFreqWindow /InputLen

    ! Inverse FFT of OutputLen-point	
    call dfftw_plan_dft_c2r_1d(iNumPlanTransform_inv, OutputLen, OutputValW, OutputVal, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(iNumPlanTransform_inv, OutputValW, OutputVal)
    call dfftw_destroy_plan(iNumPlanTransform_inv)

    OutputVal = OutputVal * Magnitude * dTaperingWindow(InputLen,cSpace%dDt,2.0_dp, 2.0_dp) ;

    DEALLOCATE(OutputValW)
    DEALLOCATE(dMultFactor)
  END SUBROUTINE LinScatterer

  SUBROUTINE NonLinScatterer(InputVal , OutputVal, cSpace, Magnitude, OrderFFT, OrderP )

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
    real(dp), intent(in) 		:: 		InputVal(:)
    real(dp), intent(out) 		:: 		OutputVal(:)
    integer(i8b),intent(in)		::		OrderFFT, OrderP
    real(dp), intent(in) 		:: 		Magnitude
    type(Space), intent(in) 	:: 		cSpace


    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   InputValW       		dpc   Frequency components of input signal
    !   OutputValW     			dpc   Frequency components of output signal
    !   InputLen   				i8b   Length of Input Signal
    !   EvenOrOdd       		i8b   Variable, 0 if even or 1 if odd 
    !   iNumPlanTransform     	i8b   plan for FFT
    !   iNumPlanTransform_inv	i8b   plan for IFFT
    !   dTaperMaxFreqWindow     dp    Tapering of Max Frequency to remove artifact
    !   dLeftBand     			dp    Left limit of banded tapering
    !   dRightBand		   		dp    Right limit of banded tapering
    !

    complex(dpc) 				::	 		InputValW(size(InputVal)+1)
    complex(dpc),allocatable	::			OutputValW(:)
    complex(dpc),allocatable	::			dMultFactor(:)
    integer(i8b) 				::			InputLen, OutputLen, EvenOrOdd, iNumPlanTransform, iNumPlanTransform_inv, iDimW

    !Filtering and windowing parameters 
    real(dp) 	::			dTaperMaxFreqWindow(size(InputVal)/2+1), dLeftBand,dRightBand, dDomega,dFFTFactor
    real(dp)	::			arBuffer1(size(InputVal)+1),arBuffer2(size(InputVal)+1)
    integer		::			i,iErr

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
    InputLen  = size(InputVal)
    EvenOrOdd = mod(InputLen,2_i8b)
    iDimW = InputLen /2+1

    InputValW = 0.0
    InputValW(1:InputLen/2) = InputVal(1:InputLen-EvenOrOdd:2) + im * InputVal(2:InputLen-EvenOrOdd:2)
    if(EvenOrOdd==1) then
       InputValW(InputLen/2+EvenOrOdd) = InputVal(InputLen)
    endif

    ! Do a FFT of InputLen-point (Initial Variable Length) with complex numbers 
    call dfftw_plan_dft_r2c_1d(iNumPlanTransform, InputLen, InputValW, InputValW, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iNumPlanTransform, InputValW, InputValW)
    call dfftw_destroy_plan(iNumPlanTransform)

    ! Tapering of the max frequency because of artifact created due to higher frequencies
    dLeftBand = 0.0
    dRightBand =0.1_dp  
    dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(InputLen*cSpace%dDt),dLeftBand,dRightBand) 

    InputValW(1:iDimW) = InputValW(1:iDimW) * dTaperMaxFreqWindow

    ! Inverse FFT of OutputLen-point	
    call dfftw_plan_dft_c2r_1d(iNumPlanTransform_inv, 2*InputLen, InputValW, InputValW, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(iNumPlanTransform_inv, InputValW, InputValW)
    call dfftw_destroy_plan(iNumPlanTransform_inv)

    arBuffer1 = real(InputValW(1:InputLen),dp)
    arBuffer2 = dimag(InputValW(1:InputLen))

    InputValW(1:InputLen) = arBuffer1**OrderP + im * arBuffer2**OrderP

    !Allocate to upsample or decimate . OutputValW has to have more or equal points to the initial signal 
    ALLOCATE(OutputValW(max(InputLen+1,OutputLen+1)))

    OutputValW = 0.0D0
    ! Do a FFT of InputLen-point (Initial Variable Length) with complex numbers 
    call dfftw_plan_dft_r2c_1d(iNumPlanTransform, 2*OutputLen, InputValW, OutputValW, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(iNumPlanTransform, InputValW, OutputValW)
    call dfftw_destroy_plan(iNumPlanTransform)    

    ALLOCATE(dMultFactor(iDimW))
    dMultFactor = 1.0D0
    dFFTFactor = 1.0_dp/(2.0_dp*real(cSpace%cGrid%iD0TL,dp)**(2+OrderP-1))
    if (OrderFFT /= 0) then
       dDOmega = two_pi /real(InputLen * cSpace%dDt,dp) 
       dMultFactor(1:iDimW) = ( (/ (i,i=0,iDimW-1) /) * dDOmega * im)**OrderFFT * dFFTFactor

    endif

    OutputValW(1:iDimW) = OutputValW(1:iDimW) * dTaperMaxFreqWindow  * dMultFactor(1:iDimW)

    ! Inverse FFT of OutputLen-point	
    call dfftw_plan_dft_c2r_1d(iNumPlanTransform_inv, OutputLen, OutputValW, OutputValW, FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(iNumPlanTransform_inv, OutputValW, OutputValW)
    call dfftw_destroy_plan(iNumPlanTransform_inv)

    OutputVal(1:OutputLen:2) = real(OutputValW(1:OutputLen/2+EvenOrOdd),dp)
    OutputVal(2:OutputLen:2) = dimag(OutputValW(1:OutputLen/2))
    OutputVal = OutputVal * Magnitude;

    DEALLOCATE(OutputValW)
    DEALLOCATE(dMultFactor)

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
    ! 	based on the values inserted in input. The slices are generated based on the
    !	number of slices in each dimension (x,y or z). Then the domain of the bubble
    !	cloud is divided in slices of a number equal to the aforementioned number.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace             io  type(space)  Space from which the data array needs to be stored.
    !   dBubbleLocationN   r   dp   		a matrix of the bubbles' location , 1st dim Bubble No., 2nd dim dimension (x,y,z)
    !   Domain_Range	   r   dp			Range of the domain that includes the microbubble cloud
    !
    ! *****************************************************************************
    real(dp), intent(inout)					::		dBubbleLocationN(:,:)
    real(dp),intent(in) 					::		Domain_Range(2,3)
    type(Space),  intent(in)				::		cSpace


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
    character(len = 1024)	::  	filename, actemp
    real(dp)				::      DimStart(3) , CountPos(3), dLambdaMM
    integer(i8b)			::		i , j(3) , slicedim, ExtraSliceNum
    character(LEN=1)		:: 		xyzslicedim(cModelParams%Numslices)		!Which dimension: x, y, z or t
    real(dp) 				:: 		xyzslicepos(cModelParams%Numslices)		!Position in [mm]
    integer(i8b)			:: 		xyzslicebeam(cModelParams%Numslices)	!Which beam the slice is located; -1 if it may be all
    integer(i8b) 			:: 		xyzsliceindex(cModelParams%Numslices)	!Index within the beam
    character(len =1)		::		allslicexyzdim(BubbleParams%ClusterSlicesN)	

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
    !02	  Spectral difference by multiplying (iw)^n, ( n = 2 for second derivative)
    !03   IFFT and return to time domain signal 
    !04	  Multiply with the constant factor A (magnitude)  	
    !
    ! =============================================================================

    ExtraSliceNum  =  BubbleParams%ClusterSlicesN
    allslicexyzdim =  BubbleParams%ClusterSliceDim
    dLambdaMM = (cMediumParams.c0*1.0D3) /cModelParams.freq0;           ! Normalization factor 

    if (ExtraSliceNum == 0) then 
       write(acTemp,"('0 Slices were chosen in the input file')") 
       call PrintToLog(acTemp,2)
       return
    endif

    xyzslicedim = (cModelParams%xyzslicedim)
    xyzslicepos = (cModelParams%xyzslicepos)
    xyzslicebeam = (cModelParams%xyzslicebeam)
    xyzsliceindex = (cModelParams%xyzsliceindex)

    ! Deallocate to change dimension
    DEALLOCATE(cModelParams%xyzslicedim)
    DEALLOCATE(cModelParams%xyzslicepos)
    DEALLOCATE(cModelParams%xyzslicebeam)
    DEALLOCATE(cModelParams%xyzsliceindex)

    ! Change dimensions and reallocate
    cModelParams%Numslices = cModelParams%Numslices + ExtraSliceNum

    ALLOCATE(cModelParams%xyzslicedim(cModelParams%Numslices))
    ALLOCATE(cModelParams%xyzslicepos(cModelParams%Numslices))
    ALLOCATE(cModelParams%xyzslicebeam(cModelParams%Numslices))
    ALLOCATE(cModelParams%xyzsliceindex(cModelParams%Numslices))

    cModelParams%xyzslicedim(1:cModelParams%Numslices-ExtraSliceNum) = xyzslicedim 
    cModelParams%xyzslicepos(1:cModelParams%Numslices-ExtraSliceNum) = xyzslicepos 
    cModelParams%xyzslicebeam(1:cModelParams%Numslices-ExtraSliceNum) = xyzslicebeam
    cModelParams%xyzsliceindex(1:cModelParams%Numslices-ExtraSliceNum) = xyzsliceindex 
    ! Count how many slices are for each dimension to divide the domain 
    cModelParams%xyzslicedim(cModelParams%Numslices-ExtraSliceNum+1:cModelParams%Numslices) = allslicexyzdim

    CountPos = (/ 	count(allslicexyzdim == 'x') , &
         count(allslicexyzdim == 'y') , &
         count(allslicexyzdim == 'z') /)

    j = 1
    DimStart = 	(/cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ/)*cSpace%dDx  
    ! Separate the domain of the bubble cloud in a number of slices equal to the number of slices given by the user in the input file
    do i = ExtraSliceNum,1,-1

       slicedim = sum(abs((/ 'x', 'y', 'z'/) ==  cModelParams%xyzslicedim(cModelParams%Numslices - i+1) ) * (/ 1,2,3 /))
       ! If the counter is 1 , save the average position of the bubble cluster in this dimension
       if (ExtraSliceNum == 1) then
          cModelParams%xyzslicepos(cModelParams%Numslices-i+1) = sum((sum(dBubbleLocationN,1)/size(dBubbleLocationN,1) &
               + DimStart) *( abs((/ 'x', 'y', 'z'/) ==  cModelParams%xyzslicedim(cModelParams%Numslices - i+1) )))	 
       else 
          cModelParams%xyzslicepos(cModelParams%Numslices-i+1) = Domain_Range(1,slicedim)/dLambdaMM &
               + (Domain_Range(2,slicedim) - Domain_Range(1,slicedim))/(CountPos(slicedim)-1) * (j(slicedim) -1)/dLambdaMM
       endif

       j(slicedim) = j(slicedim) +1

    enddo

    if((cModelParams%Slicesavespecifier==iAI_FIRSTANDLAST) .or.(cModelParams%Slicesavespecifier==iAI_ALL)) then
       ! Initialize slices, because of the addition of the extra slices
       call InitSaveSlices(cSpace)

       do i =  ExtraSliceNum,1,-1     
          filename = trim(trim(sOutputDir) // trim(sOutputRoot) // int2str(0_i8b))//'_'//cModelParams%xyzslicedim(cModelParams%Numslices-i+1)//&
               int2str(cModelParams%Numslices-i+1)//int2str(0)
          ! Export slice

          if (cModelParams%xyzslicedim(cModelParams%Numslices-i+1)=='x') then
             call ExportSlice(trim(filename), "p0", cSpace, &
                  (/ 0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices-i+1), 0_i8b,0_i8b /), &
                  (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                  cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
          elseif (cModelParams%xyzslicedim(cModelParams%Numslices-i+1)=='y') then
             call ExportSlice(trim(filename), "p0", cSpace, &
                  (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices-i+1), 0_i8b /), &
                  (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                  cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
          elseif (cModelParams%xyzslicedim(cModelParams%Numslices-i+1)=='z') then
             call ExportSlice(trim(filename), "p0", cSpace, &
                  (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(cModelParams%Numslices-i+1) /), &
                  (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                  cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
          end if
       enddo
    endif

  END SUBROUTINE SaveExtraBubbleSlices


END MODULE ParnacBubbleContrast
