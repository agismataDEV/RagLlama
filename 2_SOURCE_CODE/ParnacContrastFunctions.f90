MODULE ParnacContrastFunctions
  
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
  USE ContrastGenFilt !! L.D. 03-11-2009
  
  
  USE ParnacDataRedist, ONLY : &
  ReorderDistr0ToDistr1,ReorderDistr1ToDistr0, &
  ReorderDistr1ToDistr2,ReorderDistr2ToDistr1, &
  Distr2ObtainXYZBlock, Distr2PutXYZBlock,  &
  Distr2ObtainYMirroredXYZBlock,Distr2FillYMirroredXYZBlock, &
  Distr2PutYMirroredXYZBlock
  USE ParnacDataSpaceInit, ONLY : &
  InitSpace, InitGrid, SpaceInfo, DestructSpace, &
  GridDistr0CreateEmpty, GridDistr0Deallocate, &
  GridDistr1CreateEmpty,GridDistr2CreateEmpty, iBeamOffsetX, iBeamOffsetY, &
  InitSaveSlices
  USE ParnacDerivative, ONLY : &
  iFDMinOrder, DerivLookup, &
  DerivLookupInit, DerivLookupDestroy, &
  DerivativeReal,DerivativeComplex
  USE ParnacDataMap, ONLY : &
  MapVx, MapVy, MapVxInv, MapVyInv, &
  MapVx_Small, MapVy_Small, MapVxInv_Small, MapVyInv_Small, &
  MapVt, MapVtInv
  USE ParnacTransformFilter, ONLY : &
  TransformT, TransformTInv, &
  TransformT_sml, TransformTInv_sml, &
  TransformXYZ, TransformXYZInv, &
  TransformXYZ_sml, TransformXYZInv_sml, &
  FFilterSpaceSpatial3D, dTaperingWindow, INTERP1DFREQ
  USE ParnacOutput, ONLY : &
  ExportSlice
  USE SpecialFun,ONLY: &
  INTERP1D,LINSPACE
  USE FFTW_cons

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
  !   NonlinContrastOperator   sub   Evaluate the nonlinearity contrast operator
  !                                  without anti-aliasing procedure
  !   NonlinContrastOperator_Ali sub  Evaluate the nonlinear contrast operator
  !                                  with anti-aliasing procedure
  !   InhomContrastOperator    sub   Evaluate the inhomogeneity contrast source
  !                                  without anti-aliasing procedure
  !   InhomContrastOperator_Ali sub  Evaluate the inhomogeneity contrast source
  !                                  with anti-aliasing procedure
  !   InitContrastSpace        sub   Fills a space with a specific choice for  
  !                                  the inhomogeneity contrast function
  !   dpointScatterercontrast  fun   Return a point in space of a point scatterer      
  !   dGaussianblobcontrast    fun   Return a point for the Gaussian blob
  !   dSpherecontrast          fun   Return a point for the Sphere
  !   dLunebergLensContrast    fun   Return a point for the Luneberg lens
  !   dFattylayerContrast      fun   Return a point in space for the phase screen
  !   dVeinContrast            fun   Return a point in space for a cylinder
  !   SpreadFrequencies        sub   Spread frequencies of half a frequency axis 
  !                                  over the full frequency axis
  !   CollectFrequencies       sub   Collect frequencies of a half frequency axis
  !                                  from a full frequency axis
  !   StoreContrast            sub   Store the contrast space to a temporary file
  !   LoadContrast              sub   Get the contrast space from a temporary file
  !   ConvertSmalltoAliasingBlock  sub  Rearrange data in an XYZ slice as though
  !                                  it includes an anti-aliasing region in all
  !                                  three dimensions
  !   ConvertAliasingBlocktoSmall  sub  Rearrange data in an XYZ slice without
  !                                  the anti-aliasing region in all three 
  !                                  dimensions
  !   MoveBlock                sub   copy blocks within an XYZ slice - used by
  !                                  ConvertSmalltoAliasingBlock and
  !                                  ConvertAliasingBlocktoSmall
  !   dTaperingWindow          fun   Generate a window with cosine tapering
  !   TaperXYZslice            sub   Apply tapering to a space at the beginning 
  !                                  and end of the X/Y/Z axes
  !
  ! =============================================================================
  !
  CONTAINS
  
  SUBROUTINE NonlinearKappaOperator(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    REAL(8), dimension(:,:,:,:,:), ALLOCATABLE  :: PressureMatrix
    REAL(8), dimension(:,:,:,:), ALLOCATABLE    :: PressureMatrixGlobal, PressureMatrixGlobalX, PressureMatrixGlobalY, PressureMatrixGlobalZ,PressureMatrixGlobalDerivativeX,PressureMatrixGlobalDerivativeY,PressureMatrixGlobalDerivativeZ
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1x(cSpace.iDimX),arBuffer2x(cSpace.iDimX),arBuffer1y(cSpace.iDimY),arBuffer2y(cSpace.iDimY),arBuffer1z(cSpace.iDimZ),arBuffer2z(cSpace.iDimZ)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDX, KAPPAcontrastDY, KAPPAcontrastDZ 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          !print *, "Barrier7"
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex);
          
          call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          !print *, "Barrier8"
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=arBuffer2x(1:iLen);
          !print *, "Barrier9"
          
       end do
    end do
    
    
    ! derivative with respect to Y
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex);
          
          call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=arBuffer2y(1:iLen);
          
          
       end do
    end do
    
    ! derivative with respect to Z
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ);
          
          call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=arBuffer2z(1:iLen);
          
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDY
    !CLOSE(3, STATUS = 'keep') 
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDZ
    !CLOSE(3, STATUS = 'keep')  
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    !print *, iDimX
    !print *, (pcgrid%aid0loc(pcGrid%iD0LocN,1))
    
    
    
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1x=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1x(pcgrid%aid0loc(iindex+1,1)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2x(pcgrid%aid0loc(iindex+1,1)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    DEALLOCATE (KAPPAcontrastDX)
    DEALLOCATE (KAPPAcontrastDY)
    DEALLOCATE (KAPPAcontrastDZ)  
    
  END SUBROUTINE NonlinearKappaOperator
  
  SUBROUTINE NonlinearKappaOperatorDX(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1x(cSpace.iDimX),arBuffer2x(cSpace.iDimX)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDX
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    call TransformTInv(pcGrid)
    
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    
    call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    call ReorderDistr1ToDistr0(pcGrid)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          !print *, "Barrier7"
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex);
          
          call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          !print *, "Barrier8"
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=arBuffer2x(1:iLen);
          !print *, "Barrier9"
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1x=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1x(pcgrid%aid0loc(iindex+1,1)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2x(pcgrid%aid0loc(iindex+1,1)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    DEALLOCATE (KAPPAcontrastDX)
    
  END SUBROUTINE NonlinearKappaOperatorDX
  
  SUBROUTINE NonlinearKappaOperatorDY(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dy(dy(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1y(cSpace.iDimY),arBuffer2y(cSpace.iDimY)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDY 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)   
    ! derivative with respect to Y
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex);
          
          call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=arBuffer2y(1:iLen);
          
          
       end do
    end do
    
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDY
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1y=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1y(pcgrid%aid0loc(iindex+1,2)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2y(pcgrid%aid0loc(iindex+1,2)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    
    DEALLOCATE (KAPPAcontrastDY) 
    
  END SUBROUTINE NonlinearKappaOperatorDY
  
  SUBROUTINE NonlinearKappaOperatorDZ(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dz(dz(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1z(cSpace.iDimZ),arBuffer2z(cSpace.iDimZ)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDZ 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    ! derivative with respect to Z
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ);
          
          call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=arBuffer2z(1:iLen);
          
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDZ
    !CLOSE(3, STATUS = 'keep')  
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    !print *, iDimX
    !print *, (pcgrid%aid0loc(pcGrid%iD0LocN,1))
    
    
    
    ! derivative with respect to Z
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          do tindex=1,iDimT
             
             arBuffer1z=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   arBuffer1z(pcgrid%aid0loc(iindex+1,3)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2z(pcgrid%aid0loc(iindex+1,3)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    DEALLOCATE (KAPPAcontrastDZ)  
    
  END SUBROUTINE NonlinearKappaOperatorDZ
  
  !-------------DIRECTION
  SUBROUTINE NonlinearKappaOperatorDirectionDX(cSpace,IncSpace)
    
    ! =============================================================================
    
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1x(cSpace.iDimX),arBuffer2x(cSpace.iDimX)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDX 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperatorDirection",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid !this contains the pressure field
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          !print *, "Barrier7"
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex);
          
          call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          !print *, "Barrier8"
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=arBuffer2x(1:iLen);
          !print *, "Barrier9"
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    
    
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1x=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1x(pcgrid%aid0loc(iindex+1,1)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2x(pcgrid%aid0loc(iindex+1,1)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    
    iStart=0
    ! Now we have to multiply the two to get p(dxkx)dx(G(r(p)))
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    DEALLOCATE (KAPPAcontrastDX)
    
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDX
  !-------------DIRECTION
  
  !-------------DIRECTION
  SUBROUTINE NonlinearKappaOperatorDirectionDY(cSpace,IncSpace)
    
    ! =============================================================================
    
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1y(cSpace.iDimY),arBuffer2y(cSpace.iDimY)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDY 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperatorDirection",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid !this contains the pressure field
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)   
    ! derivative with respect to Y
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex);
          
          call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=arBuffer2y(1:iLen);
          
          
       end do
    end do
    
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)	
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1y=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1y(pcgrid%aid0loc(iindex+1,2)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2y(pcgrid%aid0loc(iindex+1,2)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    
    iStart=0
    ! Now we have to multiply the two to get p(dxkx)dx(G(r(p)))
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    DEALLOCATE (KAPPAcontrastDY)
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDY
  !-------------DIRECTION
  
  !-------------DIRECTION
  SUBROUTINE NonlinearKappaOperatorDirectionDZ(cSpace,IncSpace)
    
    ! =============================================================================
    
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1z(cSpace.iDimZ),arBuffer2z(cSpace.iDimZ)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDZ 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperatorDirection",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid !this contains the pressure field
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    ! derivative with respect to Z
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ);
          
          call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=arBuffer2z(1:iLen);
          
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)	
    ! derivative with respect to Z
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          do tindex=1,iDimT
             
             arBuffer1z=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   arBuffer1z(pcgrid%aid0loc(iindex+1,3)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2z(pcgrid%aid0loc(iindex+1,3)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup)
    
    
    iStart=0
    ! Now we have to multiply the two to get p(dxkx)dx(G(r(p)))
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    DEALLOCATE (KAPPAcontrastDZ)
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDZ
  !-------------DIRECTION
  
  SUBROUTINE NonlinearKappaOperatorLinDX(cSpace, IncSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = nabla(nabla(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)  ::	cSpace
    type(Space), target, intent(in)     ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1x(cSpace.iDimX),arBuffer2x(cSpace.iDimX)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDX
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Calculating p^2
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    call ReorderDistr0ToDistr1(pcGridIn)
    call TransformT_sml(pcGridIn)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          !print *, "Barrier7"
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex);
          
          call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          !print *, "Barrier8"
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=arBuffer2x(1:iLen);
          !print *, "Barrier9"
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDX
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1x=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1x(pcgrid%aid0loc(iindex+1,1)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1x, arBuffer2x, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,2)+1)==yindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2x(pcgrid%aid0loc(iindex+1,1)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    DEALLOCATE (KAPPAcontrastDX)
    
  END SUBROUTINE NonlinearKappaOperatorLinDX
  
  
  
  SUBROUTINE NonlinearKappaOperatorLinDY(cSpace, IncSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = nabla(nabla(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)  ::	cSpace
    type(Space), target, intent(in)     ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1y(cSpace.iDimY),arBuffer2y(cSpace.iDimY)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDY
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Calculating p^2
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    call ReorderDistr0ToDistr1(pcGridIn)
    call TransformT_sml(pcGridIn)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)   
    ! derivative with respect to Y
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex);
          
          call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=arBuffer2y(1:iLen);
          
          
       end do
    end do
    
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDY
    !CLOSE(3, STATUS = 'keep') 
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          do tindex=1,iDimT
             
             arBuffer1y=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   arBuffer1y(pcgrid%aid0loc(iindex+1,2)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1y, arBuffer2y, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,3)+1)==zindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2y(pcgrid%aid0loc(iindex+1,2)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    
    DEALLOCATE (KAPPAcontrastDY) 
    
  END SUBROUTINE NonlinearKappaOperatorLinDY
  
  
  
  SUBROUTINE NonlinearKappaOperatorLinDZ(cSpace, IncSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = nabla(nabla(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)  ::	cSpace
    type(Space), target, intent(in)     ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1z(cSpace.iDimZ),arBuffer2z(cSpace.iDimZ)
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered, KAPPAContrast, KAPPAcontrastDZ 
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = cSpace%iDimT/2 + 1
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Calculating p^2
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    call ReorderDistr0ToDistr1(pcGridIn)
    call TransformT_sml(pcGridIn)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    
    ! now we have to calculate dxkappa dykappa dzkappa and multiply each term with p^2 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !if (iProcID==1) then
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrast
    !CLOSE(3, STATUS = 'keep')  
    !end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    ! derivative with respect to Z
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ);
          
          call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=arBuffer2z(1:iLen);
          
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    call DerivLookupDestroy(cDerivLookup);
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) KAPPAcontrastDZ
    !CLOSE(3, STATUS = 'keep')  
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to thi term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    ! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    
    !print *, iDimX
    !print *, (pcgrid%aid0loc(pcGrid%iD0LocN,1))
    
    
    
    ! derivative with respect to X
    iLen	= iDimZ;
    do xindex=1,iDimX
       do yindex=1,iDimY
          do tindex=1,iDimT
             
             arBuffer1z=0
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   arBuffer1z(pcgrid%aid0loc(iindex+1,3)+1) = pcGrid%parD0(iStart+tindex)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             call DerivativeReal(arBuffer1z, arBuffer2z, iLen, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
             
             
             iStart=0
             do iIndex=0,pcGrid%iD0LocN-1
                
                if (((pcgrid%aid0loc(iindex+1,1)+1)==xindex).AND.((pcgrid%aid0loc(iindex+1,2)+1)==yindex)) then
                   pcGrid%parD0(iStart+tindex)=arBuffer2z(pcgrid%aid0loc(iindex+1,3)+1)
                end if
                
                iStart=iStart+pcGrid%iD0IS
                
             end do
             
             
          end do
       end do
    end do
    call DerivLookupDestroy(cDerivLookup);
    ! Now pcGrid contains dx(dxk p^2) 
    
    DEALLOCATE (KAPPAcontrastDZ)  
    
  END SUBROUTINE NonlinearKappaOperatorLinDZ
  
  SUBROUTINE NonlinContrastOperator(cSpace)
    
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
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^NL(p) = kappa0 * beta * d^2 / dt^2 p^2 for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the temporal dimension. The 2nd order temporal 
    !   derivative is implemented with a Finite Difference approximation; we use
    !   d^2/dt^2 p^2 = d/dt (2 p d/dt p).
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iIndex        i8b   loop counter over the positions in the space
    !   iLen          i8b   length of a time trace
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1     dp    Temporary buffer containing a time trace
    !   arBuffer2     dp    Temporary buffer containing a time trace
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    integer(i8b)::				iIndex, iLen, iStart,i
    type(DerivLookup)::			cDerivLookup
    real(dp)::				arBuffer1(cSpace.iDimT),arBuffer2(cSpace.iDimT),dCorrection;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy
    !   
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   Initialize the cDerivLookup structure with the weights of the finite 
    !     difference stencil
    !   For each position in the space, do
    !       compute the finite difference d/dt p
    !       compute 2 * p * d/dt p
    !       compute d/dt ( 2 * p * d/dt p)
    !   Apply the factor kappa_0 * beta in the contrast source
    !   Deinitialize the cDerivLookup structure
    !
    ! =============================================================================
    
    call PrintToLog("NonlinOperator",4)
    
    ! Initialize
    call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    iStart	= 0;
    iLen	= cSpace.iDimT;
    
    ! Now we process t-block
    do iIndex = 0, cSpace.cGrid.iD0LocN-1
       
       arBuffer1=cSpace.cGrid.parD0(1+iStart: iStart+iLen);
       
       call DerivativeReal(arBuffer1, &
       arBuffer2, &
       iLen, &
       cSpace.dDt, &
       cDerivLookup.arWeights, &
       cDerivLookup.aiPoints)
       
       arBuffer1		= 2 * arBuffer2(1:iLen) * arBuffer1(1:iLen);
       
       call DerivativeReal(arBuffer1, arBuffer2, &
       iLen, &
       cSpace.dDt, &
       cDerivLookup.arWeights, &
       cDerivLookup.aiPoints)
       
       cSpace.cGrid.parD0(1+iStart: iStart+iLen)=arBuffer2(1:iLen);
       
       iStart		= iStart + cSpace.cGrid.iD0TL;
       
    end do
    
    call PrintToLog("Multiply with the correction factor", 2);
    dCorrection	= cMediumParams%Kappa0 * cMediumParams%Beta
    do i = 1, cSpace%cGrid%iD0LocSize
       cSpace%cGrid%parD0(i)	= cSpace%cGrid%parD0(i) * dCorrection;
    end do
    
    ! DeInitialize
    call DerivLookupDestroy(cDerivLookup);
    
  END SUBROUTINE NonlinContrastOperator
  
  !SUBROUTINE NonlinContrastOperator_AliINHOMOGENEOUSCONTRASTINC(cSpace)
  !
  !! =============================================================================
  !!
  !!   Programmer: Koos Huijssen
  !!
  !!   Language: Fortran 90
  !!
  !!   Version Date    Comment
  !!   ------- -----   -------
  !!   1.0     090505  Original code (KH)
  !!
  !! *****************************************************************************
  !!
  !!   DESCRIPTION
  !!
  !!   The subroutine NonlinContrastOperator_Ali computes the nonlinearity 
  !!   contrast source S^NL(p) = kappa0 * beta * d^2 / dt^2 p^2 for the given 
  !!   space. The data in cSpace should be in Distribution 0. This implementation 
  !!   employs an anti-aliasing approach for the temporal dimension. The temporal
  !!   derivative is implemented in the Fourier domain with a spectral derivative.
  !!
  !! *****************************************************************************
  !!
  !!   INPUT/OUTPUT PARAMETERS
  !!
  !!   cSpace   io   type(space)  space for which the nonlinearity contrast source
  !!                              is determined
  !!
  !	type(Space), target, intent(inout)::		cSpace
  !	
  !! *****************************************************************************
  !! 
  !!   LOCAL PARAMETERS      
  !!
  !!   pcGrid    type(Grid)  pointer to the grid in cSpace
  !!   iIndex        i8b   loop counter over the positions in the space
  !!   iStart        i8b   start of each time trace
  !!   i             i8b   loop counter
  !!   iDimW         i8b   number of points in the angular frequency domain
  !!   arBuffer1     dp    Temporary buffer containing a time/frequency trace
  !!   arBuffer2     dp    Temporary buffer containing a time/frequency trace
  !!   dDOmega       dp    Step size for the frequency axis
  !!   dFFTFactor    dp    Correction factor from the FFT's
  !!   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
  !!                            end
  !!   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
  !!                            frequency components
  !!   dMultFactor   dp    Multiplication factor accounting for the derivative
  !!                       term -omega**2 in the frequency domain
  !!   dLeftBand     dp    Size of the tapering window at the beginning
  !!   dRightBand    dp    Size of the tapering window at the end
  !!   dCorrection   dp    Correction factor of the medium parameters 
  !!
  !	type(Grid), pointer :: pcGrid
  !	integer(i8b)::		   iIndex, iStart, i, iDimW;
  !	real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT)
  !	real(dp) ::			   dDOmega, dFFTFactor
  !	real(dp), allocatable :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:)
  !	real(dp) ::		dLeftBand,dRightBand,dCorrection
  !	
  !! *****************************************************************************
  !!
  !!   I/O
  !!
  !!   log file entries
  !!   
  !! *****************************************************************************
  !!
  !!   SUBROUTINES/FUNCTIONS CALLED
  !!
  !!   PrintToLog
  !!   dTaperingWindow
  !!   ReorderDistr0ToDistr1
  !!   TransformT_sml
  !!   TransformTInv
  !!   TransformT
  !!   TransformTInv_sml
  !!   ReorderDistr1ToDistr0
  !!   
  !! *****************************************************************************
  !!
  !!   PSEUDOCODE
  !!
  !!   If required, apply the tapering window to the field that removes the sharp
  !!     edge at start and end of the time trace
  !!   Reorder to distribution 1
  !!   Transform the data in the time dimension - use 'small' transform (iDimT)
  !!   If required, apply the frequency tapering window that removes the highest
  !!     frequencies - these may be unrealistically large in the iterative scheme
  !!   Inverse transform with a 'large' transform (2*iDimT) that includes the
  !!     anti-aliasing frequency range
  !!   Apply the square p**2
  !!   Transform to the frequency domain with a 'large' transform (2*iDimT)
  !!   If required, again apply the frequency tapering window
  !!   Apply the FFT factors and the derivative term -omega**2
  !!   Inverse transform with a 'small' transform (iDimT)
  !!   Reorder to distribution 0
  !!   Apply the factor kappa_0 * beta in the contrast source
  !!
  !! =============================================================================
  !	
  !	call PrintToLog("NonlinOperator_Ali",4)
  !
  !	pcGrid=>cSpace%cGrid
  !
  !	iDimW = cSpace%iDimT/2 + 1
  !
  !	allocate(dTaperSupportWindow(1:cSpace%iDimT),&
  !			dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW))
  !
  !	if (cModelParams.UseSupportTapering .EQV. .true.) then
  !		!use tapering of the field at start and end to prevent wraparound leakage
  !
  !		!build tapering window, the first two periods and the last two periods are tapered
  !		!in the contrast source
  !		dLeftBand = 2.0_dp
  !		dRightBand = 2.0_dp
  !		dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
  !
  !		!multiply the field with it
  !		iStart = 0
  !		do iIndex=0,pcGrid%iD0LocN-1
  !			pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
  !					pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
  !
  !			iStart=iStart+pcGrid%iD0IS
  !		end do
  !
  !	end if
  !
  !	!First, transform to W-domain (use small transform, no wraparound regions)
  !	call ReorderDistr0ToDistr1(pcGrid)
  !	call TransformT_sml(pcGrid)
  !
  !	if (cModelParams.UseFreqTapering .EQV. .true.) then
  !	!Tapering of the highest frequency part; otherwise the chopoff noise around the 
  !	! Nyquist frequency is being distributed over the entire frequency axis by the p**2
  !	! and is blown up by the double derivative...
  !
  !		!build tapering window, the highest .1*f0 are tapered in the contrast source
  !		dLeftBand = 0
  !		dRightBand = 0.1_dp
  !		dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
  !
  !		!multiply the field with it
  !		iStart = 0
  !		do iIndex=0,pcGrid%iD1LocN-1
  !			pcGrid%pacD1(iStart+1:iStart+iDimW) = &
  !					pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
  !
  !			iStart=iStart+pcGrid%iD1IS
  !		end do
  !
  !	end if
  !
  !	!Now, transform back to T-domain as though it was a grid with wraparound regions
  !	!however, the wraparound regions are now used as anti-aliasing regions...
  !	call TransformTInv(pcGrid)
  !
  !	!Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
  !	!Only real multiplication..
  !	iStart=0
  !	do iIndex = 0, cSpace.cGrid.iD1LocN-1
  !
  !		arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
  !		arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
  !
  !		pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
  !			arBuffer1**2 + im * arBuffer2**2;		
  !
  !		iStart		= iStart + cSpace.cGrid.iD1IS;
  !
  !	end do
  !
  !!call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr6",&
  !!						"Contrast", &
  !!						cSpace, &
  !!						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
  !!						(/ 2*cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
  !!						.true.);
  !
  !	!Transform to W-domain and multiply with -w**2 and with 
  !	! FFT normalization factor (1/N for each FORWARD transform, 
  !	! and an extra square because of the square factor above;
  !	! this gives (1/N)**2*(1/(2*N)) )
  !	call TransformT(pcGrid)
  !
  !	dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
  !	dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
  !	dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
  !    
  !	if (cModelParams.UseFreqTapering .EQV. .true.) then
  !	!Here, the max freq tapering is applied as well
  !		dMultFactor = dMultFactor * dTaperMaxFreqWindow
  !	end if
  !
  !	iStart=0
  !	do iIndex = 0, cSpace.cGrid.iD1LocN-1
  !
  !		pcGrid%pacD1(iStart+1:iStart+iDimW) = &
  !			pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
  !
  !		iStart		= iStart + cSpace.cGrid.iD1TL;
  !
  !	end do
  !
  !!call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr8",&
  !!						"Contrast", &
  !!						cSpace, &
  !!						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
  !!						(/ cSpace%iDimT+1, 1_i8b, 1_i8b, cSpace%iDimZ /), &
  !!						.false.);
  !
  !	!Finally, transform back to T-domain with small transform, 
  !	! without the anti-aliasing regions, and reorder to dist. 0
  !	call TransformTInv_sml(pcGrid)
  !	call ReorderDistr1ToDistr0(pcGrid)
  !
  !!call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr10",&
  !!						"Contrast", &
  !!						cSpace, &
  !!						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
  !!						(/ cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
  !!						.true.);
  !
  !	call PrintToLog("Multiply with the correction factor", 2);
  !	dCorrection	= cMediumParams%Kappa0 * cMediumParams%Beta
  !	do i = 1, pcGrid%iD0LocSize
  !		pcGrid%parD0(i)	= pcGrid%parD0(i) * dCorrection;
  !	end do
  !
  !	deallocate(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor)
  !
  !END SUBROUTINE NonlinContrastOperator_AliINHOMOGENEOUSCONTRASTINC
  
  
!!!L.D.
  SUBROUTINE NonlinContrastOperator_Ali(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     061009  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator_Ali computes the nonlinearity 
    !   contrast source S^NL(p) = kappa0 * beta * d^2 / dt^2 p^2 for the given 
    !   space. The data in cSpace should be in Distribution 0. This implementation 
    !   employs an anti-aliasing approach for the temporal dimension. The temporal
    !   derivative is implemented in the Fourier domain with a spectral derivative.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    type(Grid), pointer :: pcGrid
    integer(i8b)::		   iIndex, iStart, i, iDimW, iDimZ, iDimX, iDimY, INDEXFORZ;
    real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT)
    real(dp) ::			   dDOmega, dFFTFactor
    real(dp) 					:: dLeftBand,dRightBand,dCorrection
    real(dp) 					:: alpha1,alpha2,alpha3,alpha4,alpha5,b_att0,b_att1,b_att2,b_att3,b_att4,b_att5
    complex(dpc), allocatable 	:: acBuffer(:), AttenContrast1(:),AttenContrast2(:),AttenContrast3(:),AttenContrast4(:),AttenContrast5(:)
    real(dp), allocatable 		:: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:), ContrastFiltered(:,:,:) !! L.D. 03-11-2009
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv
    !   TransformT
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   If required, apply the tapering window to the field that removes the sharp
    !     edge at start and end of the time trace
    !   Reorder to distribution 1
    !   Transform the data in the time dimension - use 'small' transform (iDimT)
    !   If required, apply the frequency tapering window that removes the highest
    !     frequencies - these may be unrealistically large in the iterative scheme
    !   Inverse transform with a 'large' transform (2*iDimT) that includes the
    !     anti-aliasing frequency range
    !   Apply the square p**2
    !   Transform to the frequency domain with a 'large' transform (2*iDimT)
    !   If required, again apply the frequency tapering window
    !   Apply the FFT factors and the derivative term -omega**2
    !   Inverse transform with a 'small' transform (iDimT)
    !   Reorder to distribution 0
    !   Apply the factor kappa_0 * beta in the contrast source
    !
    ! =============================================================================
    
    call PrintToLog("NonlinOperator_Ali",4)
    
    pcGrid=>cSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    iDimZ = cSpace%iDimZ
    iDimX = cSpace%iDimX
    iDimY = cSpace%iDimY
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW), acBuffer(1:pcgrid%iD1LocSize))
	iMemAllocated=iMemAllocated + ( PRODUCT(SHAPE(dTaperSupportWindow)) + PRODUCT(SHAPE(dTaperMaxFreqWindow)) + PRODUCT(SHAPE(dMultFactor)))* dpS + PRODUCT(SHAPE(acBuffer)) * dpcS
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow	
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !!L.D. 
    ! acBuffer now contains the tapered pressure in w-domain
    acBuffer=pcGrid%pacD1
    
    !!L.D.
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1**2 + im * arBuffer2**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex
       
    end do
    
    !Transform to W-domain and multiply with -w**2 and with 
    ! FFT normalization factor (1/N for each FORWARD transform, 
    ! and an extra square because of the square factor above;
    ! this gives (1/N)**2*(1/(2*N)) )
    call TransformT(pcGrid)
    
    dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Here, the max freq tapering is applied as well
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
    end if
    
    iStart=0
	
    dCorrection	= cMediumParams%Kappa0 * cMediumParams%Beta
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor* dCorrection ;		
       
       acBuffer(iStart+1:iStart+iDimW)=acBuffer(iStart+1:iStart+iDimW)* 2 * dMultFactor*real(pcGrid%iD0TL,dp)**2;
       !pcGrid%pacD1(iStart+1:iStart+iDimW)= pcGrid%pacD1(iStart+1:iStart+iDimW) * dCorrection    ! Delete (1) this for previous version
       iStart		= iStart + cSpace.cGrid.iD1TL;    
       
    end do
    
    iMemAllocated=iMemAllocated - ( PRODUCT(SHAPE(dTaperSupportWindow)) + PRODUCT(SHAPE(dTaperMaxFreqWindow)) + PRODUCT(SHAPE(dMultFactor)))* dpS
    DEALLOCATE(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor)
    call PrintToLog("Multiply with the correction factor", 2);
    
	if (ALLOCATED(ContrastFiltered) ) then	  ! Remove (2) this if for previous version
    	!! L.D. 02-12-2009 INHOMOGENEITIES IN THE NONLINEAR PARAMETER
		ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    	iMemAllocated=iMemAllocated + PRODUCT(SHAPE(ContrastFiltered)) * dpS
		ContrastFiltered=0
	    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
		istart=0
		do iindex=0,pcgrid%id1locn-1
		   
		   !!---------------------------------------------------------------
		   if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
		      pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA1*cSourceParams.BETA1* ((cMediumParams.c0**2)/(cSourceParams.SOS1**2)) ! Insert eventual correction factor Beta/BetaBG
		   elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
		      pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA2*cSourceParams.BETA2* ((cMediumParams.c0**2)/(cSourceParams.SOS2**2)) !0.000000001692909375602988 
		   elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
		      pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA3*cSourceParams.BETA3* ((cMediumParams.c0**2)/(cSourceParams.SOS3**2))!0.000000002500320333387511 
		   elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
		      pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA4*cSourceParams.BETA4* ((cMediumParams.c0**2)/(cSourceParams.SOS4**2))!0.000000001523550004645656 
		   elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
		      pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA5*cSourceParams.BETA5* ((cMediumParams.c0**2)/(cSourceParams.SOS5**2))!* 0.000000001673303658296936   * ((cMediumParams.c0**2)/(1578**2))!0.000000001673303658296936 ! Insert eventual correction factor Beta/BetaBG
		   else
			 pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
		   end if
		   !!---------------------------------------------------------------			
		   
		   istart=istart+pcgrid%id1is
		end do
		!! L.D. 02-12-2009
		
		dFFTFactor = 2.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp))
		!call PrintToLog("AttenuationEval", 2);
		
		
		b_att1=cSourceParams.b1  !! Blood
		b_att2=cSourceParams.b2  !! Brain 
		b_att3=cSourceParams.b3  !! Breast 
		b_att4=cSourceParams.b4  !! Heart
		b_att5=cSourceParams.b5  !! Liver
		alpha1=(cSourceParams.alpha1*100)/((4*half_pi*1000000)**(b_att1))        !! Blood
		alpha2=(cSourceParams.alpha2*100)/((4*half_pi*1000000)**(b_att2))        !! Brain
		alpha3=(cSourceParams.alpha3*100)/((4*half_pi*1000000)**(b_att3))        !! Breast
		alpha4=(cSourceParams.alpha4*100)/((4*half_pi*1000000)**(b_att4))        !! Heart
		alpha5=(cSourceParams.alpha5*100)/((4*half_pi*1000000)**(b_att5))        !! Liver 
		
		
		ALLOCATE(AttenContrast1(1:iDimW),AttenContrast2(1:iDimW),AttenContrast3(1:iDimW),AttenContrast4(1:iDimW),AttenContrast5(1:iDimW))
        iMemAllocated=iMemAllocated + (PRODUCT(SHAPE(AttenContrast1)) + PRODUCT(SHAPE(AttenContrast2)) + PRODUCT(SHAPE(AttenContrast3)) + PRODUCT(SHAPE(AttenContrast4)) + PRODUCT(SHAPE(AttenContrast5))) * dpcS
		
		do i=0, idimw-1
		   
		   !! ATTENUATION MODELLED ONLY VIA CONTRAST SOURCE APPROACH	
		   !	AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1);
		   !	AttenContrast2(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha2*(abs(i*dDOmega*cModelParams%freq0)**b_att2)+1+cMediumParams%c0*alpha2*tan(half_pi*b_att2)*(abs(i*dDOmega*cModelParams%freq0)**(b_att2-1)-(two_pi*cMediumParams%fc0)**(b_att2-1)))**2-1);
		   !	AttenContrast3(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha3*(abs(i*dDOmega*cModelParams%freq0)**b_att3)+1+cMediumParams%c0*alpha3*tan(half_pi*b_att3)*(abs(i*dDOmega*cModelParams%freq0)**(b_att3-1)-(two_pi*cMediumParams%fc0)**(b_att3-1)))**2-1);
		   !	AttenContrast4(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha4*(abs(i*dDOmega*cModelParams%freq0)**b_att4)+1+cMediumParams%c0*alpha4*tan(half_pi*b_att4)*(abs(i*dDOmega*cModelParams%freq0)**(b_att4-1)-(two_pi*cMediumParams%fc0)**(b_att4-1)))**2-1);
		   !	AttenContrast5(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha5*(abs(i*dDOmega*cModelParams%freq0)**b_att5)+1+cMediumParams%c0*alpha5*tan(half_pi*b_att5)*(abs(i*dDOmega*cModelParams%freq0)**(b_att5-1)-(two_pi*cMediumParams%fc0)**(b_att5-1)))**2-1);
		   
		   !! WITH LOSSY GREEN'S FUNCTION
		   AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
		   
		   AttenContrast2(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha2*(abs(i*dDOmega*cModelParams%freq0)**b_att2)+1+cMediumParams%c0*alpha2*tan(half_pi*b_att2)*(abs(i*dDOmega*cModelParams%freq0)**(b_att2-1)-(two_pi*cMediumParams%fc0)**(b_att2-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
		   
		   AttenContrast3(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha3*(abs(i*dDOmega*cModelParams%freq0)**b_att3)+1+cMediumParams%c0*alpha3*tan(half_pi*b_att3)*(abs(i*dDOmega*cModelParams%freq0)**(b_att3-1)-(two_pi*cMediumParams%fc0)**(b_att3-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
		   
		   AttenContrast4(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha4*(abs(i*dDOmega*cModelParams%freq0)**b_att4)+1+cMediumParams%c0*alpha4*tan(half_pi*b_att4)*(abs(i*dDOmega*cModelParams%freq0)**(b_att4-1)-(two_pi*cMediumParams%fc0)**(b_att4-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
		   
		   AttenContrast5(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha5*(abs(i*dDOmega*cModelParams%freq0)**b_att5)+1+cMediumParams%c0*alpha5*tan(half_pi*b_att5)*(abs(i*dDOmega*cModelParams%freq0)**(b_att5-1)-(two_pi*cMediumParams%fc0)**(b_att5-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
		   
		end do
		
		AttenContrast1(1)= 0;
		AttenContrast2(1)= 0;
		AttenContrast3(1)= 0;
		AttenContrast4(1)= 0;	
		AttenContrast5(1)= 0;	
		
		iStart = 0
		INDEXFORZ=0
		ContrastFiltered=0
	!    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009

		do iindex=0,pcgrid%id1locn-1
		   !! HEART MODEL -------------------------------
		   if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
		      acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS1**2) * AttenContrast1 - (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS1**2))))
		      
		   else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
		      acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS2**2) * AttenContrast2 - (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS2**2))))
		      
		   else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
		      acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS3**2) * AttenContrast3 - (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS3**2))))
		      
		   else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
		      acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS4**2) * AttenContrast4 - (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS4**2))))
		      
		   else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
		      acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS5**2) * AttenContrast5 - (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS5**2))))
		   else
		      acbuffer(istart+1:istart+idimw) = 0
		   end if
		   
		   pcgrid%pacd1(istart+1:istart+idimw) =  pcgrid%pacd1(istart+1:istart+idimw) - acbuffer(istart+1:istart+idimw) 
		   istart=istart+pcgrid%id1is
		   
		end do
   	 	iMemAllocated=iMemAllocated - (PRODUCT(SHAPE(AttenContrast1)) + PRODUCT(SHAPE(AttenContrast2)) + PRODUCT(SHAPE(AttenContrast3)) + PRODUCT(SHAPE(AttenContrast4)) + PRODUCT(SHAPE(AttenContrast5)) ) * dpcS
    	iMemAllocated=iMemAllocated - PRODUCT(SHAPE(ContrastFiltered)) * dpS
    	DEALLOCATE( AttenContrast1,AttenContrast2,AttenContrast3,AttenContrast4,AttenContrast5)
    	DEALLOCATE(ContrastFiltered)
	endif
    iMemAllocated=iMemAllocated - PRODUCT(SHAPE(acBuffer)) * dpcS
    DEALLOCATE(acBuffer)
    
    !Finally, transform back to T-domain with small transform, 
    ! without the anti-aliasing regions, and reorder to dist. 0
    !call PrintToLog("AttenuationEval4", 2);
    call TransformTInv_sml(pcGrid) 
    call ReorderDistr1ToDistr0(pcGrid)
  END SUBROUTINE NonlinContrastOperator_Ali
  
  !----------------
  SUBROUTINE Theta2CG_Ali(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     061009  Original code (LD)
    !
    ! *****************************************************************************
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    type(Grid), pointer :: pcGrid
    integer(i8b)::		   iIndex, iStart, i, iDimW, iDimZ, iDimX, iDimY, INDEXFORZ;
    real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT)
    real(dp) ::			   dDOmega, dFFTFactor
    real(dp), allocatable :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::		dLeftBand,dRightBand,dCorrection
    real(dp) :: alpha1,alpha2,alpha3,alpha4,alpha5,b_att1,b_att2,b_att3,b_att4,b_att5
    complex(dpc), allocatable :: acBuffer(:), acBuffer2(:),AttenContrast1(:),AttenContrast2(:),AttenContrast3(:),AttenContrast4(:),AttenContrast5(:)
    real(8), dimension(:,:,:), allocatable :: ContrastFiltered !! L.D. 03-11-2009
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv
    !   TransformT
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    ! =============================================================================
    
    call PrintToLog("NonlinOperator_Ali",4)
    
    pcGrid=>cSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    iDimZ = cSpace%iDimZ
    iDimX = cSpace%iDimX
    iDimY = cSpace%iDimY
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW),&
    acBuffer(1:pcgrid%iD1LocSize),acBuffer2(1:pcgrid%iD1LocSize),&
    AttenContrast1(1:iDimW),AttenContrast2(1:iDimW),AttenContrast3(1:iDimW),AttenContrast4(1:iDimW),AttenContrast5(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          !pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          !		pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) !No Tapering
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          !pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          !		pcGrid%pacD1(iStart+1:iStart+iDimW) ! NO TAPERING		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !!L.D. 
    ! acBuffer now contains the tapered pressure in w-domain
    acBuffer=pcGrid%pacD1
    acBuffer2=pcGrid%pacD1
    
    
    !!L.D.
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1**2 + im * arBuffer2**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex
       
    end do
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr6",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ 2*cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !Transform to W-domain and multiply with -w**2 and with 
    ! FFT normalization factor (1/N for each FORWARD transform, 
    ! and an extra square because of the square factor above;
    ! this gives (1/N)**2*(1/(2*N)) )
    call TransformT(pcGrid)
    
    dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Here, the max freq tapering is applied as well
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
       !dMultFactor = dMultFactor ! NO TAPERING
    end if
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       acBuffer2(iStart+1:iStart+iDimW)=acBuffer2(iStart+1:iStart+iDimW)* 2 *dMultFactor*real(pcGrid%iD0TL,dp)**2;;
       
       iStart		= iStart + cSpace.cGrid.iD1TL;
       
    end do
    
    call PrintToLog("Multiply with the correction factor", 2);
    dCorrection	= cMediumParams%Kappa0 * cMediumParams%Beta
    !print *,cMediumParams%Beta
    
    !do i = 1, pcGrid%iD1LocSize
    !	pcGrid%pacD1(i)	= pcGrid%pacD1(i) * dCorrection;
    !end do
    
    !! L.D. 02-12-2009 INHOMOGENEITIES IN THE NONLINEAR PARAMETER
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    ContrastFiltered=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    do iindex=0,pcgrid%id1locn-1
       !!---------------------------------------------------------------
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA1*cSourceParams.BETA1* ((cMediumParams.c0**2)/(cSourceParams.SOS1**2)) ! * ((cMediumParams.c0**2)/(1584**2)) ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA2*cSourceParams.BETA2* ((cMediumParams.c0**2)/(cSourceParams.SOS2**2)) !0.000000001692909375602988  ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA3*cSourceParams.BETA3* ((cMediumParams.c0**2)/(cSourceParams.SOS3**2))!0.000000002500320333387511 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA4*cSourceParams.BETA4* ((cMediumParams.c0**2)/(cSourceParams.SOS4**2))!0.000000001523550004645656 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) *  cSourceParams.KAPPA5*cSourceParams.BETA5* ((cMediumParams.c0**2)/(cSourceParams.SOS5**2))!* 0.000000001673303658296936   * ((cMediumParams.c0**2)/(1578**2))!0.000000001673303658296936 ! Insert eventual correction factor Beta/BetaBG
       else
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
       end if
       !!---------------------------------------------------------------		
       
       istart=istart+pcgrid%id1is
    end do
    !! L.D. 02-12-2009		
    
    !dKangular is the complex wavenumber (normalized wrt c0/f0) that includes losses
    !dKangular = cMediumParams%c0/cModelParams%freq0 * (dOmega*cModelParams%freq0/cMediumParams%c0	&
    !+ cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*dOmega*cModelParams%freq0		&
    !* (abs(dOmega*cModelParams%freq0)**(cMediumParams%b_att-1) -				&
    !(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1))						&
    !- im*cMediumParams%alpha0_att*abs(dOmega*cModelParams%freq0)**cMediumParams%b_att)
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr8",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT+1, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.false.);
    
    !Finally, transform back to T-domain with small transform, 
    ! without the anti-aliasing regions, and reorder to dist. 0
    !call PrintToLog("AttenuationEval4", 2);
    call TransformTInv_sml(pcGrid) 
    !call PrintToLog("AttenuationEval5", 2);
    call ReorderDistr1ToDistr0(pcGrid)
    !call PrintToLog("AttenuationEval6", 2);
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr10",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !L.D.
    !				call ExportSlice(trim(sOutputDir) // trim(sOutputRoot) // int2str(int(cModelParams%xyzslicepos(i) *1000* cMediumParams.c0/cModelParams.freq0 * 1.0D3)) //"ctr10",&
    !						"Contrast", cSpace, &
    !					(/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
    !					(/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
    !					.true.);
    !L.D.
    deallocate(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor, acBuffer, acBuffer2,AttenContrast1,AttenContrast2,AttenContrast3,AttenContrast4,AttenContrast5)
    
  END SUBROUTINE Theta2CG_Ali
  !----------------
  
  
!!!L.D. LCSM L.D. 12-07-2010
  SUBROUTINE NonlinContrastOperatorlin_Ali(cSpace,IncSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     061009  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperatorlin_Ali computes the nonlinearity 
    !   contrast source S^NL(p) = kappa0 * beta * d^2 / dt^2 p^2 for the given 
    !   space. The data in cSpace should be in Distribution 0. This implementation 
    !   employs an anti-aliasing approach for the temporal dimension. The temporal
    !   derivative is implemented in the Fourier domain with a spectral derivative.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    type(Space), target, intent(in) ::                IncSpace
    type(Grid), pointer :: pcGrid,pcGridStatic
    integer(i8b)::		   iIndex, iStart, i, iDimW, iDimZ, iDimX, iDimY, iDimW2, iDimZ2, iDimX2, iDimY2, INDEXFORZ;
    real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT),arBuffer1Static(cSpace%iDimT),arBuffer2Static(cSpace%iDimT)
    real(dp) ::			   dDOmega, dFFTFactor
    real(dp), allocatable :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::		dLeftBand,dRightBand,dCorrection
    real(dp) :: alpha1, b_att1, alpha2, b_att2, alpha3, b_att3, alpha4, b_att4, alpha5, b_att5
    complex(dpc), allocatable :: acBuffer(:),acBuffer2(:),  AttenContrast1(:), AttenContrast2(:), AttenContrast3(:), AttenContrast4(:), AttenContrast5(:)
    real(8), dimension(:,:,:), allocatable :: ContrastFiltered !! L.D. 03-11-2009
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv
    !   TransformT
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    ! =============================================================================
    
    call PrintToLog("NonlinOperator_Ali_Linearized",4)
    
    pcGrid=>cSpace%cGrid
    
    pcGridStatic=>IncSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    iDimZ = cSpace%iDimZ
    iDimX = cSpace%iDimX
    iDimY = cSpace%iDimY
    iDimW2 = IncSpace%iDimT/2 + 1
    iDimZ2 = IncSpace%iDimZ
    iDimX2 = size(pcGridStatic%parD0)
    iDimY2 = size(pcGrid%parD0)
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW),&
    acBuffer(1:pcgrid%iD1LocSize),acBuffer2(1:pcgrid%iD1LocSize),&
    AttenContrast1(1:iDimW),AttenContrast2(1:iDimW),AttenContrast3(1:iDimW),AttenContrast4(1:iDimW),AttenContrast5(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow	
          pcGridStatic%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGridStatic%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow			
          !acBuffer2%parD0(iStart+1:iStart+cSpace%iDimT) = &
          !		acBuffer2%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow	
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call ReorderDistr0ToDistr1(pcGridStatic)
    
    call TransformT_sml(pcGrid)
    call TransformT_sml(pcGridStatic)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          pcGridStatic%pacD1(iStart+1:iStart+iDimW) = &
          pcGridStatic%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          !acBuffer2%pacD1(iStart+1:iStart+iDimW) = &
          !		acBuffer2%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !!L.D. 
    ! acBuffer and acBuffer2 now contain the tapered pressure in w-domain
    acBuffer=pcGrid%pacD1
    acBuffer2=pcGrid%pacD1
    
    
    !!L.D.
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridStatic)
    !call TransformTInv(acBuffer2)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       arBuffer1Static=real(pcGridStatic%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2Static=dimag(pcGridStatic%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       !pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &  !! If normal Neumann
       !	arBuffer1**2 + im * arBuffer2**2;		
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &  !! If normal Neumann
            arBuffer1*arBuffer1Static + im * arBuffer2*arBuffer2Static;	
       
       iStart		= iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex
       
    end do
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr6",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ 2*cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !Transform to W-domain and multiply with -w**2 and with 
    ! FFT normalization factor (1/N for each FORWARD transform, 
    ! and an extra square because of the square factor above;
    ! this gives (1/N)**2*(1/(2*N)) )
    call TransformT(pcGrid)
    call TransformT(pcGridStatic)
    
    
    dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Here, the max freq tapering is applied as well
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
       !dMultFactor = dMultFactor ! NO TAPERING
    end if
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       acBuffer2(iStart+1:iStart+iDimW)=acBuffer2(iStart+1:iStart+iDimW)* 2 *dMultFactor*real(pcGrid%iD0TL,dp)**2;;
       
       iStart		= iStart + cSpace.cGrid.iD1TL;
       
    end do
    
    call PrintToLog("Multiply with the correction factor", 2);
    dCorrection	= 2*cMediumParams%Kappa0 * cMediumParams%Beta
    
    
    !do i = 1, pcGrid%iD1LocSize
    !	pcGrid%pacD1(i)	= pcGrid%pacD1(i) * dCorrection;
    !end do
    
    !! L.D. 02-12-2009 INHOMOGENEITIES IN THE NONLINEAR PARAMETER
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    ContrastFiltered=0
    call PrintToLog("PreMatrix", 2);
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    call PrintToLog("PostMatrix", 2);
    istart=0
    
    do iindex=0,pcgrid%id1locn-1
       !!---------------------------------------------------------------
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2*   cSourceParams.KAPPA1*cSourceParams.BETA1* ((cMediumParams.c0**2)/(cSourceParams.SOS1**2)) ! 0.0000000015040 * ((cMediumParams.c0**2)/(1584**2)) ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2*   cSourceParams.KAPPA2*cSourceParams.BETA2* ((cMediumParams.c0**2)/(cSourceParams.SOS2**2)) !0.000000001692909375602988  ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2*   cSourceParams.KAPPA3*cSourceParams.BETA3* ((cMediumParams.c0**2)/(cSourceParams.SOS3**2))!0.000000002500320333387511 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2*   cSourceParams.KAPPA4*cSourceParams.BETA4* ((cMediumParams.c0**2)/(cSourceParams.SOS4**2))!0.000000001523550004645656 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2*   cSourceParams.KAPPA5*cSourceParams.BETA5* ((cMediumParams.c0**2)/(cSourceParams.SOS5**2))!*  0.000000001673303658296936   * ((cMediumParams.c0**2)/(1578**2))!0.000000001673303658296936 ! Insert eventual correction factor Beta/BetaBG
       else
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
       end if
       !!---------------------------------------------------------------		
       
       !						if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001503985923654306 ! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001692909375602988  ! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2  * 0.000000002500320333387511 ! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001523550004645656 ! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001673303658296936 ! Insert eventual correction factor Beta/BetaBG
       !						else
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
       !						end if	
       !						!----------------------------------------------------------------------
       !**********************************************************************
       istart=istart+pcgrid%id1is
    end do
    !! L.D. 02-12-2009
    call PrintToLog("PostMatrix 2", 2);
    dFFTFactor = 2.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp))
    
    
    b_att1=cSourceParams.b1  !! Blood
    b_att2=cSourceParams.b2  !! Brain 
    b_att3=cSourceParams.b3  !! Breast 
    b_att4=cSourceParams.b4  !! Heart
    b_att5=cSourceParams.b5  !! Liver
    alpha1=(cSourceParams.alpha1*100)/((4*half_pi*1000000)**(b_att1))        !! Blood
    alpha2=(cSourceParams.alpha2*100)/((4*half_pi*1000000)**(b_att2))        !! Brain
    alpha3=(cSourceParams.alpha3*100)/((4*half_pi*1000000)**(b_att3))        !! Breast
    alpha4=(cSourceParams.alpha4*100)/((4*half_pi*1000000)**(b_att4))        !! Heart
    alpha5=(cSourceParams.alpha5*100)/((4*half_pi*1000000)**(b_att5))        !! Liver  
    
    call PrintToLog("PostMatrix 3", 2);
    do i=0, idimw-1
       !! Only contrast source technique used
       !	AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1);
       
       !! LOSSY GREEN's function included
       AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast2(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha2*(abs(i*dDOmega*cModelParams%freq0)**b_att2)+1+cMediumParams%c0*alpha2*tan(half_pi*b_att2)*(abs(i*dDOmega*cModelParams%freq0)**(b_att2-1)-(two_pi*cMediumParams%fc0)**(b_att2-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast3(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha3*(abs(i*dDOmega*cModelParams%freq0)**b_att3)+1+cMediumParams%c0*alpha3*tan(half_pi*b_att3)*(abs(i*dDOmega*cModelParams%freq0)**(b_att3-1)-(two_pi*cMediumParams%fc0)**(b_att3-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast4(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha4*(abs(i*dDOmega*cModelParams%freq0)**b_att4)+1+cMediumParams%c0*alpha4*tan(half_pi*b_att4)*(abs(i*dDOmega*cModelParams%freq0)**(b_att4-1)-(two_pi*cMediumParams%fc0)**(b_att4-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast5(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha5*(abs(i*dDOmega*cModelParams%freq0)**b_att5)+1+cMediumParams%c0*alpha5*tan(half_pi*b_att5)*(abs(i*dDOmega*cModelParams%freq0)**(b_att5-1)-(two_pi*cMediumParams%fc0)**(b_att5-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       
    end do
    call PrintToLog("PostMatrix 3", 2);
    AttenContrast1(1)= 0;
    AttenContrast2(1)= 0;
    AttenContrast3(1)= 0;
    AttenContrast4(1)= 0;	
    AttenContrast5(1)= 0;		
    
    iStart = 0
    INDEXFORZ=0
    ContrastFiltered=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    
    do iindex=0,pcgrid%id1locn-1
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then					
          acbuffer(istart+1:istart+idimw) = &
          acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS1**2)) * attencontrast1
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS1**2)))
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then	
          acbuffer(istart+1:istart+idimw) = &
          acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS2**2)) * attencontrast2
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS2**2)))
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then	
          acbuffer(istart+1:istart+idimw) = &
          acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS3**2)) * attencontrast3
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS3**2)))
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then	
          acbuffer(istart+1:istart+idimw) = &
          acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS4**2)) * attencontrast4
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS4**2))) !&
          !acbuffer2(istart+1:istart+idimw) * -( i * dDOmega )**2*dFFTFactor* 1/2 * ((1500**2-cMediumParams.c0**2)/(1500**2))
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then	
          acbuffer(istart+1:istart+idimw) = &
          acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS5**2)) * attencontrast5
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS5**2)))
          !acbuffer2(istart+1:istart+idimw) * -( i * dDOmega )**2*dFFTFactor* 1/2 * ((1500**2-cMediumParams.c0**2)/(1500**2))
       else
          acbuffer(istart+1:istart+idimw)  = 0
          acbuffer2(istart+1:istart+idimw) = 0
       endif
       !acbuffer(istart+1:istart+idimw) * attencontrast
       
       pcgrid%pacd1(istart+1:istart+idimw) = &
       pcgrid%pacd1(istart+1:istart+idimw) - acbuffer(istart+1:istart+idimw) + acbuffer2(istart+1:istart+idimw)	
       istart=istart+pcgrid%id1is
       
    end do
    
    !dKangular is the complex wavenumber (normalized wrt c0/f0) that includes losses
    !dKangular = cMediumParams%c0/cModelParams%freq0 * (dOmega*cModelParams%freq0/cMediumParams%c0	&
    !+ cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*dOmega*cModelParams%freq0		&
    !* (abs(dOmega*cModelParams%freq0)**(cMediumParams%b_att-1) -				&
    !(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1))						&
    !- im*cMediumParams%alpha0_att*abs(dOmega*cModelParams%freq0)**cMediumParams%b_att)
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr8",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT+1, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.false.);
    
    !Finally, transform back to T-domain with small transform, 
    ! without the anti-aliasing regions, and reorder to dist. 0
    !call PrintToLog("AttenuationEval4", 2);
    call TransformTInv_sml(pcGrid) 
    call TransformTInv_sml(pcGridStatic) 
    !call TransformTInv_sml(acbuffer2)
    
    !call PrintToLog("AttenuationEval5", 2);
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridStatic)
    !call ReorderDistr1ToDistr0(acbuffer2)
    !call PrintToLog("AttenuationEval6", 2);
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr10",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !L.D.
    !				call ExportSlice(trim(sOutputDir) // trim(sOutputRoot) // int2str(int(cModelParams%xyzslicepos(i) *1000* cMediumParams.c0/cModelParams.freq0 * 1.0D3)) //"ctr10",&
    !						"Contrast", cSpace, &
    !					(/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
    !					(/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
    !					.true.);
    !L.D.
    deallocate(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor, acBuffer, acBuffer2, AttenContrast1, AttenContrast2, AttenContrast3, AttenContrast4, AttenContrast5)
    
  END SUBROUTINE NonlinContrastOperatorlin_Ali
!!!!L.D.
  
  !============================================================C.G.===L.D. 10-08-2011============
  SUBROUTINE NonlinContrastOperatorCG_Ali(cSpace)
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     061009  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperatorCG_Ali computes the nonlinearity 
    !   contrast source S^NL(p) = kappa0 * beta * d^2 / dt^2 p^2 for the given 
    !   space. The data in cSpace should be in Distribution 0. This implementation 
    !   employs an anti-aliasing approach for the temporal dimension. The temporal
    !   derivative is implemented in the Fourier domain with a spectral derivative.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    
    type(Grid), pointer :: pcGrid
    integer(i8b)::		   iIndex, iStart, i, iDimW, iDimZ, iDimX, iDimY, iDimW2, iDimZ2, iDimX2, iDimY2, INDEXFORZ;
    real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT)
    real(dp) ::			   dDOmega, dFFTFactor
    real(dp), allocatable :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::		dLeftBand,dRightBand,dCorrection
    real(dp) :: alpha1, b_att1, alpha2, b_att2
    complex(dpc), allocatable :: acBuffer(:), acBuffer2(:), AttenContrast1(:), AttenContrast2(:)
    real(8), dimension(:,:,:), allocatable :: ContrastFiltered !! L.D. 03-11-2009
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv
    !   TransformT
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   If required, apply the tapering window to the field that removes the sharp
    !     edge at start and end of the time trace
    !   Reorder to distribution 1
    !   Transform the data in the time dimension - use 'small' transform (iDimT)
    !   If required, apply the frequency tapering window that removes the highest
    !     frequencies - these may be unrealistically large in the iterative scheme
    !   Inverse transform with a 'large' transform (2*iDimT) that includes the
    !     anti-aliasing frequency range
    !   Apply the square p**2
    !   Transform to the frequency domain with a 'large' transform (2*iDimT)
    !   If required, again apply the frequency tapering window
    !   Apply the FFT factors and the derivative term -omega**2
    !   Inverse transform with a 'small' transform (iDimT)
    !   Reorder to distribution 0
    !   Apply the factor kappa_0 * beta in the contrast source
    !
    ! =============================================================================
    
    call PrintToLog("NonlinOperator_Ali_Linearized",4)
    
    pcGrid=>cSpace%cGrid
    
    
    iDimW = cSpace%iDimT/2 + 1
    iDimZ = cSpace%iDimZ
    iDimX = cSpace%iDimX
    iDimY = cSpace%iDimY
    
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW),&
    acBuffer(1:pcgrid%iD1LocSize),acBuffer2(1:pcgrid%iD1LocSize),&
    AttenContrast1(1:iDimW),AttenContrast2(1:iDimW))
    !print *, "blocco numero 1"
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow			
          !pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          !		pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) !No Tapering
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
       
    end if
    !print *, "blocco numero 2"
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    call TransformT_sml(pcGrid)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          !pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          !		pcGrid%pacD1(iStart+1:iStart+iDimW) ! NO TAPERING		
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !!L.D. 
    ! acBuffer and acBuffer2 now contain the tapered pressure in w-domain
    acBuffer=pcGrid%pacD1
    acBuffer2=pcGrid%pacD1
    
    
    !!L.D.
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       
       !pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &  !! If normal Neumann
       !	arBuffer1**2 + im * arBuffer2**2;		
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &  !! If normal Neumann
            arBuffer1 + im * arBuffer2;	
       
       iStart		= iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex
       
    end do
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr6",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ 2*cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !Transform to W-domain and multiply with -w**2 and with 
    ! FFT normalization factor (1/N for each FORWARD transform, 
    ! and an extra square because of the square factor above;
    ! this gives (1/N)**2*(1/(2*N)) )
    call TransformT(pcGrid)
    
    
    dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Here, the max freq tapering is applied as well
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
       !dMultFactor = dMultFactor ! NO TAPERING
    end if
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       acBuffer2(iStart+1:iStart+iDimW)=acBuffer2(iStart+1:iStart+iDimW)* 2 *dMultFactor*real(pcGrid%iD0TL,dp)**2;;
       
       iStart		= iStart + cSpace.cGrid.iD1TL;
       
    end do
    
    call PrintToLog("Multiply with the correction factor", 2);
    dCorrection	= (2*cMediumParams%Beta)/( ((cMediumParams%c0)**4)*cMediumParams.rho0 )
    !print *,cMediumParams%Beta
    
    !do i = 1, pcGrid%iD1LocSize
    !	pcGrid%pacD1(i)	= pcGrid%pacD1(i) * dCorrection;
    !end do
    
    !! L.D. 02-12-2009 INHOMOGENEITIES IN THE NONLINEAR PARAMETER
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    ContrastFiltered=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    do iindex=0,pcgrid%id1locn-1
       
       
       !!---------------------------------------------------------------
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * cSourceParams.KAPPA1*cSourceParams.BETA1* (1.0_dp/(cSourceParams.SOS1**2))!*  0.000000001503985923654306 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * cSourceParams.KAPPA2*cSourceParams.BETA2* (1.0_dp/(cSourceParams.SOS2**2))  !0.000000001692909375602988  ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * cSourceParams.KAPPA3*cSourceParams.BETA3* (1.0_dp/(cSourceParams.SOS3**2))!0.000000002500320333387511 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * cSourceParams.KAPPA4*cSourceParams.BETA4* (1.0_dp/(cSourceParams.SOS4**2)) !0.000000001523550004645656 ! Insert eventual correction factor Beta/BetaBG
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * cSourceParams.KAPPA5*cSourceParams.BETA5* (1.0_dp/(cSourceParams.SOS5**2)) !0.000000001673303658296936 !0.000000001673303658296936 ! Insert eventual correction factor Beta/BetaBG
       else
          pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
       end if
       !!---------------------------------------------------------------
       
       !						if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001503985923654306 * (1.0_dp/(cMediumParams%c0**2))! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 *0.000000001692909375602988 * (1.0_dp/(cMediumParams%c0**2)) ! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000002500320333387511 * (1.0_dp/(cMediumParams%c0**2))! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001523550004645656 * (1.0_dp/(cMediumParams%c0**2))! Insert eventual correction factor Beta/BetaBG
       !						elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * 2 * 0.000000001673303658296936 * (1.0_dp/(cMediumParams%c0**2))! Insert eventual correction factor Beta/BetaBG
       !						else
       !						pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
       !						end if	
       
       istart=istart+pcgrid%id1is
    end do
    !! L.D. 02-12-2009    
    
    
    !dKangular is the complex wavenumber (normalized wrt c0/f0) that includes losses
    !dKangular = cMediumParams%c0/cModelParams%freq0 * (dOmega*cModelParams%freq0/cMediumParams%c0	&
    !+ cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*dOmega*cModelParams%freq0		&
    !* (abs(dOmega*cModelParams%freq0)**(cMediumParams%b_att-1) -				&
    !(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1))						&
    !- im*cMediumParams%alpha0_att*abs(dOmega*cModelParams%freq0)**cMediumParams%b_att)
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr8",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT+1, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.false.);
    
    !Finally, transform back to T-domain with small transform, 
    ! without the anti-aliasing regions, and reorder to dist. 0
    !call PrintToLog("AttenuationEval4", 2);
    call TransformTInv_sml(pcGrid) 
    
    !call PrintToLog("AttenuationEval5", 2);
    call ReorderDistr1ToDistr0(pcGrid)
    
    !call PrintToLog("AttenuationEval6", 2);
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "ctr10",&
    !						"Contrast", &
    !						cSpace, &
    !						(/ 0_i8b, 45_i8b, 0_i8b, 0_i8b /), &
    !						(/ cSpace%iDimT, 1_i8b, 1_i8b, cSpace%iDimZ /), &
    !						.true.);
    
    !L.D.
    !				call ExportSlice(trim(sOutputDir) // trim(sOutputRoot) // int2str(int(cModelParams%xyzslicepos(i) *1000* cMediumParams.c0/cModelParams.freq0 * 1.0D3)) //"ctr10",&
    !						"Contrast", cSpace, &
    !					(/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
    !					(/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
    !					.true.);
    !L.D.
    deallocate(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor, acBuffer, acBuffer2, AttenContrast1, AttenContrast2)
    
  END SUBROUTINE NonlinContrastOperatorCG_Ali
  !============================================================C.G.===L.D. 10-08-2011============
  
  !============================================================C.G.===L.D. 10-08-2011============
  SUBROUTINE ConversionCG_Ali(cSpace)
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     061009  Original code (LD)
    !
    ! *****************************************************************************
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the nonlinearity contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::		cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   i             i8b   loop counter
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    type(Grid), pointer :: pcGrid
    integer(i8b)::		   iIndex, iStart, i, iDimW, iDimZ, iDimX, iDimY, INDEXFORZ;
    real(dp) ::			   arBuffer1(cSpace%iDimT),arBuffer2(cSpace%iDimT)
    real(dp) ::			   dDOmega, dFFTFactor
    real(dp), allocatable :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::		dLeftBand,dRightBand,dCorrection
    real(dp) :: alpha1,alpha2,alpha3,alpha4,alpha5,b_att1,b_att2,b_att3,b_att4,b_att5
    complex(dpc), allocatable :: acBuffer(:), acBuffer2(:), AttenContrast1(:),AttenContrast2(:),AttenContrast3(:),AttenContrast4(:),AttenContrast5(:)
    real(8), dimension(:,:,:), allocatable :: ContrastFiltered !! L.D. 03-11-2009
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv
    !   TransformT
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    !
    ! =============================================================================
    
    call PrintToLog("NonlinOperator_Ali",4)
    
    pcGrid=>cSpace%cGrid
    !acBuffer2=>cSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    iDimZ = cSpace%iDimZ
    iDimX = cSpace%iDimX
    iDimY = cSpace%iDimY
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW),&
    acBuffer(1:pcgrid%iD1LocSize),acBuffer2(1:pcgrid%iD1LocSize),&
    AttenContrast1(1:iDimW),AttenContrast2(1:iDimW),AttenContrast3(1:iDimW),AttenContrast4(1:iDimW),AttenContrast5(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          !acBuffer2%parD0(iStart+1:iStart+cSpace%iDimT) = &
          !		acBuffer2%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    !	call ReorderDistr0ToDistr1(acBuffer2)
    
    call TransformT_sml(pcGrid)
    !	call TransformT_sml(acBuffer2)
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow
          !acBuffer2%pacD1(iStart+1:iStart+iDimW) = &
          !		acBuffer2%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow	
          
          iStart=iStart+pcGrid%iD1IS
       end do
       
    end if
    
    !!L.D. 
    ! acBuffer now contains the tapered pressure in w-domain
    acBuffer=pcGrid%pacD1
    acBuffer2=pcGrid%pacD1
    
    
    !!L.D.
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    call TransformTInv(pcGrid)
    !	call TransformTInv(acBuffer2)
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1 + im * arBuffer2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS; ! Even if Time domain D1 Because REal and Complex
       
    end do
    
    call TransformT(pcGrid)
    !	call TransformT(acBuffer2)
    
    
    dFFTFactor = 1.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp)**3)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Here, the max freq tapering is applied as well
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
       !dMultFactor = dMultFactor ! NO TAPERING
    end if
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       acBuffer2(iStart+1:iStart+iDimW)=acBuffer2(iStart+1:iStart+iDimW)* 2 *dMultFactor*real(pcGrid%iD0TL,dp)**2;
       
       iStart		= iStart + cSpace.cGrid.iD1TL;
       
    end do
    
    call PrintToLog("Multiply with the correction factor", 2);
    dCorrection	= 1.0_dp/(cMediumParams%c0**2)
    
    !! L.D. 02-12-2009 INHOMOGENEITIES 
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    ContrastFiltered=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    do iindex=0,pcgrid%id1locn-1
       !! HEART MODEL INTRODUCTION
       !!---------------------------------------------------------------
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS1**2) ! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS1**2) ! Insert eventual correction factor Beta/BetaBG
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS2**2)  ! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS2**2)  ! Insert eventual correction factor Beta/BetaBG
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS3**2) ! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS3**2) ! Insert eventual correction factor Beta/BetaBG
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS4**2)! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS4**2)! Insert eventual correction factor Beta/BetaBG
          
       elseif (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS5**2) ! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * 1.0_dp/(cSourceParams.SOS5**2) ! Insert eventual correction factor Beta/BetaBG
          
       else
          acBuffer(istart+1:istart+idimw)= acBuffer(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
          acBuffer2(istart+1:istart+idimw)= acBuffer2(istart+1:istart+idimw) * dCorrection     ! Insert eventual correction factor Beta/BetaBG
          
       end if
       !!---------------------------------------------------------------			
       !pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection *(contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,3)+1))
       !pcGrid%pacD1(istart+1:istart+idimw)= pcGrid%pacD1(istart+1:istart+idimw) * dCorrection !Homogeneous Case
       istart=istart+pcgrid%id1is
    end do
    !! L.D. 02-12-2009
    
    dFFTFactor = 2.0_dp/(2.0_dp*real(pcGrid%iD0TL,dp))
    !call PrintToLog("AttenuationEval", 2);
    
    b_att1=cSourceParams.b1  !! Blood
    b_att2=cSourceParams.b2  !! Brain 
    b_att3=cSourceParams.b3  !! Breast 
    b_att4=cSourceParams.b4  !! Heart
    b_att5=cSourceParams.b5  !! Liver
    alpha1=(cSourceParams.alpha1*100)/((4*half_pi*1000000)**(b_att1))        !! Blood
    alpha2=(cSourceParams.alpha2*100)/((4*half_pi*1000000)**(b_att2))        !! Brain
    alpha3=(cSourceParams.alpha3*100)/((4*half_pi*1000000)**(b_att3))        !! Breast
    alpha4=(cSourceParams.alpha4*100)/((4*half_pi*1000000)**(b_att4))        !! Heart
    alpha5=(cSourceParams.alpha5*100)/((4*half_pi*1000000)**(b_att5))        !! Liver 
    
    do i=0, idimw-1
       
       !! ATTENUATION MODELLED ONLY VIA CONTRAST SOURCE APPROACH	
       !	AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1);
       !	AttenContrast2(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha2*(abs(i*dDOmega*cModelParams%freq0)**b_att2)+1+cMediumParams%c0*alpha2*tan(half_pi*b_att2)*(abs(i*dDOmega*cModelParams%freq0)**(b_att2-1)-(two_pi*cMediumParams%fc0)**(b_att2-1)))**2-1);
       !	AttenContrast3(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha3*(abs(i*dDOmega*cModelParams%freq0)**b_att3)+1+cMediumParams%c0*alpha3*tan(half_pi*b_att3)*(abs(i*dDOmega*cModelParams%freq0)**(b_att3-1)-(two_pi*cMediumParams%fc0)**(b_att3-1)))**2-1);
       !	AttenContrast4(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha4*(abs(i*dDOmega*cModelParams%freq0)**b_att4)+1+cMediumParams%c0*alpha4*tan(half_pi*b_att4)*(abs(i*dDOmega*cModelParams%freq0)**(b_att4-1)-(two_pi*cMediumParams%fc0)**(b_att4-1)))**2-1);
       !	AttenContrast5(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha5*(abs(i*dDOmega*cModelParams%freq0)**b_att5)+1+cMediumParams%c0*alpha5*tan(half_pi*b_att5)*(abs(i*dDOmega*cModelParams%freq0)**(b_att5-1)-(two_pi*cMediumParams%fc0)**(b_att5-1)))**2-1);
       
       !! WITH LOSSY GREEN'S FUNCTION
       AttenContrast1(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha1*(abs(i*dDOmega*cModelParams%freq0)**b_att1)+1+cMediumParams%c0*alpha1*tan(half_pi*b_att1)*(abs(i*dDOmega*cModelParams%freq0)**(b_att1-1)-(two_pi*cMediumParams%fc0)**(b_att1-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast2(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha2*(abs(i*dDOmega*cModelParams%freq0)**b_att2)+1+cMediumParams%c0*alpha2*tan(half_pi*b_att2)*(abs(i*dDOmega*cModelParams%freq0)**(b_att2-1)-(two_pi*cMediumParams%fc0)**(b_att2-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast3(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha3*(abs(i*dDOmega*cModelParams%freq0)**b_att3)+1+cMediumParams%c0*alpha3*tan(half_pi*b_att3)*(abs(i*dDOmega*cModelParams%freq0)**(b_att3-1)-(two_pi*cMediumParams%fc0)**(b_att3-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast4(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha4*(abs(i*dDOmega*cModelParams%freq0)**b_att4)+1+cMediumParams%c0*alpha4*tan(half_pi*b_att4)*(abs(i*dDOmega*cModelParams%freq0)**(b_att4-1)-(two_pi*cMediumParams%fc0)**(b_att4-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       AttenContrast5(1+i)= -( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*alpha5*(abs(i*dDOmega*cModelParams%freq0)**b_att5)+1+cMediumParams%c0*alpha5*tan(half_pi*b_att5)*(abs(i*dDOmega*cModelParams%freq0)**(b_att5-1)-(two_pi*cMediumParams%fc0)**(b_att5-1)))**2-1)-(-( i * dDOmega )**2*dFFTFactor*(((cMediumParams%c0/(i*im*dDOmega*cModelParams%freq0))*cMediumParams%alpha0_att*(abs(i*dDOmega*cModelParams%freq0)**cMediumParams%b_att)+1+cMediumParams%c0*cMediumParams%alpha0_att*tan(half_pi*cMediumParams%b_att)*(abs(i*dDOmega*cModelParams%freq0)**(cMediumParams%b_att-1)-(two_pi*cMediumParams%fc0)**(cMediumParams%b_att-1)))**2-1));
       
    end do
    
    AttenContrast1(1)= 0;
    AttenContrast2(1)= 0;
    AttenContrast3(1)= 0;
    AttenContrast4(1)= 0;	
    AttenContrast5(1)= 0;	
    
    iStart = 0
    INDEXFORZ=0
    
    ContrastFiltered=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    
    do iindex=0,pcgrid%id1locn-1
       !! HEART MODEL -------------------------------
       if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==1) then
          acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS1**2)) * CONJG(attencontrast1)
          !acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(1584**2)) * (attencontrast1)
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS1**2))) !
          
       else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==2) then
          acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS2**2)) * CONJG(attencontrast2)
          !acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(1562**2)) * (attencontrast2) 
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS2**2))) !
          
       else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==3) then
          acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS3**2)) * CONJG(attencontrast3)
          !acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(1510**2)) * (attencontrast3)
          acBuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS3**2))) !
          
       else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==4) then
          acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS4**2)) * CONJG(attencontrast4) 
          !acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(1554**2)) * (attencontrast4)
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS4**2))) ! !acbuffer2(istart+1:istart+idimw) * -( i * dDOmega )**2*dFFTFactor * 1/2 * ((1482**2-cMediumParams.c0**2)/(1482**2))
          
       else if (contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,2)+1,pcgrid%aid1loc(iindex+1,3)+1)==5) then
          acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(cSourceParams.SOS5**2)) * CONJG(attencontrast5) 
          !acbuffer(istart+1:istart+idimw) = acbuffer(istart+1:istart+idimw) * ((cMediumParams.c0**2)/(1578**2)) * (attencontrast5) 
          acbuffer2(istart+1:istart+idimw) = acBuffer2(istart+1:istart+idimw)* (1.0_dp-((cMediumParams.c0**2)/(cSourceParams.SOS5**2))) !acbuffer2(istart+1:istart+idimw) * -( i * dDOmega )**2*dFFTFactor * 1/2 * ((1482**2-cMediumParams.c0**2)/(1482**2))
          
       else
          acbuffer(istart+1:istart+idimw) = 0
          acBuffer2(istart+1:istart+idimw) = 0
       end if
       !!  --------------------------------------------		
       
       !acbuffer(istart+1:istart+idimw) = &
       !acbuffer(istart+1:istart+idimw) * attencontrast1*contrastfiltered(pcgrid%aid1loc(iindex+1,1)+1,pcgrid%aid1loc(iindex+1,3)+1)
       !acbuffer(istart+1:istart+idimw) * attencontrast
       
       pcgrid%pacd1(istart+1:istart+idimw) = &
       -acbuffer(istart+1:istart+idimw) + acBuffer2(istart+1:istart+idimw)	
       istart=istart+pcgrid%id1is
       
    end do
    
    
    !call PrintToLog("AttenuationEval4", 2);
    call TransformTInv_sml(pcGrid) 
    !call TransformTInv_sml(acBuffer2) 
    !call PrintToLog("AttenuationEval5", 2);
    call ReorderDistr1ToDistr0(pcGrid)
    !call ReorderDistr1ToDistr0(acBuffer2)
    !call PrintToLog("AttenuationEval6", 2);
    
    deallocate(dTaperSupportWindow,dTaperMaxFreqWindow, dMultFactor, acBuffer, acBuffer2, AttenContrast1,AttenContrast2,AttenContrast3,AttenContrast4,AttenContrast5)
    
  END SUBROUTINE ConversionCG_Ali
  !============================================================C.G.===L.D. 10-08-2011============
  
  SUBROUTINE InhomContrastOperator(cSpace,cInhomContrast)
    
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
    !   The subroutine InhomContrastOperator computes the inhomogeneity contrast
    !   source S^C(p) = ( 1 - c_0**2/c**2 ) d^2 / dt^2 p for the given space.
    !   The data in cSpace should be in Distribution 0. This implementation does 
    !   not prevent aliasing in the spatial dimension. The 2nd order temporal 
    !   derivative is implemented with a spectral derivative.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace          io   type(space)   Space for which the inhomogeneity 
    !                                      contrast source is determined
    !   cInhomContrast   i   type(space)   Inhomogeneity contrast space
    !
    type(Space), intent(inout), target ::	cSpace
    type(Space), intent(in), target ::		cInhomContrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid    type(Grid)  pointer to the grid in cSpace
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   iLen          i8b   length of each time trace
    !   i             i8b   loop counter
    !   iContrastIndex i8b  index in the inhomogeneity contrast space 
    !   iDimW         i8b   number of points in the angular frequency domain
    !   arBuffer1     dp    Temporary buffer containing a time/frequency trace
    !   arBuffer2     dp    Temporary buffer containing a time/frequency trace
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !
    type(Grid), pointer :: pcGrid
    integer(i8b)::		iIndex, iLen, iStart, i, iContrastIndex, iDimW
    real(dp)::			arBuffer1(cSpace.iDimT),arBuffer2(cSpace.iDimT);
    real(dp) ::			dDOmega, dFFTFactor
    real(dp), allocatable ::dTaperSupportWindow(:), dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::			dLeftBand,dRightBand,dCorrection
    
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
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   If required, apply the tapering window to the field that removes the sharp
    !     edge at start and end of the time trace
    !   Reorder to distribution 1
    !   Transform the data in the time dimension - use 'small' transform (iDimT)
    !   If required, apply the frequency tapering window that removes the highest
    !     frequencies - these may be unrealistically large in the iterative scheme
    !   Apply the FFT factors and the derivative term -omega**2
    !   Inverse transform with a 'small' transform (iDimT)
    !   Reorder to distribution 0
    !   Apply the inhomogeneity contrast factor ( 1 - c_0**2/c**2 )
    !
    ! =============================================================================
    
    call PrintToLog("InhomContrastOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW))
    
    !First, calculate dt^2 p with tapered spectral derivative
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    dFFTFactor = 1.0_dp/real(pcGrid%iD0TL,dp)
    dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(cSpace%iDimT,dp)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
    end if
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       iStart		= iStart + cSpace%cGrid%iD1TL;
       
    end do
    
    !Finally, transform back to T-domain with small transform, 
    ! and reorder to dist. 0
    call TransformTInv_sml(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    !Next, we multiply with the inhomogeneity contrast function
    iStart	= 0;
    iLen	= cSpace.iDimT;
    
    ! Now we process t-block, loop over each XYZ-position
    do iIndex = 0, pcGrid%iD0LocN-1
       
       iContrastIndex = 1 + iIndex*cInhomContrast%cGrid%iD0IS
       
       pcGrid.parD0(iStart+1: iStart+iLen)= cInhomContrast%cGrid%parD0(iContrastIndex) * &
       pcGrid.parD0(iStart+1: iStart+iLen);
       
       iStart	= iStart + pcGrid%iD0TL;
       
    end do
    
  END SUBROUTINE InhomContrastOperator
  
  SUBROUTINE InhomContrastOperator_Ali(cSpace)
    
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
    !   The subroutine InhomContrastOperator_Ali computes the inhomogeneity contrast
    !   source S^C(p) = ( 1 - c_0**2/c**2 ) d^2 / dt^2 p for the given space.
    !   The data in cSpace should be in Distribution 0. This implementation employs
    !   an anti-aliasing approach for the spatial dimensions. The temporal
    !   derivative is implemented in the Fourier domain with a spectral derivative.
    
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  Space for which the inhomogeneity contrast 
    !                              source is determined
    !
    type(Space), intent(inout), target ::	cSpace
    
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid          type(Grid)   pointer to the grid in cSpace
    !   cContrastSmall  type(Space)  Space with the inhomogeneity contrast
    !   cXYZslice       type(Space)  XYZ slice space for the multiplication
    !                                with the contrast per frequency 
    !   cContrastslice  type(Space)  Contrast slice for the multiplication
    !                                per frequency
    !   iIndex        i8b   loop counter over the positions in the space
    !   iStart        i8b   start of each time trace
    !   iLen          i8b   length of each time trace
    !   i             i8b   loop counter
    !   iContrastIndex i8b  index in the inhomogeneity contrast space 
    !   iIndexX       i8b   beam index in the X-dimension
    !   iIndexY       i8b   beam index in the Y-dimension
    !   iIndexZ       i8b   beam index in the Z-dimension
    !   iOmega        i8b   index in the angular frequency dimension
    !   iDimW         i8b   number of points in the angular frequency domain
    !   dDOmega       dp    Step size for the frequency axis
    !   dFFTFactor    dp    Correction factor from the FFT's
    !   dTaperSupportWindow  dp  Temporal window with tapering at beginning and
    !                            end
    !   dTaperMaxFreqWindow  dp  Frequency window with tapering of the highest
    !                            frequency components
    !   dMultFactor   dp    Multiplication factor accounting for the derivative
    !                       term -omega**2 in the frequency domain
    !   dLeftBand     dp    Size of the tapering window at the beginning
    !   dRightBand    dp    Size of the tapering window at the end
    !   dCorrection   dp    Correction factor of the medium parameters 
    !   acTemp       char   temporary char array for output log messages
    !
    type(Grid), pointer :: pcGrid
    type(Space) ::		cXYZslice,cContrastSmall,cContrastslice
    integer(i8b)::		iIndex, iLen, iStart, i, iContrastIndex;
    integer(i8b) ::		iIndexX, iIndexY, iIndexZ, iOmega, iDimW
    real(dp) ::			dDOmega, dFFTFactor
    real(dp), allocatable ::dTaperSupportWindow(:), dTaperMaxFreqWindow(:), dMultFactor(:)
    real(dp) ::			dLeftBand,dRightBand,dCorrection
    character(len=1024)::       acTemp
    logical(lgt),parameter 		:: false =.false.
    
    ! *****************************************************************************
    !
    !   I/O
    !
    !   log file entries, and if necessary load the inhomogeneity contrast from
    !   temporary file
    !   
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   PrintToLog
    !   dTaperingWindow
    !   ReorderDistr0ToDistr1
    !   TransformT_sml
    !   SpreadFrequencies
    !   ReorderDistr1ToDistr2
    !   InitSpace
    !   InitGrid
    !   GridDistr2CreateEmpty
    !   LoadContrast
    !   DestructSpace
    !   Distr2ObtainXYZBlock
    !   ConvertSmalltoAliasingBlock
    !   TransformXYZ
    !   TransformXYZInv
    !   TaperXYZslice
    !   mapVx_Small
    !   mapVy_Small
    !   TransformXYZ_sml
    !   ConvertAliasingtoSmallBlock
    !   TransformXYZInv_sml
    !   mapVxInv_Small
    !   mapVyInv_Small
    !   Distr2PutXYZBlock
    !   ReorderDistr2ToDistr1
    !   CollectFrequencies
    !   TransformTInv_sml
    !   ReorderDistr1ToDistr0
    !   
    ! *****************************************************************************
    !
    !   PSEUDOCODE
    !
    !   If required, apply the tapering window to the field that removes the sharp
    !     edge at start and end of the time trace
    !   Reorder the field space to distribution 1
    !   Transform the field space in the time dimension - use 'small' 
    !     transform (iDimT)
    !   If required, apply the frequency tapering window that removes the highest
    !     frequencies - these may be unrealistically large in the iterative scheme
    !   Apply the FFT factors and the derivative term -omega**2
    !   Spread the frequencies from the 'small' frequency axis [(iDimT+1)/2] over the 
    !     entire frequency axis (iDimT+1)
    !   Reorder to distribution 2
    !   Initialize an XYZ field slice for the multiplication with the inhomogeneity
    !     contrast per frequency
    !   Initialize an XYZ contrast slice for the multiplication per frequency
    !   Initialize the inhomogeneity contrast space
    !   Obtain the inhomogeneity contrast from temporary file; the contrast is 
    !     given in the spatial Fourier domain
    !   Copy the contrast to the contrast slice including anti-aliasing regions in 
    !     all spatial dimensions
    !   Inverse transform the contrast slice to the spatial domain 
    !   De-initialize the inhomogeneity contrast space
    !   Apply FFT correction factors due to the transforms
    !   For each frequency do
    !       Obtain the XYZ field slice for this frequency from the field space
    !       If required, apply tapering in all spatial dimensions
    !       Map XYZ field slice in the x- and y-dimensions - use 'small' maps
    !       Transform XYZ field slice to spatial Fourier domain - use 'small' 
    !         transform (iDimX*iDimY*iDimZ)
    !       Convert 'small' spatial frequency axes to 'large' spatial frequency
    !         axes including the anti-aliasing regions
    !       Inverse transform XYZ-slice with a 'large' transform (8*iDimX*iDimY
    !         *iDimZ)
    !       Multiply XYZ field slice with XYZ inhomogeneity contrast slice
    !       Transform XYZ field slice to the spatial Fourier domain - use 'large'
    !         transform (8*iDimX*iDimY*iDimZ)
    !       Convert 'large' XYZ-slice with anti-aliasing regions to 'small' XYZ-
    !         slice
    !       Inverse transform XYZ-slice with a 'small' transform 
    !         (iDimX*iDimY*iDimZ)
    !       Inverse Map XYZ field slice in the x- and y-dimensions - 'small' maps
    !       Store XYZ field slice for this frequency in the field space
    !   De-initialize slices
    !   Reorder the field space to distribution 1 (field space now contains the
    !     contrast source)
    !   Collect the frequencies from the entire frequency axis (iDimT+1) to
    !     the 'small' frequency axis [(iDimT+1)/2]
    !   Inverse transform the field space in the time dimension - use 'small'
    !     transform (iDimT) 
    !   Reorder field space to distribution 0
    !
    ! =============================================================================
    
    call PrintToLog("InhomContrastOperator_Ali",4)
    
    pcGrid=>cSpace%cGrid
    
    iDimW = cSpace%iDimT/2 + 1
    
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW),dMultFactor(1:iDimW))
    
    !First, calculate dt^2 p with tapered spectral derivative
    
    !If requested, apply a support tapering to the field
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    !First, transform to W-domain (use small transform)
    call ReorderDistr0ToDistr1(pcGrid)
    call TransformT_sml(pcGrid)
    
    dFFTFactor = 1.0_dp/real(pcGrid%iD0TL,dp)
    dDOmega = two_pi / (cSpace%dDt * cSpace%iDimT)
    dMultFactor = - ( (/ (i,i=0,iDimW-1) /) * dDOmega )**2 * dFFTFactor
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       dMultFactor = dMultFactor * dTaperMaxFreqWindow
    end if
    
    !Now, perform the d^2/dt^2 operation in W-domain 
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       pcGrid%pacD1(iStart+1:iStart+iDimW) = &
       pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor;		
       
       iStart		= iStart + cSpace%cGrid%iD1TL;
       
    end do
    
    ! Next, apply the contrast employing an anti-aliasing region in X,Y and Z
    
    ! Spread the Frequencies in such a way that evaluation is efficient in the 
    ! XYZ-local distribution - See the header of SpreadFrequencies for more
    ! explanation.
    call SpreadFrequencies(pcGrid)
    
    ! redistribute to Dist 2
    call ReorderDistr1toDistr2(pcGrid)
    
    ! define XYZslice (one Omega value, oversampled XYZ axes)
    call InitSpace(cXYZSlice, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
    cSpace%cGrid%iProcN-1, 2*cSpace%iDimX, 2*cSpace%iDimY, 2*cSpace%iDimZ, &
	0_i8b, 2*cSpace%iStartX, 2*cSpace%iStartY, 2*cSpace%iStartZ, &
	0_i8b, 0_i8b, 0_i8b, &
	2.0_dp*cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
    call InitGrid(cXYZSlice, pcGrid%iProcN, pcGrid%iProcID, (/ false,false,false,false/));
    
    ! define cContrastslicesmall (one Omega value, normal XYZ axes)
    call InitSpace(cContrastSmall, iSI_INHOMCONTRAST, cModelParams.UseYSymmetry, &
    cSpace%cGrid%iProcN-1, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ, &
	0_i8b, cSpace%iStartX, cSpace%iStartY, cSpace%iStartZ, &
	0_i8b, 0_i8b, 0_i8b, &
    cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
    call InitGrid(cContrastSmall, pcGrid%iProcN, pcGrid%iProcID, (/ false,false,false,false/));
    
    ! define cContrastslice (one Omega value, oversampled XYZ axes)
    call InitSpace(cContrastslice, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
    cSpace%cGrid%iProcN-1, 2*cSpace%iDimX, 2*cSpace%iDimY, 2*cSpace%iDimZ, &
	0_i8b, 2*cSpace%iStartX, 2*cSpace%iStartY, 2*cSpace%iStartZ, &
	0_i8b, 0_i8b, 0_i8b, &
	2.0_dp*cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
    call InitGrid(cContrastslice, pcGrid%iProcN, pcGrid%iProcID, (/ false,false,false,false/));
    
    ! read Contrast from disk; stored Contrast is already mapped and transformed
    call GridDistr2CreateEmpty(cContrastSmall%cGrid);
    call LoadContrast(cContrastSmall,"Contrast")
    call GridDistr2CreateEmpty(cContrastslice%cGrid);
    
    ! Copy Contrast into large Slice and ifft to the oversampled XYZ axes
    !	print *, " Punto di Interesse Odierno 0"
    call Distr2ObtainXYZBlock(cContrastSmall%cGrid, cContrastslice%cGrid, 0_i8b)
    !	print *, " Punto di Interesse Odierno 1"
    call ConvertSmalltoAliasingBlock(cContrastSlice)
    call TransformXYZInv(cContrastslice%cGrid,.true.)
    call DestructSpace(cContrastSmall)
    
    !FFT normalization factor for the multiplication with dt2p: for each FORWARD transform (i.e. not the BACKWARD transform..) 
    ! 1/(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ)
    ! times 1/(2*cSpace%iDimX*2*cSpace%iDimY*2*cSpace%iDimZ)
    dFFTFactor = 1.0_dp/(8.0_dp*real(cSpace%iDimX*cSpace%iDimY*cSpace%iDimZ)**2)
    cContrastSlice%cGrid%pacD2 = cContrastSlice%cGrid%pacD2 * dFFTFactor
    
    !now for each omega, multiply the (double derivative of the) field with the contrast in
    !a large, oversampled slice
    call GridDistr2CreateEmpty(cXYZslice%cGrid)
    !do only half of the omega's, since we have used SpreadFrequencies
    ! now it is only a little bit too much..
    do iOmega = 0,(pcGrid%iD2locN+1)/2-1
       
       write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iOmega, (pcGrid%iD2locN+1)/2-1
       call PrintToLog(acTemp,3)
       
       !	get XYZ block into slice
       !		print *, " Punto di Interesse Odierno 2 New"
       call Distr2ObtainXYZBlock(cSpace%cGrid, cXYZslice%cGrid, iOmega)
       !        print *, " Punto di Interesse Odierno 3 New"
       !if(pcGrid%aiD2Loc(1+iOmega)==115) then
       !
       !call ExportSlice(   trim(sOutputDir) // trim(sOutputRoot) // "CTRSlice1",&
       !					"cSlice", &
       !					cXYZslice, &
       !					(/ cXYZslice%cGrid%iProcID, 0_i8b, 0_i8b, 0_i8b /), &
       !					(/ 1_i8b, cXYZslice%iDimX, cXYZslice%iDimY, cXYZslice%iDimZ /), &
       !					.false.);
       !end if
       
       !apply tapering before OverSampling
       if (cModelParams.UseSupportTapering .EQV. .true.) then
          !			print *, " Punto di Interesse Odierno 7 New"
          call TaperXYZslice(cXYZSlice,1.0_dp)
          !			print *, " Punto di Interesse Odierno 8 New"
       end if
       !        print *, " Punto di Interesse Odierno 9 New"
       !	map XYZ - small map
       call mapVx_Small(cXYZslice)
       !		print *, " Punto di Interesse Odierno 10 New"
       call mapVy_Small(cXYZslice)
       !        print *, " Punto di Interesse Odierno 11 New"
       !	transform small in XYZ, and inverse transform large with anti-aliasing regions
       call TransformXYZ_sml(cXYZslice%cGrid)
       call ConvertSmalltoAliasingBlock(cXYZslice)
       call TransformXYZInv(cXYZslice%cGrid,.true.)
       
       !	multiply with contrast function, including renormalization due to transforms
       cXYZslice%cGrid%pacD2=cXYZslice%cGrid%pacD2 * cContrastslice%cGrid%pacD2
       
       !	transform large and inverse transform small in XYZ
       !		print *, "Punto di Interesse Odierno"
       call TransformXYZ(cXYZslice%cGrid,.true.)
       call ConvertAliasingBlocktoSmall(cXYZslice)
       call TransformXYZInv_sml(cXYZslice%cGrid)
       
       !	inverse map XYZ - small map
       call mapVxInv_Small(cXYZslice)
       call mapVyInv_Small(cXYZslice)
       
       !	put XYZ block into field
       call Distr2PutXYZBlock(cXYZslice%cGrid, cSpace%cGrid, iOmega)
       
    end do
    
    !call ExportSlice(        trim(sOutputDir) // trim(sOutputRoot) // "CTR5",&
    !                        "PLin", &
    !                        cSpace, &
    !						(/ 0_i8b, 0_i8b, -cSpace.iStartY, 0_i8b /), &
    !						(/ cSpace.iDimT+1, cSpace.iDimX, 1_i8b, cSpace.iDimZ /), &
    !						.false.);
    
    ! kill XYZslice and Contrastslice
    call DestructSpace(cXYZslice)
    call DestructSpace(cContrastslice)
    
    ! redistribute to Distr 1
    call ReorderDistr2toDistr1(pcGrid)
    
    !perform the inverse of SpreadFrequencies
    call CollectFrequencies(pcGrid)
    
    !Transform back to T-domain with small transform, 
    ! and reorder to dist. 0
    call TransformTInv_sml(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    deallocate(dTaperSupportWindow,dTaperMaxFreqWindow,dMultFactor)
    
  END SUBROUTINE InhomContrastOperator_Ali
  
  SUBROUTINE InitContrastSpace(cInhomContrast, iContrastID)
    
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
    !   The subroutine InitContrastSpace initializes the inhomogeneity contrast in
    !   function the wave speed. The type of inhomogeneity is determined by the 
    !   iContrastID identifier, which may be e.g. a sphere, a blob or a more 
    !   complex contrast function. The method of storage is determined by the 
    !   modelparams.UseAntiAliasing parameter.
    
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cInhomContrast   io   type(space)  Inhomogeneity contrast space
    !   iContrastID      i        i8b      Contrast Identifier for the choice of 
    !                                      contrast function.
    !
    type(Space), intent(inout), target ::	cInhomContrast
    integer(i8b), intent(in) ::				iContrastID
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   cContrastLarge  type(Space) Inhomogeneity contrast space (oversampled)
    !   cContrastSmall  type(Space) Inhomogeneity contrast space
    !   pcGridLarge     type(Grid)  pointer to grid (oversampled)
    !   pcGridSmall     type(Grid)  pointer to grid
    !   iIndex              i8b     index to position in data array 
    !   iIndexLarge         i8b     index to position in data array (overs.)
    !   iIndexX             i8b     beam index in the x-dimension
    !   iIndexY             i8b     beam index in the y-dimension
    !   iIndexZ             i8b     beam index in the z-dimension
    !   dFFTFactor          dp      correction factor for the FFT's
    !   iOversampling       i8b     oversampling factor
    !   dOversampling       dp      oversampling factor as dp
    !   iErr                i8b     error number
    !   acTemp              char    temporary char array for output log messages
    !
    type(Space), target ::		cContrastLarge, cContrastSmall
    type(Grid), pointer :: pcGridLarge, pcGridSmall
    integer(i8b)::		iIndex, iIndexLarge, iIndexX, iIndexY, iIndexZ
    real(dp) ::			dFFTFactor 
    real(dp), parameter ::		dOversampling = 2.0_dp
    integer(i8b), parameter ::	iOversampling = dOversampling
    integer(i8b):: 	iErr;
    character(len=1024)::       acTemp
    logical(lgt),parameter 	:: false = .false.
    
    ! *****************************************************************************
    !
    !   I/O
    !
    !   log file entries, and if necessary store the inhomogeneity contrast to
    !   temporary file
    !   
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   PrintToLog
    !   InitSpace
    !   InitGrid
    !   GridDistr2CreateEmpty
    !   dGaussianblobcontrast
    !   dSpherecontrast
    !   dLunebergLensContrast
    !   dGaussianblobcontrast
    !   TaperXYZSlice
    !   MapVx
    !   MapVy
    !   FFilterSpaceSpatial3D
    !   MapVxInv   
    !   MapVyInv
    !   TransformXYZ
    !   storecontrast
    !   DestructSpace
    !
    ! =============================================================================
    
    call PrintToLog("InitContrastSpace",3)
    
    !Initialize a space with oversampling - time dimension is needed for each processor to have one slice in Dist 2...
    call InitSpace(cContrastLarge, cInhomContrast%iSpaceIdentifier, cInhomContrast%bYSymm, &
    cInhomContrast.cGrid.iProcN-1, iOversampling*cInhomContrast.iDimX, iOversampling*cInhomContrast.iDimY, iOversampling*cInhomContrast.iDimZ,  &
	0_i8b, iOversampling*cInhomContrast.iStartX, iOversampling*cInhomContrast.iStartY, iOversampling*cInhomContrast.iStartZ,  &
	0_i8b, 0_i8b, 0_i8b, &
    dOversampling*cInhomContrast%dFnyq, cInhomContrast.dTanX, cInhomContrast.dTanY, cInhomContrast.dTanT )
    call InitGrid(cContrastLarge, cInhomContrast.cGrid.iProcN, cInhomContrast.cGrid.iProcID, (/ false,false,false,false/));                                ! Define the structure to store the data on the near-field
    
    pcGridLarge => cContrastLarge%cGrid
    
    !Create grid in Distr 2
    call GridDistr2CreateEmpty(pcGridLarge)
    
    ! Fill the grid with a simple form or a combination of simple forms
    ! NB The contrast space contains the contrast relative to c0: (c-c0)/c0
    ! In the contrast source term, a different term is needed, namely (1 - (c0/c)**2)
    ! Conversion occurs in the LinContrastOperator routines
    call PrintToLog("Fill the large contrast space",4)
    if (iContrastID == iCI_COMPLEXCONTRAST) then
       
       do iIndexX = 0, pcGridLarge%iD2XL-1
          do iIndexY = 0, pcGridLarge%iD2YL-1
             do iIndexZ = 0, pcGridLarge%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS
                
                !semi-realistic test
                pcGridLarge%pacD2(iIndex) = dGaussianblobcontrast(cContrastLarge, (/ iIndexX, iIndexY, iIndexZ /), &
		1.10_dp, (/ 3.0_dp, 0.0_dp, 40.0_dp /), (/ 2.5_dp, 5.0_dp, 2.5_dp /)) + &
                dFattylayerContrast(cContrastLarge,(/ iIndexX, iIndexY, iIndexZ /), &
		0.95_dp, (/ -15.0_dp, -7.5_dp, 6.5_dp /), (/ 15.0_dp, 7.5_dp, 7.5_dp /), &
		0.5_dp, 1.75_dp, 0.0_dp, 0.0_dp) + &
                dVeinContrast(cContrastLarge,(/ iIndexX, iIndexY, iIndexZ /), &
		0.95_dp, 1.00_dp, 0_i8b, (/ -10.0_dp, 0.0_dp, 20.0_dp /), 4.0_dp, 5.0_dp)
                
             end do
          end do
       end do
       
    elseif (iContrastID == iCI_SPHERE) then
       
       do iIndexX = 0, pcGridLarge%iD2XL-1
          do iIndexY = 0, pcGridLarge%iD2YL-1
             do iIndexZ = 0, pcGridLarge%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS
                
                ! Sphere with radius R=1 mm, centered at z=10 mm and contrast 1.01*c0
                pcGridLarge%pacD2(iIndex) = dSpherecontrast(cContrastLarge, &
                (/ iIndexX, iIndexY, iIndexZ /), 1.01_dp,  &
                (/ 0.0_dp, 0.0_dp, 10.0_dp /), (/ 1.0_dp, 1.0_dp, 1.0_dp /))
                
             end do
          end do
       end do
       
    elseif (iContrastID == iCI_LUNEBERG) then
       
       do iIndexX = 0, pcGridLarge%iD2XL-1
          do iIndexY = 0, pcGridLarge%iD2YL-1
             do iIndexZ = 0, pcGridLarge%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS
                
                ! Luneberg lens
                pcGridLarge%pacD2(iIndex) = dLunebergLensContrast(cContrastLarge, (/ iIndexX, iIndexY, iIndexZ /), &
                (/ 0.0_dp, 0.0_dp, 10.0_dp /), 3.0_dp)
                
             end do
          end do
       end do
       
    elseif (iContrastID == iCI_BLOB) then
       
       do iIndexX = 0, pcGridLarge%iD2XL-1
          do iIndexY = 0, pcGridLarge%iD2YL-1
             do iIndexZ = 0, pcGridLarge%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS
                
                ! Gaussian blob
                pcGridLarge%pacD2(iIndex) = dGaussianblobcontrast(cContrastLarge, (/ iIndexX, iIndexY, iIndexZ /), &
		1.1_dp, (/ 0.0_dp, 0.0_dp, 10.0_dp /), (/ 2.0_dp, 2.0_dp, 2.0_dp /))
                
             end do
          end do
       end do
    elseif (iContrastID == iCI_BUBBLE) then
       do iIndexX = 0, pcGridLarge%iD2XL-1
          do iIndexY = 0, pcGridLarge%iD2YL-1
             do iIndexZ = 0, pcGridLarge%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS
                
                ! Gaussian blob
                pcGridLarge%pacD2(iIndex) = dpointScatterercontrast(cContrastLarge, (/ iIndexX, iIndexY, iIndexZ /), &
		1.1_dp, (/ 0.0_dp, 0.0_dp, 10.0_dp /))
                
             end do
          end do
       end do
    end if
    
    !Translate the contrast -- now defined as (c-c0)/c0 -- to the 
    ! term in the contrast source -- which is (1-(c0/c)**2)
    ! This has to be done before the filtering takes place
    pcGridLarge%pacD2 = 1.0_dp - (1.0_dp/(1.0_dp+real(pcGridLarge%pacD2,dp)))**2
    
    ! Taper the grid at the start and end in X/Y/Z
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       call TaperXYZSlice(cContrastLarge,1.0_dp)
    end if
    
    ! Map and filter the grid in all dimensions, either by an ideal filter or by a half-band filter
    call MapVx(cContrastLarge)
    call MapVy(cContrastLarge)
    
    call PrintToLog("Filter the large contrast space",4)
    call FFilterSpaceSpatial3D(cContrastLarge,iOversampling)
    
    call MapVyInv(cContrastLarge)
    call MapVxInv(cContrastLarge)
    
    !For the anti-aliasing modus, we need the contrast to be stored on disk in Distr 2.
    if (cModelParams%UseAntiAliasing) then
       call PrintToLog("Reduce large to small contrast space",4)
       call InitSpace(cContrastSmall, cInhomContrast%iSpaceIdentifier, cInhomContrast%bYSymm, &
       cInhomContrast.cGrid.iProcN-1, cInhomContrast.iDimX, cInhomContrast.iDimY, cInhomContrast.iDimZ,  &
	0_i8b, cInhomContrast.iStartX, cInhomContrast.iStartY, cInhomContrast.iStartZ,  &
	0_i8b, 0_i8b, 0_i8b, &
       cInhomContrast%dFnyq, cInhomContrast.dTanX, cInhomContrast.dTanY, cInhomContrast.dTanT )
       call InitGrid(cContrastSmall, cInhomContrast.cGrid.iProcN, cInhomContrast.cGrid.iProcID, (/ false,false,false,false/));                                ! Define the structure to store the data on the near-field
       
       pcGridSmall => cContrastSmall%cGrid
       
       !Create grid in Distr 2
       call GridDistr2CreateEmpty(pcGridSmall)
       !copy the appropriate values from the large to the small grid in Distr 2
       do iIndexX=0,pcGridSmall%iD2XL-1
          do iIndexY=0,pcGridSmall%iD2YL-1
             do iIndexZ=0,pcGridSmall%iD2ZL-1
                
                iIndex = 1 + iIndexX*pcGridSmall%iD2XS + iIndexY*pcGridSmall%iD2YS + iIndexZ*pcGridSmall%iD2ZS
                iIndexLarge = 1 + &
                ((iIndexX+cContrastSmall%iStartX)*iOversampling-cContrastLarge%iStartX)*pcGridLarge.iD2XS + &
                ((iIndexY+cContrastSmall%iStartY)*iOversampling-cContrastLarge%iStartY)*pcGridLarge.iD2YS + &
                ((iIndexZ+cContrastSmall%iStartZ)*iOversampling-cContrastLarge%iStartZ)*pcGridLarge.iD2ZS
                
                pcGridSmall%pacD2(iIndex) = &
                pcGridLarge%pacD2(iIndexLarge)
                
             end do
          end do
       end do
       
       call mapVx(cContrastSmall)
       call mapVy(cContrastSmall)
       !		print *, "Punto di Interesse Odierno 2"
       call TransformXYZ(cContrastSmall%cGrid,.true.)
       
       !FFT normalization factor for the contrast: for each FORWARD transform (i.e. not the BACKWARD transform..) 
       ! 1/(cInhomContrast%iDimX*cInhomContrast%iDimY*cInhomContrast%iDimZ)
       dFFTFactor = 1.0_dp/(real(cInhomContrast%iDimX*cInhomContrast%iDimY*cInhomContrast%iDimZ))
       
       ! Apply FFT normalization factor to the inhomogeneity contrast
       cContrastSmall%cGrid%pacD2 = &
       real(cContrastSmall%cGrid%pacD2,dp) * dFFTFactor 
       
       call PrintToLog("Store contrast space on disk",4)
       
       !store small grid
       call storecontrast(cContrastSmall,"Contrast")
       
       call DestructSpace(cContrastSmall)
    end if
    
    ! Create a grid in the cInhomContrast space and copy the appropriate values into it
    ! be sure that the spacing is correct, that -Istart is taken correctly (i.e. at zero, not at 1/2 or so)
    ! This is distr. 0; for the anti-aliasing procedure, this step isn't really necessary, but it's nice to
    ! be able to store the contrastspace in InitInhomContrast().
    call PrintToLog("Fill cInhomContrast",4)
    pcGridSmall => cInhomContrast%cGrid
    call GridDistr0CreateEmpty(pcGridSmall)
    do iIndex = 0, pcGridSmall.iD0LocN-1
       
       !find the correct index in the Large grid
       iIndexX = (pcGridSmall%aiD0Loc(1+iIndex,1)+cInhomContrast%iStartX)*iOversampling-cContrastLarge%iStartX
       iIndexY = (pcGridSmall%aiD0Loc(1+iIndex,2)+cInhomContrast%iStartY)*iOversampling-cContrastLarge%iStartY
       iIndexZ = (pcGridSmall%aiD0Loc(1+iIndex,3)+cInhomContrast%iStartZ)*iOversampling-cContrastLarge%iStartZ
       
       pcGridSmall%parD0(1 + iIndex*pcGridSmall%iD0IS) = &
       real(pcGridLarge%pacD2(1 + iIndexX*pcGridLarge%iD2XS + iIndexY*pcGridLarge%iD2YS + iIndexZ*pcGridLarge%iD2ZS),dp)		
       
    end do
    
    ! Destroy the large grid
    call DestructSpace(cContrastLarge)
    
  END SUBROUTINE InitContrastSpace
  
  FUNCTION dpointScatterercontrast(cSpace, iBeamIndex, dContrast, dLocation)
    
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
    !   The function dpointScatterercontrast returns a point in space of a 
    !   (filtered) point scatterer.
    
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dContrast   i       dp       Magnitude of the contrast in c/c0
    !   dLocation   i       dp       Location of the point scatterer in the global
    !                                grid in mm.
    !   dpointScatterercontrast   o   dp   value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) ::	cSpace
    integer(i8b), intent(in) ::	iBeamIndex(3)	
    real(dp), intent(in) ::	dContrast			
    real(dp), intent(in) ::	dLocation(3)	    
    real(dp) :: dpointScatterercontrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dLocationN   dp   Global coordinate of the location, normalized with 
    !                     respect to the wavelength (as employed everywhere in the 
    !                     code)
    !   dRadius      dp   Distance of the requested point in space to the location
    !                     of the point scatterer
    !   dGlobx       dp   x-coordinate of the requested point
    !   dGloby       dp   y-coordinate of the requested point
    !   dGlobz       dp   z-coordinate of the requested point
    !   dK_c         dp   spatial angular cutoff frequency of the filter
    !
    real(dp) ::			dLocationN(3),dRadius
    real(dp) ::			dGlobX,dGlobY,dGlobZ
    real(dp) ::			dK_c
    
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
    
    !normalized and derived params
    dLocationN = dLocation * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dK_c = pi / cSpace%dDx
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1)  + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2)  + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3)  + cSpace.iStartZ ) * cSpace.dDx
    dRadius = sqrt((dGlobX-dLocationN(1))**2 + (dGlobY-dLocationN(2))**2 + (dGlobZ-dLocationN(3))**2)
    
    !Filtered point scatterer, from vd Veeken...
    if (dRadius<1.0E-10) then
       dpointScatterercontrast = (dContrast) * dK_c**3 / (6.0_dp * pi**2)
    else
       dpointScatterercontrast = (dContrast ) * 1.0_dp / (2.0_dp * pi**2)  * &
       ( sin(dK_c * dRadius) / dRadius**3 - dK_c * cos(dK_c * dRadius) / dRadius**2)
    end if
    
  END FUNCTION dpointScatterercontrast
  
  FUNCTION dGaussianblobcontrast(cSpace, iBeamIndex, dContrast, dBlobLocation, dFWHM)
    
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
    !   The function dGaussianblobcontrast returns a point in space of a 
    !   Gaussian blob contrast - the size of the blob may vary in each of the three
    !   dimensions.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dContrast   i       dp       Magnitude of the contrast, max(c/c0)
    !   dBlobLocation   i   dp       Location of the blob in the global grid in mm.
    !   dFWHM       i       dp       Full Width at Half Maximum of the blob in each
    !                                dimension, in mm.
    !   dGaussianblobcontrast   o   dp   value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) ::	cSpace
    integer(i8b), intent(in) ::	iBeamIndex(3)	! Location of point, given as a beam index 
    real(dp), intent(in) ::	dContrast			! max(abs(c)) = dContrast * c_0
    real(dp), intent(in) ::	dBlobLocation(3)	! Location of the blob (mm)
    real(dp), intent(in) ::	dFWHM(3)			! Full Width at Half Maximum (mm)
    real(dp) :: dGaussianblobcontrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dBlobLocationN   dp   Global coordinate of the location of the blob, 
    !                     normalized with respect to the wavelength (as employed 
    !                     everywhere in the code)
    !   dSigma       dp   standard deviation used to determine the size of the 
    !                     Gaussian beam in each dimension, in wavelengths
    !   dRadius2     dp   Square of the distance of the requested point in space 
    !                     to the location of the point scatterer, normalized with
    !                     sigma
    !   dGlobx       dp   x-coordinate of the requested point
    !   dGloby       dp   y-coordinate of the requested point
    !   dGlobz       dp   z-coordinate of the requested point
    !
    real(dp) ::			dBlobLocationN(3),dSigma(3),dRadius2
    real(dp) ::			dGlobX,dGlobY,dGlobZ
    
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
    
    !normalized and derived params
    dBlobLocationN = dBlobLocation * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dSigma = dFWHM * cModelParams.freq0/(cMediumParams.c0 * 1.0d3) / sqrt(2.0 * log(2.0))
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1) + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2) + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3) + cSpace.iStartZ) * cSpace.dDx
    dRadius2 = ((dGlobX-dBlobLocationN(1))/dSigma(1))**2 + &
    ((dGlobY-dBlobLocationN(2))/dSigma(2))**2 + &
    ((dGlobZ-dBlobLocationN(3))/dSigma(3))**2
    
    dGaussianblobcontrast = (dContrast - 1.0_dp) * exp(-(0.5_dp*dRadius2))
    
  END FUNCTION dGaussianblobcontrast
  
  FUNCTION dSpherecontrast(cSpace, iBeamIndex, dContrast, dSphereLocation, dSphereRadius)
    
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
    !   The function dSpherecontrast returns a point in space of a sphere contrast.
    !   The radius may vary in each dimension (so: it is not necessarily a sphere).
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dContrast   i       dp       Magnitude of the contrast, c/c0
    !   dSphereLocation  i  dp       Location of the sphere in the global grid in 
    !                                mm.
    !   dSphereRadius    i  dp       Radius of the sphere in each dimension, in mm.
    !   dSpherecontrast   o   dp     value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) ::	cSpace
    integer(i8b), intent(in) ::	iBeamIndex(3)	! Location of point, given as a beam index 
    real(dp), intent(in) ::	dContrast			! max(abs(c-c0)) = dContrast * c_0
    real(dp), intent(in) ::	dSphereLocation(3)	! Location of the sphere (mm)
    real(dp), intent(in) ::	dSphereRadius(3)	! Sphere radius (mm)
    real(dp) :: dSpherecontrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dSphereLocationN   dp   Global coordinate of the location of the sphere, 
    !                     normalized with respect to the wavelength (as employed 
    !                     everywhere in the code)
    !   dSphereRadiusN     dp   Distance of the requested point in space to the 
    !                     center of the sphere, normalized wrt the wavelength.
    !   dRadius2     dp   square of the radius, each dimension normalized with
    !                     respect to the sphere radius
    !   dGlobx       dp   x-coordinate of the requested point
    !   dGloby       dp   y-coordinate of the requested point
    !   dGlobz       dp   z-coordinate of the requested point
    !
    real(dp) ::			dSphereLocationN(3),dSphereRadiusN(3),dRadius2
    real(dp) ::			dGlobX,dGlobY,dGlobZ
    
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
    
    !normalized params
    dSphereLocationN = dSphereLocation * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dSphereRadiusN = dSphereRadius * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1) + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2) + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3) + cSpace.iStartZ) * cSpace.dDx
    dRadius2 =  ((dGlobX-dSphereLocationN(1))/dSphereRadiusN(1))**2 + &
    ((dGlobY-dSphereLocationN(2))/dSphereRadiusN(2))**2 + &
    ((dGlobZ-dSphereLocationN(3))/dSphereRadiusN(3))**2
    
    if (dRadius2<=1.0_dp) then
       dSpherecontrast = (dContrast - 1.0_dp)
    else
       dSpherecontrast = 0.0_dp
    end if
    
  END FUNCTION dSpherecontrast
  
  FUNCTION dLunebergLensContrast(cSpace,iBeamIndex, dLensLocation, dLensRadius)
    
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
    !   The function dLunebergLensContrast returns a point in space of a contrast
    !   function that emulates a Luneberg lens.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dLensLocation  i  dp         Location of the lens in the global grid in 
    !                                mm.
    !   dLensRadius    i  dp         Radius of the lens in each dimension, in mm.
    !   dLunebergLensContrast   o   dp   value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) :: cSpace
    integer(i8b), intent(in) :: iBeamIndex(3)	! Location of point, given as a beam index 
    real(dp), intent(in) :: dLensLocation(3)! Center location of the lens (mm)
    real(dp), intent(in) :: dLensRadius		! Lens Radius (mm)
    real(dp) :: dLunebergLensContrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dLensLocationN   dp   Global coordinate of the location of the lens, 
    !                     normalized with respect to the wavelength (as employed 
    !                     everywhere in the code)
    !   dRadiusN     dp   Distance of the requested point in space to the center of
    !                     the lens, normalized wrt the wavelength.
    !   dRadius2     dp   square of the radius normalized with respect to the lens
    !                     radius
    !   dGlobx       dp   x-coordinate of the requested point
    !   dGloby       dp   y-coordinate of the requested point
    !   dGlobz       dp   z-coordinate of the requested point
    !   dRefractIndex  dp refraction index c0/c within the lens
    !
    real(dp) ::			dLensLocationN(3),dLensRadiusN,dRadius2
    real(dp) ::			dGlobX,dGlobY,dGlobZ
    real(dp) ::			dRefractIndex
    
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
    
    !normalized and derived params
    dLensLocationN = dLensLocation * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dLensRadiusN = dLensRadius * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1) + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2) + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3) + cSpace.iStartZ) * cSpace.dDx
    dRadius2 = (dGlobX-dLensLocationN(1))**2 + (dGlobY-dLensLocationN(2))**2 + (dGlobZ-dLensLocationN(3))**2
    
    if (dRadius2<=dLensRadiusN**2) then
       
       dRefractIndex = sqrt(2-dRadius2/dLensRadiusN**2)
       
       dLunebergLensContrast = 1.0_dp/dRefractIndex -1.0_dp
       
    else
       
       dLunebergLensContrast = 0.0_dp
       
    end if
    
  END FUNCTION dLunebergLensContrast
  
  FUNCTION dFattylayerContrast(cSpace,iBeamIndex, dContrast, dStart, dEnd, dAx, dLamx, dAy, dLamy)
    
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
    !   The function dFattylayerContrast returns a point in space of a phase screen
    !   contrast in the xy-plane, at a certain position in z, with a layer 
    !   thickness with a certain average value and modulations at the larger 
    !   z-coordinate in x and y. It is used to emulate a fatty layer, hence the 
    !   name.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dContrast   i       dp       value of the contrast in c/c0
    !   dStart      i       dp       coordinate in global space in mm giving the
    !                                start of the layer in all three dimensions
    !   dEnd        i       dp       coordinate in global space in mm giving the
    !                                end of the layer in all three dimensions. In
    !                                the z-dimension it gives the average end,
    !                                thickness modulation is imposed on it
    !   dAx         i       dp       magnitude of the layer thickness modulation 
    !                                in the x-dimension.
    !   dLamx       i       dp       wavelength of the layer thickness modulation
    !                                in the x-dimension
    !   dAy         i       dp       magnitude of the layer thickness modulation 
    !                                in the y-dimension.
    !   dLamy       i       dp       wavelength of the layer thickness modulation
    !                                in the y-dimension
    !   dFattylayerContrast   o  dp  value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) :: cSpace
    integer(i8b), intent(in) :: iBeamIndex(3)	! Location of point, given as a beam index
    real(dp), intent(in) :: dContrast			! max(abs(c/c0)) = dContrast
    real(dp), intent(in) :: dStart(3)			! start of the layer (mm)
    real(dp), intent(in) :: dEnd(3)				! end of the layer (mm)
    real(dp), intent(in) :: dAx					! Magnitude of variation in dEndZ in x-direction (mm)
    real(dp), intent(in) :: dLamx				! Wavelength of variation in x (lambda0)
    real(dp), intent(in) :: dAy					! Magnitude of variation in dEndZ in y-direction (mm)
    real(dp), intent(in) :: dLamy				! Wavelength of variation in y (lambda0)
    real(dp) :: dFattylayerContrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dStartN    dp   Global coordinate of the starting location of the layer,
    !                   normalized with respect to the wavelength (as employed 
    !                   everywhere in the code)
    !   dEndN      dp   Global coordinate of the starting location of the layer,
    !                   normalized with respect to the wavelength
    !   dAxN       dp   Magnitude of the layer thickness modulation in x, 
    !                   normalized   
    !   dAyN       dp   Magnitude of the layer thickness modulation in y, 
    !                   normalized
    !   dGlobx     dp   x-coordinate of the requested point
    !   dGloby     dp   y-coordinate of the requested point
    !   dGlobz     dp   z-coordinate of the requested point
    !   dBoundary  dp thickness of the layer at this x and y
    !
    real(dp) ::			dStartN(3),dEndN(3),dAxN,dAyN
    real(dp) ::			dGlobX,dGlobY,dGlobZ,dBoundary
    
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
    
    !normalized and derived params
    dStartN = dStart * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dEndN = dEnd * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dAxN = dAx * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dAyN = dAy * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1) + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2) + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3) + cSpace.iStartZ) * cSpace.dDx
    
    if(dLamX/= 0) then
       if (dLamY/=0) then
          dBoundary = dEndN(3) + dAxN * sin(two_pi*dGlobX/dLamx) + dAyN * sin(two_pi*dGlobY/dLamY)
       else
          dBoundary = dEndN(3) + dAxN * sin(two_pi*dGlobX/dLamx)
       end if
    else
       dBoundary = dEndN(3) + dAyN * sin(two_pi*dGlobY/dLamY)
    end if
    
    if ((dGlobX >= dStartN(1) .and. dGlobX <= dEndN(1)) .and. &
       (dGlobY >= dStartN(2) .and. dGlobY <= dEndN(2)) .and. &
       (dGlobZ >= dStartN(3) .and. dGlobZ <= dBoundary)) then
       
       dFattylayerContrast = dContrast - 1.0_dp
       
    else
       
       dFattylayerContrast = 0.0_dp
       
    end if
    
  END FUNCTION dFattylayerContrast
  
  FUNCTION dVeinContrast(cSpace,iBeamIndex, dWallContrast, dBloodContrast, iVeinOrientation, &
    dVeinLocation, dVeinInnerRadius, dVeinOuterRadius)
    
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
    !   The function dFattylayerContrast returns a point in space of a phase screen
    !   contrast in the xy-plane, at a certain position in z, with a layer 
    !   thickness with a certain average value and modulations at the larger 
    !   z-coordinate in x and y. It is used to emulate a fatty layer, hence the 
    !   name.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace      i   type(space)  The space in which the point has to be given
    !   iBeamIndex  i       i8b      Beam index of the point which has to be 
    !                                returned
    !   dWallContrast   i   dp       value of the contrast of the vein wall in c/c0
    !   dBloodContrast  i   dp       value of the contrast of the lumen in c/c0
    !   iVeinorientation   i   dp    orientation of the vein - 0 -> orientation is 
    !                                in y dimension; >0 -> orientation is in x 
    !                                dimension
    !                                start of the layer in all three dimensions
    !   dVeinLocation   i   dp       Location of the center of the vein in global
    !                                coordinates, in mm
    !   dVeinInnerRadius   i   dp    radius of the lumen, in mm
    !   dVeinOuterRadius   i   dp    outer radius of the vein wall, in mm
    !   dVeinContrast      o  dp     value of the contrast in (c-c0)/c0 at the 
    !                                requested point in space
    !
    type(Space), intent(in) :: cSpace
    integer(i8b), intent(in) :: iBeamIndex(3)	! Location of point, given as a beam index 
    real(dp), intent(in) :: dWallContrast		! Contrast of vein wall
    real(dp), intent(in) :: dBloodContrast		! Contrast of lumen (blood)
    integer(i8b), intent(in) :: iVeinOrientation ! 0 -> orientation is in y dimension; >0 -> orientation is in x dimension
    real(dp), intent(in) :: dVeinLocation(3)! Center location of the lens (mm)
    real(dp), intent(in) :: dVeinInnerRadius	! Vein inner wall Radius (mm)
    real(dp), intent(in) :: dVeinOuterRadius	! Vein outer wall Radius (mm)
    real(dp) :: dVeinContrast
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dVeinLocationN(3)    dp   Global coordinate of the center location of the 
    !                   vein, normalized with respect to the wavelength (as 
    !                   employed everywhere in the code)
    !   dVeinInnerRadiusN   dp    radius of the lumen, in wavelengths
    !   dVeinOuterRadiusN   dp    outer radius of the vein wall, in wavelengths
    !   dRadius2   dp   square of the perpendicular distance of the point in space 
    !                   to the center axis of the vein
    !   dGlobx     dp   x-coordinate of the requested point
    !   dGloby     dp   y-coordinate of the requested point
    !   dGlobz     dp   z-coordinate of the requested point
    !
    real(dp) ::			dVeinLocationN(3),dVeinInnerRadiusN,dVeinOuterRadiusN,dRadius2
    real(dp) ::			dGlobX,dGlobY,dGlobZ
    
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
    
    !normalized and derived params
    dVeinLocationN = dVeinLocation * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dVeinInnerRadiusN = dVeinInnerRadius * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    dVeinOuterRadiusN = dVeinOuterRadius * cModelParams.freq0/(cMediumParams.c0 * 1.0d3)
    
    !KH determine value of rContrast for this position
    ! Determine current location in global grid
    dGlobX = (iBeamIndex(1) + cSpace.iStartX + iBeamOffsetX(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobY = (iBeamIndex(2) + cSpace.iStartY + iBeamOffsetY(cSpace,iBeamIndex(3))) * cSpace.dDx
    dGlobZ = (iBeamIndex(3) + cSpace.iStartZ) * cSpace.dDx
    if (iVeinOrientation==0) then
       dRadius2 = (dGlobX-dVeinLocationN(1))**2 + (dGlobZ-dVeinLocationN(3))**2
    else
       dRadius2 = (dGlobY-dVeinLocationN(2))**2 + (dGlobZ-dVeinLocationN(3))**2
    end if
    
    if (dRadius2<=dVeinInnerRadiusN**2) then
       
       dVeinContrast = dBloodContrast - 1.0_dp
       
    else if (dRadius2<=dVeinOuterRadiusN**2) then
       
       dVeinContrast = dWallContrast - 1.0_dp
       
    else
       
       dVeinContrast = 0.0_dp
       
    end if
    
  END FUNCTION dVeinContrast
  
  SUBROUTINE SpreadFrequencies(cGrid)
    
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
    !   The subroutine SpreadFrequencies distributes the frequencies in the first
    !   half of the frequency range in such a way over the total range that each
    !   processor gets an equal share when the Grid is redistributed from T-local
    !   (Distribution 0 and 1) to XYZ-local (Distribution 2). This subroutine is 
    !   employed in InhomContrastOperator_Ali to prevent an unbalanced workload
    !   between the processors for a space where only the first half of the
    !   frequency axis is actually used. The current implementation has a fixed 
    !   time/frequency dimension in Distr. 1 of 2*iDimT time instants to include
    !   the wraparound region necessary for the convolution operation, resulting
    !   in (iDimT+1) frequencies in the temporal Fourier domain (only the positive
    !   frequencies). However, for the current function this is not optimal, as 
    !   we have only (iDimT+1)/2 time instants/frequencies; to ensure that in 
    !   Distr. 2 the workload is evenly distributed over all processors, we spread
    !   the frequencies from the first half of the 2*iDimT range over the entire 
    !   range by splitting the first half of the frequency range in iProcN/2 
    !   blocks (iProcN being the total number of processors), and copying the 
    !   second half of each block to the equivalent block in the second half of 
    !   the frequency range.
    !   This approach at least distributes the workload evenly, but of course it is
    !   not very memory efficient. We may consider a less fixed definition of the 
    !   dimension lengths for each distribution in a future revision of the code.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cGrid   io   type(grid)   The grid which frequencies have to be spread. The
    !                             data is supposed to be in Distr. 1.
    !
    type(Grid), intent(inout) :: cGrid
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iSrcIndex    i8b   source index in the data array
    !   iDestIndex   i8b   destination index in the data array
    !   iI           i8b   index over the blocks within the first half of the
    !                      frequency range
    !   iJ           i8b   index over the frequencies within the second half of 
    !                      each block
    !   iK           i8b   index over the XYZ-positions within each frequency
    !
    integer(i8b) :: iSrcIndex, iDestIndex, iI, iJ, iK	
    
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
    
    call PrintToLog("Spread frequencies",4)
    
    !for each 0...floor(iProcN/2)-1 block of iD2LocN omega rows,
    ! take the upper half ceil(iD2LocN/2)..iD2LocN-1 omega rows,
    ! and copy them to the first half 0...ceil(iD2LocN/2)-1 of 
    ! the block of iD2LocN omega rows starting from ceil(iProcN/2)
    do iI=0,cGrid%iProcN/2-1
       do iJ=(cGrid%iD2LocN+1)/2,cGrid%iD2LocN-1
          
          iSrcIndex = 1 + (iI*cGrid%iD2LocN + iJ)*cGrid%iD1TS
          iDestIndex = 1 + (((cGrid%iProcN+1)/2+iI)*cGrid%iD2LocN + (iJ - (cGrid%iD2LocN+1)/2))*cGrid%iD1TS
          
          do iK = 0,cGrid%iD1LocN-1
             
             cGrid%pacD1(iDestIndex) = cGrid%pacD1(iSrcIndex)
             
             iSrcIndex = iSrcIndex + cGrid%iD1IS
             iDestIndex = iDestIndex + cGrid%iD1IS
             
          end do
          
       end do
    end do
    
  END SUBROUTINE SpreadFrequencies
  
  SUBROUTINE CollectFrequencies(cGrid)
    
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
    !   The subroutine CollectFrequencies gathers the frequencies over the total
    !   frequency range to the first half of the frequency range that were spread
    !   by the subroutine SpreadFrequencies, i.e. it performs the inverse operation
    !   of SpreadFrequencies.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cGrid   io   type(grid)   The grid which frequencies have to be spread. The
    !                             data is supposed to be in Distr. 1.
    !
    type(Grid), intent(inout) :: cGrid
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iSrcIndex    i8b   source index in the data array
    !   iDestIndex   i8b   destination index in the data array
    !   iI           i8b   index over the blocks within the first half of the
    !                      frequency range
    !   iJ           i8b   index over the frequencies within the second half of 
    !                      each block
    !   iK           i8b   index over the XYZ-positions within each frequency
    !
    integer(i8b) :: iSrcIndex, iDestIndex, iI, iJ, iK	
    
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
    
    call PrintToLog("Collect frequencies",4)
    
    !Do the exact inverse of the spread_frequencies operation
    ! Only the iSrcIndex and iDestIndex lines are interchanged
    do iI=0,cGrid%iProcN/2-1
       do iJ=(cGrid%iD2LocN+1)/2,cGrid%iD2LocN-1
          
          iSrcIndex = 1 + (((cGrid%iProcN+1)/2+iI)*cGrid%iD2LocN + (iJ - (cGrid%iD2LocN+1)/2))*cGrid%iD1TS
          iDestIndex = 1 + (iI*cGrid%iD2LocN + iJ)*cGrid%iD1TS
          
          do iK = 0,cGrid%iD1LocN-1
             
             cGrid%pacD1(iDestIndex) = cGrid%pacD1(iSrcIndex)
             
             iSrcIndex = iSrcIndex + cGrid%iD1IS
             iDestIndex = iDestIndex + cGrid%iD1IS
             
          end do
          
       end do
    end do
    
  END SUBROUTINE CollectFrequencies
  
  SUBROUTINE StoreContrast(cSpace,acStoredContrastName)
    
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
    !   The subroutine StoreContrast stores the contrast space to a temporary file.
    !   only the data array is stored. It is assumed to be in Distribution 2, and 
    !   is stored as such.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)   The space whose data array needs to be stored
    !   acStoredContrastName   i   char   filename of the temporary file,
    !                               excluding the path to the temp dir
    !
    type(Space), intent(inout)::            cSpace
    character(*), intent(in):: acStoredContrastName
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iErr            i8b    Error number
    !   iLi             i8b    Loop counter over the blocks to store
    !   iLast           i8b    Index of the last full block
    !   iBlockL         i8b    Size of the blocks
    !   acFilename      char   Total filename including path and suffix
    !   acTemp          char   Temporary char array for output log messages
    !
    integer(i8b)::                          iErr
    integer(i8b)::                          iLi, iLast;
    character(len=1024)::                   acFileName,acTemp;
    integer(i8b), parameter::               iBlockL                = 2**18;
    
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
    !   SWStartAndCount
    !   SWStop
    !   int2str
    !   PrintToLog
    !
    ! =============================================================================
    
    call PrintToLog("Store Contrast",4)
    
    if (cSpace%cGrid%iDistr /= 2) then
       write (acTemp, '("Error occurred during the execution of StoreContrast, aborting program")');
       call PrintToLog(acTemp, -1);
       write (acTemp, '("the grid provided as input is not in distribution 2, but in distribution ", I3, "")') cSpace%cGrid%iDistr;
       call PrintToLog(acTemp, -1);
       stop
    end if
    
    !write the contrast blockwise to prevent a large temporary array to be created on stack by Fortran
    call SwStartAndCount(cswDisk)
    call PrintToLog('Write P_lin  to disk', 2);
    acFileName        = trim(sOutputDir)//acStoredContrastName//int2str(cSpace%cGrid%iProcID);
    open (unit=iExportUNIT,file=trim(acFileName),status="REPLACE",form="UNFORMATTED",iostat=iErr)
    if(iErr/=0) then
       write(acTemp,"('Error in StoreContrast, file ',A20,'does not exist')") acStoredContrastName
       call PrintToLog(acTemp,-1)
       stop
    end if
    iLast = int(cSpace%cGrid%iD2LocSize/iBlockL);
    do iLi = 0, iLast-1
       write (unit=iExportUNIT, iostat=iErr) cSpace%cGrid%pacD2(1+iLi*iBlockL:(iLi+1)*iBlockL);
    end do
    write (unit=iExportUNIT, iostat=iErr) cSpace%cGrid%pacD2(1+iLast*iBlockL:cSpace%cGrid%iD2LocSize);
    close (unit=iExportUNIT);
    
    call SWStop(cswDisk)
    
  END SUBROUTINE StoreContrast
  
  SUBROUTINE LoadContrast(cSpace,acStoredContrastName)
    
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
    !   The subroutine LoadContrast loads the inhomogeneity contrast from a 
    !   temporary file in the temp dir to a given space. The data array is in 
    !   Distr. 2, it is assumed that all other information in the space is correct. 
    !   The storing was done in a blockwise fashion.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace             io   type(space)  Space to which the data array needs  
    !                                        to be loaded.
    !   acStoredContrastName  io   char      Filename of the temp storage file.
    !
    type(Space), intent(inout)::            cSpace
    character(*), intent(in):: acStoredContrastName
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iErr            i8b    Error number
    !   iLi             i8b    Loop counter over the blocks to load
    !   iLast           i8b    Index of the last full block
    !   iBlockL         i8b    Size of the blocks
    !   pacBuffer       dpc    Temporary buffer for a block load
    !   acFilename      char   Total filename including path and suffix
    !   acTemp          char   Temporary char array for output log messages
    !
    integer(i8b)::                          iErr
    integer(i8b)::                          iLi, iLast;
    integer(i8b), parameter::               iBlockL                = 2**18;
    complex(dpc), allocatable::             pacBuffer(:);
    character(len=1024)::                   acFileName
    character(len=1024)::                   acTemp;
    
    ! *****************************************************************************
    !
    !   I/O
    !
    !   log file entries and loading of the data array from temporary file
    !   
    ! *****************************************************************************
    !
    !   SUBROUTINES/FUNCTIONS CALLED
    !
    !   SwStartAndCount
    !   PrintToLog
    !   int2str
    !   SWStop
    !   
    ! =============================================================================
    
    call PrintToLog("Load Contrast",4)
    
    if (cSpace%cGrid%iDistr /= 2) then
       write (acTemp, '("Error occurred during the execution of LoadContrast, aborting program")');
       call PrintToLog(acTemp, -1);
       write (acTemp, '("the grid provided as input is not in distribution 2, but in distribution ", I3, "")') cSpace%cGrid%iDistr;
       call PrintToLog(acTemp, -1);
       stop
    end if
    
    !read the contrast blockwise to prevent a large temporary array to be created on stack by Fortran
    allocate(pacBuffer(iBlockL));
    call PrintToLog("Buffer allocated", 3)
    call PrintToLog('Read Contrast from disk', 4);
    call SwStartAndCount(cswDiskAcces)
    acFileName        = trim(sOutputDir)//acStoredContrastName//int2str(cSpace%cGrid%iProcID);
    open (unit=iExportUNIT,file=trim(acFileName),status="OLD",form="UNFORMATTED",iostat=iErr)
    if(iErr/=0) then
       write(acTemp,"('Error in LoadContrast, file ',A20,'does not exist')") acStoredContrastName
       call PrintToLog(acTemp,-1)
       stop
    end if
    iLast = int(cSpace%cGrid%iD2LocSize/iBlockL);
    do iLi = 0, iLast-1
       read (unit=iExportUNIT, iostat=iErr) pacBuffer;
       cSpace%cGrid%pacD2(1+iLi*iBlockL:(iLi+1)*iBlockL) = pacBuffer
    end do
    read (unit=iExportUNIT, iostat=iErr) pacBuffer(1:cSpace%cGrid%iD2LocSize-iLast*iBlockL);
    cSpace%cGrid%pacD2(1+iLast*iBlockL:cSpace%cGrid%iD2LocSize) = pacBuffer(1:cSpace%cGrid%iD2LocSize-iLast*iBlockL)
    close (unit=iExportUNIT);
    call SWStop(cswDiskAcces)
    
    deallocate(pacBuffer);
    call PrintToLog("Free buffer", 3)
    
    call SWStop(cswDisk)
    
  END SUBROUTINE LoadContrast
  
  SUBROUTINE ConvertSmalltoAliasingBlock(cXYZslice)
    
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
    !   The subroutine ConvertSmalltoAliasingBlock copies all frequency ranges in 
    !   the small block stored in 1/8 of the cXYZslice slice to the correct 
    !   regions: keep one part with all positive frequencies at its place, but copy
    !   all other seven parts with negative frequencies to the far ends of the large 
    !   block and put zeros in their original places. This is because the order of
    !   the frequencies with an FFT is first the positive, and then the negative 
    !   from largest negative to smallest negative: [0 1 2 3 4 ... ... -3 -2 -1]
    !   The anti-aliasing regions are always the largest positive and negative
    !   regions, which explains the necessity to copy the negative axes to the far
    !   ends of the ranges.
    !   We assume that, since cXYZslice includes the anti-aliasing regions which
    !   double the dimensions, the total dimensions have an even length.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cXYZslice   io   type(space)   XYZ slice whose spatial frequencies need
    !                                  to be reshuffled in all three dimensions 
    !
    type(Space), intent(inout), target :: cXYZslice
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid      type(Grid)   pointer to the Grid structure within cXYZslice
    !   iIndexX      i8b   loop counter in x
    !   iIndexY      i8b   loop counter in y
    !   iIndexZ      i8b   loop counter in z
    !   iStartX      i8b   start index of the block to move in x
    !   iEndX        i8b   end index of the block to move in x
    !   iDisplX      i8b   distance over which to displace the block in x
    !   iStartY      i8b   start index of the block to move in y
    !   iEndY        i8b   end index of the block to move in y
    !   iDisplY      i8b   distance over which to displace the block in y
    !   iStartZ      i8b   start index of the block to move in z
    !   iEndZ        i8b   end index of the block to move in z
    !   iDisplZ      i8b   distance over which to displace the block in z
    !
    type(Grid), pointer :: pcGrid
    integer(i8b) :: iIndexX,iIndexY,iIndexZ
    integer(i8b) :: iStartX,iEndX,iDisplX,iStartY,iEndY,iDisplY,iStartZ,iEndZ,iDisplZ
    
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
    !   MoveBlock
    !   
    ! =============================================================================
    
    pcGrid => cXYZslice%cGrid
    
    iStartX = (cXYZslice%iDimX/2+1)/2
    iEndX   = cXYZslice%iDimX/2-1
    iDisplX = cXYZslice%iDimX/2  
    
    iStartY = (cXYZslice%iDimY/2+1)/2
    iEndY   = cXYZslice%iDimY/2-1
    iDisplY = cXYZslice%iDimY/2  
    
    iStartZ = (cXYZslice%iDimZ/2+1)/2
    iEndZ   = cXYZslice%iDimZ/2-1
    iDisplZ = cXYZslice%iDimZ/2
    
    !If IDim -without anti-aliasing region -is even, we make the largest 
    ! negative frequency zero otherwise we don't keep real contrast when 
    ! transformed in 3D
    if (mod(cXYZslice%iDimX,4_i8b)==0) then
       do iIndexY = 0, iEndY
          do iIndexZ = 0, iEndZ
             pcGrid%pacD2(1 + iStartX*pcGrid%iD2XS + iIndexY*pcGrid%iD2YS + iIndexZ*pcGrid%iD2ZS) = 0.0_dp
          end do
       end do
    end if
    
    if (mod(cXYZslice%iDimY,4_i8b)==0) then
       do iIndexX = 0, iEndX
          do iIndexZ = 0, iEndZ
             pcGrid%pacD2(1 + iIndexX*pcGrid%iD2XS + iStartY*pcGrid%iD2YS + iIndexZ*pcGrid%iD2ZS) = 0.0_dp
          end do
       end do
    end if
    
    if (mod(cXYZslice%iDimZ,4_i8b)==0) then
       do iIndexX = 0, iEndX
          do iIndexY = 0, iEndY
             pcGrid%pacD2(1 + iIndexX*pcGrid%iD2XS + iIndexY*pcGrid%iD2YS + iStartZ*pcGrid%iD2ZS) = 0.0_dp
          end do
       end do
    end if
    
    ! layer with positive z frequencies
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    ! layer with negative z frequencies
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
  END SUBROUTINE ConvertSmalltoAliasingBlock
  
  SUBROUTINE ConvertAliasingBlocktoSmall(cXYZslice)
    
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
    !   The subroutine ConvertSmalltoAliasingBlock performs the inverse operation
    !   of ConvertSmalltoAliasingBlock.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cXYZslice   io   type(space)   XYZ slice whose spatial frequencies need
    !                                  to be reshuffled in all three dimensions 
    !
    type(Space), intent(inout), target :: cXYZslice
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   pcGrid      type(Grid)   pointer to the Grid structure within cXYZslice
    !   iStartX      i8b   start index of the block to move in x
    !   iEndX        i8b   end index of the block to move in x
    !   iDisplX      i8b   distance over which to displace the block in x
    !   iStartY      i8b   start index of the block to move in y
    !   iEndY        i8b   end index of the block to move in y
    !   iDisplY      i8b   distance over which to displace the block in y
    !   iStartZ      i8b   start index of the block to move in z
    !   iEndZ        i8b   end index of the block to move in z
    !   iDisplZ      i8b   distance over which to displace the block in z
    !
    type(Grid), pointer :: pcGrid
    integer(i8b) :: iStartX,iEndX,iDisplX,iStartY,iEndY,iDisplY,iStartZ,iEndZ,iDisplZ
    
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
    !   MoveBlock
    !   
    ! =============================================================================
    
    pcGrid => cXYZslice%cGrid
    
    iStartX = cXYZslice%iDimX/2 + (cXYZslice%iDimX/2+1)/2
    iEndX   = cXYZslice%iDimX-1
    iDisplX = -cXYZslice%iDimX/2  
    
    iStartY = cXYZslice%iDimY/2 + (cXYZslice%iDimY/2+1)/2
    iEndY   = cXYZslice%iDimY-1
    iDisplY = -cXYZslice%iDimY/2  
    
    iStartZ = cXYZslice%iDimZ/2 + (cXYZslice%iDimZ/2+1)/2
    iEndZ   = cXYZslice%iDimZ-1
    iDisplZ = -cXYZslice%iDimZ/2
    
    ! layer with positive z frequencies
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ 0_i8b  , iStartZ-1, 0_i8b   /))
    
    ! layer with negative z frequencies
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ 0_i8b  , iStartY-1, 0_i8b   /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ 0_i8b  , iStartX-1, 0_i8b   /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
    call MoveBlock(pcGrid,	(/ iStartX, iEndX    , iDisplX /), &
    (/ iStartY, iEndY    , iDisplY /), &
    (/ iStartZ, iEndZ    , iDisplZ /))
    
  END SUBROUTINE ConvertAliasingBlocktoSmall
  
  SUBROUTINE MoveBlock(cGrid,iXind,iYind,iZind)
    
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
    !   The subroutine MoveBlock copies a block in a Distr. 2 data array from one 
    !   position to another and fill the original position with zeros
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cGrid   io   type(grid)   grid with distr.2 data array in which a block of
    !                             data is repositioned
    !   iXind    i     i8b        3-element vector with start index, end index and
    !                             displacement in the x-dimension
    !   iYind    i     i8b        3-element vector with start index, end index and
    !                             displacement in the y-dimension
    !   iZind    i     i8b        3-element vector with start index, end index and
    !                             displacement in the z-dimension
    !
    type(Grid), intent(inout) :: cGrid
    integer(i8b), intent(in) :: iXind(3)
    integer(i8b), intent(in) :: iYind(3)
    integer(i8b), intent(in) :: iZind(3)
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   iIndexX      i8b   loop counter in x
    !   iIndexY      i8b   loop counter in y
    !   iIndexZ      i8b   loop counter in z
    !
    integer(i8b) :: iIndexX,iIndexY,iIndexZ
    
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
    
    do iIndexX = iXind(1),iXind(2)
       do iIndexY = iYind(1),iYind(2)
          do iIndexZ = iZind(1),iZind(2)
             
             cGrid%pacD2(1 + (iIndexX + iXind(3))*cGrid%iD2XS + &
             (iIndexY + iYind(3))*cGrid%iD2YS + &
             (iIndexZ + iZind(3))*cGrid%iD2ZS) = &
             cGrid%pacD2(1 + iIndexX*cGrid%iD2XS + iIndexY*cGrid%iD2YS + iIndexZ*cGrid%iD2ZS)
             
             cGrid%pacD2(1 + iIndexX*cGrid%iD2XS + iIndexY*cGrid%iD2YS + iIndexZ*cGrid%iD2ZS) = 0.0_dp
             
          end do
       end do
    end do
    
  END SUBROUTINE MoveBlock
  
  SUBROUTINE TaperXYZslice(cSlice,dTaperwidth)
    
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
    !   The function TaperXYZslice applies a tapering on an XYZ slice in all three
    !   spatial dimensions. The tapering width is the same on all sides.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSlice       io  type(Space)   XYZ slice which needs to be tapered
    !   dTaperWidth  i      dp         Width of the taperings with respect to 
    !                                  the wavelength
    !
    type(Space), intent(inout) :: cSlice
    real(dp), intent(in) :: dTaperWidth
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   dTaperX   dp    tapering window in the x-dimension
    !   dTaperY   dp    tapering window in the y-dimension
    !   dTaperZ   dp    tapering window in the z-dimension
    !   iIndex    i8b   index to the elements in the data array
    !   iIndexX   i8b   index in the x-dimension
    !   iIndexY   i8b   index in the y-dimension
    !   iIndexZ   i8b   index in the z-dimension
    !
    real(dp) :: dTaperX(cSlice%iDimX),dTaperY(cSlice%iDimY),dTaperZ(cSlice%iDimZ)
    integer(i8b) :: iIndex,iIndexX,iIndexY,iIndexZ
    
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
    !   dTaperingWindow
    !   
    ! =============================================================================
    
    dTaperX = dTaperingWindow(cSlice%iDimX,cSlice%dDx,dTaperwidth,dTaperwidth)
    dTaperY = dTaperingWindow(cSlice%iDimY,cSlice%dDx,dTaperwidth,dTaperwidth)
    dTaperZ = dTaperingWindow(cSlice%iDimZ,cSlice%dDx,dTaperwidth,dTaperwidth)
    
    do iIndexX = 0, cSlice%cGrid%iD2XL-1
       do iIndexY = 0, cSlice%cGrid%iD2YL-1
          do iIndexZ = 0, cSlice%cGrid%iD2ZL-1
             
             iIndex = 1 + iIndexX*cSlice%cGrid%iD2XS + iIndexY*cSlice%cGrid%iD2YS + iIndexZ*cSlice%cGrid%iD2ZS
             
             cSlice%cGrid%pacD2(iIndex) = cSlice%cGrid%pacD2(iIndex) * (dTaperX(1+iIndexX) * dTaperY(1+iIndexY) * dTaperZ(1+iIndexZ))
             
          end do
       end do
    end do
    
  END SUBROUTINE TaperXYZslice
  
  ! Added for NKC
  SUBROUTINE NonlinearKappaOperatorDXALi(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. The spatial 
    !   derivative is implemented in kx space.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1x,arBuffer2x
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1xC,arBuffer2xC,dKvector
    real(dp)   , dimension(cSpace.iDimX)        :: dummydKvector1,dummydKvector2,dKvectorReal
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDX
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do xindex=1,iDimX-1
       if (xindex.le.((iDimX)/2) +1) then
          dKvectorReal(xindex)=((xindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       else
          dKvectorReal(xindex)=((xindex-iDimX)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    call TransformT(pcGrid)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2xC=arBuffer2xC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=(arBuffer1xC)*(1.0_dp/iLen); ! Normalization is done with respect to the dimension of the forward Fourier transformation 
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDX
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    !    if (IProcID==1) then
    !        
    !        print *, "pcGrid%iD0IS"
    !        print *, pcGrid%iD0IS
    !        print *, "pcGrid%iD0TL"
    !        print *, pcGrid%iD0TL
    !        print *, "pcGrid%iD0TS"
    !        print *, pcGrid%iD0TS
    !    
    !    end if 
    
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDX)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    !here we perform the fft transform, multply with Kx and go back to time domain. Normalization has to be performed with respect to the forward fft
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    iLen=pcGrid%iD2XL
    do iIndex=0,(pcGrid%iD2LocSize/pcGrid%iD2XL)-1
       
       arBuffer1x=pcGrid%pacD2(iStart+1:iStart+pcGrid%iD2XL)
       
       ! fft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       arBuffer2xC=arBuffer2xC*dKvector    
       
       ! ifft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       
       !   		            pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*cSpace%dDx*(1.0_dp/iLen); 
       pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*(1.0_dp/iLen); 
       
       
       iStart=iStart+pcGrid%iD2XL
    end do
    
    !
    !if (iProcID==1) then
    !    print *, "iD2LocSize"
    !    print *, pcGrid%iD2LocSize
    !    print *, "iD2GlobN"
    !    print *, pcGrid%iD2GlobN
    !    print *, "iD2LocN"
    !    print *, pcGrid%iD2LocN
    !    print *, "iD2XS"
    !    print *, pcGrid%iD2XS
    !    print *, "iD2YS"
    !    print *, pcGrid%iD2YS
    !    print *, "iD2ZS"
    !    print *, pcGrid%iD2ZS
    !    print *, "iD2IS"
    !    print *, pcGrid%iD2IS
    !    print *, "iD2XL"
    !    print *, pcGrid%iD2XL
    !    print *, "iD2YL"
    !    print *, pcGrid%iD2YL
    !    print *, "iD2ZL"
    !    print *, pcGrid%iD2ZL  
    !end if 
    
    ! we we back transform and go back to distribution 0
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorDXALi
  
  SUBROUTINE NonlinearKappaOperatorDYALi(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1y,arBuffer2y
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1yC,arBuffer2yC,dKvector
    real(dp)   , dimension(cSpace.iDimY)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDY
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do yindex=1,iDimY-1
       if (yindex.le.((iDimY)/2) +1) then
          dKvectorReal(yindex)=((yindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       else
          dKvectorReal(yindex)=((yindex-iDimY)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    call TransformT(pcGrid)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2yC=arBuffer2yC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=(arBuffer1yC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    !	
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDY
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDY)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart=0
    iStart2=0
    Timestart=0
    iLen=pcGrid%iD2YL
    
    do tindex=0, pcGrid%iD2LocN-1 ! This is the loop for time instants stored locally
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex ! This is the index that takes into account for the different starting point corresponding to a ginve time index
       
       do xindex=0, pcGrid%iD2XL-1 ! This is the loop for x points
          
          iStart=0+xindex+Timestart  ! This is the index
          iStart2=0+xindex+Timestart ! This is the index 2
          
          ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
          
          do iIndex=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points
                
                arBuffer1y(iIndex2+1) = pcGrid%pacD2(iStart+1)	
                
                iStart=iStart+pcGrid%iD2XL
                
             end do
             
             ! fft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             arBuffer2yC=arBuffer2yC*dKvector ! Here the transformed vector is multiplied with the k vector   
             
             ! ifft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points to put the outcome into the pcGrid
                
                pcGrid%pacD2(iStart2+1)=arBuffer1y(iIndex2+1)*(1.0_dp/iLen) 
                
                iStart2=iStart2+pcGrid%iD2XL
                
             end do
             
             
          end do
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
  END SUBROUTINE NonlinearKappaOperatorDYALi
  
  SUBROUTINE NonlinearKappaOperatorDZALi(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. This implementation does not
    !   prevent aliasing in the spatial dimension. The spatial 
    !   derivative is implemented with a Finite Difference approximation
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1z,arBuffer2z
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1zC,arBuffer2zC,dKvector
    real(dp)   , dimension(cSpace.iDimZ)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDZ
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do zindex=1,iDimZ-1
       if (zindex.le.((iDimZ)/2) +1) then
          dKvectorReal(zindex)=((zindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       else
          dKvectorReal(zindex)=((zindex-iDimZ)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    call TransformT(pcGrid)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square**2 + im * arBuffer2square**2;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimZ
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=(arBuffer1zC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDZ
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDZ)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart2=0
    iStart=0
    Timestart=0
    iLen=pcGrid%iD2ZL
    
    do tindex=0, pcGrid%iD2LocN-1
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex
       
       do xindex=0, (pcGrid%iD2XL*pcGrid%iD2YL)-1
          
          iStart=0+xindex+Timestart
          iStart2=0+xindex+Timestart
          
          
          ! Now we calculate dz(dzk p^2) so that afterwards we can apply the derivative to this term
          
          
          !            do iIndex2=0,pcGrid%iD2ZL-1
          !        
          !			    pcGrid%pacD2(iStart+1) = pcGrid%pacD2(iStart+1)*dKvector(iIndex2+1)*(cSpace%dDx**3/( real(2*cSpace%iDimT * cSpace%iDimX,dp) * cSpace%iDimY * cSpace%iDimZ))
          !			    
          !			
          !			    iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
          !					
          !		    end do
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             arBuffer1z(iIndex2+1) = pcGrid%pacD2(iStart+1)	
             
             iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
             
          end do
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector ! Here the transformed vector is multiplied with the k vector   
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points to put the outcome into the pcGrid
             
             pcGrid%pacD2(iStart2+1)=arBuffer1z(iIndex2+1)*(1.0_dp/iLen) 
             
             iStart2=iStart2+pcGrid%iD2XL*pcGrid%iD2YL
             
             
          end do
          
          
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorDZALi
  !------------------------------------------------------
  SUBROUTINE NonlinearKappaOperatorDirectionDXAli(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1x,arBuffer2x
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1xC,arBuffer2xC,dKvector
    real(dp)   , dimension(cSpace.iDimX)        :: dummydKvector1,dummydKvector2,dKvectorReal
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDX
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do xindex=1,iDimX-1
       if (xindex.le.((iDimX)/2) +1) then
          dKvectorReal(xindex)=((xindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       else
          dKvectorReal(xindex)=((xindex-iDimX)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    
    call TransformT(pcGrid)
    
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2xC=arBuffer2xC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=(arBuffer1xC)*(1.0_dp/iLen); ! Normalization is done with respect to the dimension of the forward Fourier transformation 
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDX
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    !here we perform the fft transform, multply with Kx and go back to time domain. Normalization has to be performed with respect to the forward fft
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    iLen=pcGrid%iD2XL
    do iIndex=0,(pcGrid%iD2LocSize/pcGrid%iD2XL)-1
       
       arBuffer1x=pcGrid%pacD2(iStart+1:iStart+pcGrid%iD2XL)
       
       ! fft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       arBuffer2xC=arBuffer2xC*dKvector    
       
       ! ifft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       
       !   		            pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*cSpace%dDx*(1.0_dp/iLen); 
       pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*(1.0_dp/iLen); 
       
       
       iStart=iStart+pcGrid%iD2XL
    end do
    
    !
    !if (iProcID==1) then
    !    print *, "iD2LocSize"
    !    print *, pcGrid%iD2LocSize
    !    print *, "iD2GlobN"
    !    print *, pcGrid%iD2GlobN
    !    print *, "iD2LocN"
    !    print *, pcGrid%iD2LocN
    !    print *, "iD2XS"
    !    print *, pcGrid%iD2XS
    !    print *, "iD2YS"
    !    print *, pcGrid%iD2YS
    !    print *, "iD2ZS"
    !    print *, pcGrid%iD2ZS
    !    print *, "iD2IS"
    !    print *, pcGrid%iD2IS
    !    print *, "iD2XL"
    !    print *, pcGrid%iD2XL
    !    print *, "iD2YL"
    !    print *, pcGrid%iD2YL
    !    print *, "iD2ZL"
    !    print *, pcGrid%iD2ZL  
    !end if 
    
    ! we we back transform and go back to distribution 0
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    iStart=0
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    DEALLOCATE(KAPPAcontrastDX)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDXAli
  
  SUBROUTINE NonlinearKappaOperatorDirectionDYALi(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1y,arBuffer2y
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1yC,arBuffer2yC,dKvector
    real(dp)   , dimension(cSpace.iDimY)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDY
    type(Space)                                 :: cSliceS
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do yindex=1,iDimY-1
       if (yindex.le.((iDimY)/2) +1) then
          dKvectorReal(yindex)=((yindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       else
          dKvectorReal(yindex)=((yindex-iDimY)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    
    call TransformT(pcGrid)
    
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2yC=arBuffer2yC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=(arBuffer1yC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    !	
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDY
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart=0
    iStart2=0
    Timestart=0
    iLen=pcGrid%iD2YL
    
    do tindex=0, pcGrid%iD2LocN-1 ! This is the loop for time instants stored locally
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex ! This is the index that takes into account for the different starting point corresponding to a ginve time index
       
       do xindex=0, pcGrid%iD2XL-1 ! This is the loop for x points
          
          iStart=0+xindex+Timestart  ! This is the index
          iStart2=0+xindex+Timestart ! This is the index 2
          
          ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
          
          do iIndex=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points
                
                arBuffer1y(iIndex2+1) = pcGrid%pacD2(iStart+1)	
                
                iStart=iStart+pcGrid%iD2XL
                
             end do
             
             ! fft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             arBuffer2yC=arBuffer2yC*dKvector ! Here the transformed vector is multiplied with the k vector   
             
             ! ifft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points to put the outcome into the pcGrid
                
                pcGrid%pacD2(iStart2+1)=arBuffer1y(iIndex2+1)*(1.0_dp/iLen) 
                
                iStart2=iStart2+pcGrid%iD2XL
                
             end do
             
             
          end do
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    DEALLOCATE(KAPPAcontrastDY)
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDYALi
  
  SUBROUTINE NonlinearKappaOperatorDirectionDZALi(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr;
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1z,arBuffer2z
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1zC,arBuffer2zC,dKvector
    real(dp)   , dimension(cSpace.iDimZ)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDZ
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do zindex=1,iDimZ-1
       if (zindex.le.((iDimZ)/2) +1) then
          dKvectorReal(zindex)=((zindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       else
          dKvectorReal(zindex)=((zindex-iDimZ)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    
    
    
    call TransformT(pcGrid)
    
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    
    
    
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimZ
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=(arBuffer1zC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDZ
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart2=0
    iStart=0
    Timestart=0
    iLen=pcGrid%iD2ZL
    
    do tindex=0, pcGrid%iD2LocN-1
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex
       
       do xindex=0, (pcGrid%iD2XL*pcGrid%iD2YL)-1
          
          iStart=0+xindex+Timestart
          iStart2=0+xindex+Timestart
          
          
          ! Now we calculate dz(dzk p^2) so that afterwards we can apply the derivative to this term
          
          
          !            do iIndex2=0,pcGrid%iD2ZL-1
          !        
          !			    pcGrid%pacD2(iStart+1) = pcGrid%pacD2(iStart+1)*dKvector(iIndex2+1)*(cSpace%dDx**3/( real(2*cSpace%iDimT * cSpace%iDimX,dp) * cSpace%iDimY * cSpace%iDimZ))
          !			    
          !			
          !			    iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
          !					
          !		    end do
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             arBuffer1z(iIndex2+1) = pcGrid%pacD2(iStart+1)	
             
             iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
             
          end do
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector ! Here the transformed vector is multiplied with the k vector   
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points to put the outcome into the pcGrid
             
             pcGrid%pacD2(iStart2+1)=arBuffer1z(iIndex2+1)*(1.0_dp/iLen) 
             
             iStart2=iStart2+pcGrid%iD2XL*pcGrid%iD2YL
             
             
          end do
          
          
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDZ)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorDirectionDZALi
  !------------------------------------------------------
  
  SUBROUTINE NonlinearKappaOperatorLinDXALi(cSpace, IncSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr
    integer(i8b)                                :: iLen,iStart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1x,arBuffer2x
    complex(dp), dimension(cSpace.iDimX)        :: arBuffer1xC,arBuffer2xC,dKvector
    real(dp)   , dimension(cSpace.iDimX)        :: dummydKvector1,dummydKvector2,dKvectorReal
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDX
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do yindex=1,iDimX-1
       if (yindex.le.((iDimX)/2) +1) then
          dKvectorReal(yindex)=((yindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       else
          dKvectorReal(yindex)=((yindex-iDimY)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimX,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call ReorderDistr0ToDistr1(pcGridIn)
    
    
    call TransformT(pcGrid)
    call TransformT(pcGridIn)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGridIn%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))		
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDX(iDimX,iDimY,iDimZ))
    ! derivative with respect to X
    iLen	= iDimX;
    do yindex=1,iDimY
       do zindex=1,iDimZ
          
          
          arBuffer1x=KAPPAcontrast(1:iDimX,yindex,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2xC=arBuffer2xC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDX(1:iDimX,yindex,zindex)=(arBuffer1xC)*(1.0_dp/iLen); ! Normalization is done with respect to the dimension of the forward Fourier transformation 
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADX.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDX
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDX(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       
       
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDX)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    !here we perform the fft transform, multply with Kx and go back to time domain. Normalization has to be performed with respect to the forward fft
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    iLen=pcGrid%iD2XL
    do iIndex=0,(pcGrid%iD2LocSize/pcGrid%iD2XL)-1
       
       arBuffer1x=pcGrid%pacD2(iStart+1:iStart+pcGrid%iD2XL)
       
       ! fft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1x,arBuffer2xC,fftw_forward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       arBuffer2xC=arBuffer2xC*dKvector    
       
       ! ifft
       call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2xC,arBuffer1xC,fftw_backward,fftw_estimate+fftw_unaligned)
       call dfftw_execute(plan1d)
       call dfftw_destroy_plan(plan1d)
       
       
       !   		            pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*cSpace%dDx*(1.0_dp/iLen); 
       pcGrid%pacD2(iStart+1:iStart+cSpace%iDimX)=(arBuffer1xC)*(1.0_dp/iLen); 
       
       
       iStart=iStart+pcGrid%iD2XL
    end do
    
    !
    if (iProcID==1) then
       !    print *, "cSpace%iDimX"
       !    print *, cSpace%iDimX
       !    print *, "pcGrid%iD2XL"
       !    print *, pcGrid%iD2XL
       !    print *, "iD2LocN"
       !    print *, pcGrid%iD2LocN
       !    print *, "iD2XS"
       !    print *, pcGrid%iD2XS
       !    print *, "iD2YS"
       !    print *, pcGrid%iD2YS
       !    print *, "iD2ZS"
       !    print *, pcGrid%iD2ZS
       !    print *, "iD2IS"
       !    print *, pcGrid%iD2IS
       !    print *, "iD2XL"
       !    print *, pcGrid%iD2XL
       !    print *, "iD2YL"
       !    print *, pcGrid%iD2YL
       !    print *, "iD2ZL"
       !    print *, pcGrid%iD2ZL  
    end if
    
    ! we we back transform and go back to distribution 0
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorLinDXALi
  
  SUBROUTINE NonlinearKappaOperatorLinDYALi(cSpace,IncSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1y,arBuffer2y
    complex(dp), dimension(cSpace.iDimY)        :: arBuffer1yC,arBuffer2yC,dKvector
    real(dp)   , dimension(cSpace.iDimY)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDY
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do yindex=1,iDimY-1
       if (yindex.le.((iDimY)/2) +1) then
          dKvectorReal(yindex)=((yindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       else
          dKvectorReal(yindex)=((yindex-iDimY)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimY,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call ReorderDistr0ToDistr1(pcGridIn)
    
    
    call TransformT(pcGrid)
    call TransformT(pcGridIn)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGridIn%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))		
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDY(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimY;
    do xindex=1,iDimX
       do zindex=1,iDimZ
          
          
          arBuffer1y=KAPPAcontrast(xindex,1:iDimY,zindex)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2yC=arBuffer2yC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDY(xindex,1:iDimY,zindex)=(arBuffer1yC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    !	
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADY.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDY
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDY(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDY)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart=0
    iStart2=0
    Timestart=0
    iLen=pcGrid%iD2YL
    
    do tindex=0, pcGrid%iD2LocN-1 ! This is the loop for time instants stored locally
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex ! This is the index that takes into account for the different starting point corresponding to a ginve time index
       
       do xindex=0, pcGrid%iD2XL-1 ! This is the loop for x points
          
          iStart=0+xindex+Timestart  ! This is the index
          iStart2=0+xindex+Timestart ! This is the index 2
          
          ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
          
          do iIndex=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points
                
                arBuffer1y(iIndex2+1) = pcGrid%pacD2(iStart+1)	
                
                iStart=iStart+pcGrid%iD2XL
                
             end do
             
             ! fft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1y,arBuffer2yC,fftw_forward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             arBuffer2yC=arBuffer2yC*dKvector ! Here the transformed vector is multiplied with the k vector   
             
             ! ifft
             call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2yC,arBuffer1yC,fftw_backward,fftw_estimate+fftw_unaligned)
             call dfftw_execute(plan1d)
             call dfftw_destroy_plan(plan1d)
             
             
             do iIndex2=0,pcGrid%iD2YL-1 ! This is the loop for y points to put the outcome into the pcGrid
                
                pcGrid%pacD2(iStart2+1)=arBuffer1y(iIndex2+1)*(1.0_dp/iLen) 
                
                iStart2=iStart2+pcGrid%iD2XL
                
             end do
             
             
          end do
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
  END SUBROUTINE NonlinearKappaOperatorLinDYALi
  
  SUBROUTINE NonlinearKappaOperatorLinDZALi(cSpace,IncSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
    type(Space), target, intent(in)   ::	IncSpace
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
    INTEGER(8)                                  :: plan1D
    integer(i4b)                                :: iErr
    integer(i8b)                                :: iLen,iStart,iStart2,Timestart,iDimW, iDimT,iDimZ,iDimX,iDimY,iProcN,iProcID,iindex,iIndex2,iD0IS, xindex, yindex, zindex, tindex
    type(Grid), pointer                         :: pcGrid, pcGridIn
    type(DerivLookup)                           :: cDerivLookup
    real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT),arBufferIn1square(cSpace%iDimT),arBufferIn2square(cSpace%iDimT)
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1z,arBuffer2z
    complex(dp), dimension(cSpace.iDimZ)        :: arBuffer1zC,arBuffer2zC,dKvector
    real(dp)   , dimension(cSpace.iDimZ)        :: dummydKvector1,dummydKvector2,dKvectorReal
    !	real(dp)                                    :: arBuffer1square(cSpace%iDimT),arBuffer2square(cSpace%iDimT)
    real(dp), allocatable                       :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                    :: dLeftBand,dRightBand
    real(8), dimension(:,:,:), allocatable      :: ContrastFiltered
    complex(dp), dimension(:,:,:), allocatable  :: KAPPAContrast, KAPPAcontrastDZ
    type(Space)                                 :: cSliceS
    
    real(dp)                                    :: dOmega, dCorrection, dKcutoff
    character(len=1024)                         :: acTemp;
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("NonlinearKappaOperator",4)
    
    pcGrid=>cSpace%cGrid
    pcGridIn=>IncSpace%cGrid
    
    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    
    
    iD0IS   = pcgrid%iD0IS   !the stride of the disctibution
    iProcN  = pcgrid%iProcN  !number of processors in use
    iProcID = pcgrid%iProcID !processor identifier
    
    ! Create the K vector in Kx space
    do zindex=1,iDimZ-1
       if (zindex.le.((iDimZ)/2) +1) then
          dKvectorReal(zindex)=((zindex-1)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       else
          dKvectorReal(zindex)=((zindex-iDimZ)*two_pi*2.0_dp*cSpace%dFnyq/(real(cSpace%iDimZ,dp)*cMediumParams.c0))
       end if
    end do
    dKvector= im * dKvectorReal;
    
    !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dKvector.bin',STATUS='unknown',FORM='binary')
    !WRITE(3) dKvector
    !CLOSE(3, STATUS = 'keep')
    !dKvector=dKvectorReal*(dcmplx (0.0d0, 1.0d0))
    
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(1:cSpace%iDimT),&
    dTaperMaxFreqWindow(1:iDimW))
    
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       !use tapering of the field at start and end to prevent wraparound leakage
       
       !build tapering window, the first two periods and the last two periods are tapered
       !in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(cSpace%iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow
          
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) = &
          pcGridIn%parD0(iStart+1:iStart+cSpace%iDimT) * dTaperSupportWindow		
          
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    ! --------------------------
    
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)
    call ReorderDistr0ToDistr1(pcGridIn)
    
    
    call TransformT(pcGrid)
    call TransformT(pcGridIn)
    
    
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(cSpace%iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
          pcGrid%pacD1(iStart+1:iStart+iDimW) = &
          pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))
          
          pcGridIn%pacD1(iStart+1:iStart+iDimW) = &
          pcGridIn%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * 1.0_dp/(real(pcGrid%iD0TL,dp))		
          
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
    
    call TransformTInv(pcGrid)
    call TransformTInv(pcGridIn)
    
    
    !Calculate the square of p with anti-aliasing regions, but take care of the complex multiplication!!
    !Only real multiplication..
    
    iStart=0
    do iIndex = 0, cSpace.cGrid.iD1LocN-1
       
       arBuffer1square=real(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBuffer2square=dimag(pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       arBufferIn1square=real(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT),dp);
       arBufferIn2square=dimag(pcGridIn%pacD1(iStart+1:iStart+cSpace%iDimT));
       
       pcGrid%pacD1(iStart+1:iStart+cSpace%iDimT) = &
       arBuffer1square*arBufferIn1square + im * arBuffer2square*arBufferIn2square;		
       
       iStart		= iStart + cSpace.cGrid.iD1IS;
       
    end do
    
    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(pcGridIn)
    
    ! Now pcGrid contains p^2 in space-time domain in Distribution 0
    !-----------------------------------------------------------------------------------------
    
    ! now we have to calculate dxkappa and multiply each term with p^2 
    
    ALLOCATE(ContrastFiltered(iDimX,iDimY,iDimZ))
    
    ALLOCATE(KAPPAContrast(iDimX,iDimY,iDimZ))
    
    ContrastFiltered=0
    KAPPAContrast=0
    call ContrastGenerationandFiltering(iDimX,iDimY,iDimZ,ContrastFiltered) !! L.D. 03-11-2009
    istart=0
    
    !----------Distr0
    do xindex=1,iDimX
       do yindex=1,iDimY
          do zindex=1,iDimZ
             
             if (ContrastFiltered(xindex,yindex,zindex)==1) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA1 ! Blood
                
             elseif (ContrastFiltered(xindex,yindex,zindex)==2) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA2 ! Brain
                
             elseif (contrastfiltered(xindex,yindex,zindex)==3) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA3 ! Breast
                
             elseif (contrastfiltered(xindex,yindex,zindex)==4) then
                
                KAPPAcontrast(xindex,yindex,zindex)=cSourceParams.KAPPA4 ! Heart
                
             elseif (contrastfiltered(xindex,yindex,zindex)==5) then
                
                KAPPAcontrast(xindex,yindex,zindex)= cSourceParams.KAPPA5 ! Liver
                
             else
                
                KAPPAcontrast(xindex,yindex,zindex)=1.0_dp/(cMediumParams.rho0*(cMediumParams.c0**2))
                
             end if
             
          end do
       end do
    end do
    
    
    ! now KAPPAcontrast contains the kappa values relative to every medium
    
    !-----------------------------------------------------------------------------------------
    ! Enable this part if you want to write out the contrastkappa matrix
    
    !call MPI_BARRIER(MPI_COMM_WORLD, iErr); !Wait until each processor has checked in
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPA.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrast
    !    CLOSE(3, STATUS = 'keep')  
    !    end if
    !-----------------------------------------------------------------------------------------
    
    DEALLOCATE(ContrastFiltered)
    !------------------------------------------------------------------------- KAPPA DERIVATIVES START
    ! Initialize
    !call DerivLookupInit(cModelParams%FDTOrder, iFDMinOrder, cDerivLookup, .false.)
    
    ALLOCATE(KAPPAContrastDZ(iDimX,iDimY,iDimZ))
    ! derivative with respect to Y
    iLen	= iDimZ
    do xindex=1,iDimX
       do yindex=1,iDimY
          
          
          arBuffer1z=KAPPAcontrast(xindex,yindex,1:iDimZ)
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector    
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          !                    if ((iProcID==1).AND.(yindex==int(iDimY/2)).AND.(zindex==int(iDimZ/2))) then
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\dkVector.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) dKvector
          !                    CLOSE(3, STATUS = 'keep')
          !                                        
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1x.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1x
          !                    CLOSE(3, STATUS = 'keep') 
          !                    	                                   
          !	                OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer1xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer1xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\arBuffer2xC.bin',STATUS='unknown',FORM='binary')
          !                    WRITE(3) arBuffer2xC
          !                    CLOSE(3, STATUS = 'keep') 
          !                    
          !                    end if 
          
          KAPPAcontrastDZ(xindex,yindex,1:iDimZ)=(arBuffer1zC)*(1.0_dp/iLen);
          
       end do
    end do
    
    DEALLOCATE(KAPPAcontrast)
    
    !-----------------------------------------------------------------------------------------
    ! This part enables the export of the matrices containing the derivatives of the kappacontrast
    
    !	if (iProcID==1) then
    !	OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTKAPPADZ.bin',STATUS='unknown',FORM='binary')
    !    WRITE(3) KAPPAcontrastDZ
    !    CLOSE(3, STATUS = 'keep') 
    !    end if
    
    !-----------------------------------------------------------------------------------------
    iStart=0
    ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
    do iIndex=0,pcGrid%iD0LocN-1
       pcGrid%parD0(iStart+1:iStart+cSpace%iDimT) = pcGrid%parD0(iStart+1:iStart+cSpace%iDimT)*KAPPAcontrastDZ(pcgrid%aid0loc(iindex+1,1)+1,pcgrid%aid0loc(iindex+1,2)+1,pcgrid%aid0loc(iindex+1,3)+1)
       iStart=iStart+pcGrid%iD0IS
    end do
    
    !	! Now pcGrid contains dxk p^2 and we have to perform the derivative of this function
    DEALLOCATE (KAPPAcontrastDZ)
    !
    ! we go from ditr 0 to 1
    call ReorderDistr0ToDistr1(pcGrid)
    ! we go from ditr 1 to 2
    call ReorderDistr1ToDistr2(pcGrid)
    
    iStart2=0
    iStart=0
    Timestart=0
    iLen=pcGrid%iD2ZL
    
    do tindex=0, pcGrid%iD2LocN-1
       
       Timestart=0+(pcGrid%iD2XL*pcGrid%iD2YL*pcGrid%iD2ZL)*tindex
       
       do xindex=0, (pcGrid%iD2XL*pcGrid%iD2YL)-1
          
          iStart=0+xindex+Timestart
          iStart2=0+xindex+Timestart
          
          
          ! Now we calculate dz(dzk p^2) so that afterwards we can apply the derivative to this term
          
          
          !            do iIndex2=0,pcGrid%iD2ZL-1
          !        
          !			    pcGrid%pacD2(iStart+1) = pcGrid%pacD2(iStart+1)*dKvector(iIndex2+1)*(cSpace%dDx**3/( real(2*cSpace%iDimT * cSpace%iDimX,dp) * cSpace%iDimY * cSpace%iDimZ))
          !			    
          !			
          !			    iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
          !					
          !		    end do
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points
             
             arBuffer1z(iIndex2+1) = pcGrid%pacD2(iStart+1)	
             
             iStart=iStart+pcGrid%iD2XL*pcGrid%iD2YL
             
          end do
          
          ! fft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer1z,arBuffer2zC,fftw_forward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          arBuffer2zC=arBuffer2zC*dKvector ! Here the transformed vector is multiplied with the k vector   
          
          ! ifft
          call dfftw_plan_dft_1d(plan1d,iLen,arBuffer2zC,arBuffer1zC,fftw_backward,fftw_estimate+fftw_unaligned)
          call dfftw_execute(plan1d)
          call dfftw_destroy_plan(plan1d)
          
          
          do iIndex2=0,pcGrid%iD2ZL-1 ! This is the loop for z points to put the outcome into the pcGrid
             
             pcGrid%pacD2(iStart2+1)=arBuffer1z(iIndex2+1)*(1.0_dp/iLen) 
             
             iStart2=iStart2+pcGrid%iD2XL*pcGrid%iD2YL
             
             
          end do
          
          
          
          
       end do
    end do
    
    
    call ReorderDistr2ToDistr1(pcGrid)
    call ReorderDistr1ToDistr0(pcGrid)
    
    
    
  END SUBROUTINE NonlinearKappaOperatorLinDZALi
  !-Added 
  
    SUBROUTINE LagrangianDensity_Ali(cSpace)
  
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. The spatial 
    !   derivative is implemented in kx space.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace 
	type(space), target :: phiSliceMirror, phiSlice,cSpaceTemp
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
	type(Grid), pointer                          :: pcGrid
	type(Space), target 						 ::	phi
    type(DerivLookup) , target                   :: cDerivLookup
	
    integer(i8b)                                 :: plan1D, plan1D_inv
    integer(i4b)                                 :: iErr;
    character(len=1024)                          :: acTemp, filename
	
    integer(i8b)                                 :: iDimX, iDimY, iDimZ, iDimT, iDimW, iDimX_Wrap, iDimY_Wrap, iDimZ_Wrap
	integer(i8b)								 :: iLx  , iLy  , iLz  , iLt  , i, iStart, Timestart, iIndex , YBlock
	
    integer(i8b)                      			 :: wrap_error_padX , wrap_error_padY, wrap_error_padZ, wrap_error_padFD, wrap_error_pad_MirrorY
    
    complex(dpc), allocatable                    :: arBuffer(:) , dKvectorX(:), dKvectorY(:), dKvectorZ(:)
    
	real(dp)    , allocatable  					 :: dKvectorRealX(:) , dKvectorRealY(:) , dKvectorRealZ(:)
	real(dp)    , allocatable			         :: arBuffer1square(:),arBuffer2square(:)
    
    complex(dpc), allocatable   		      	 :: dKVectorRealXYZ(:), dMultFactor(:)
    
    real(dp), allocatable                        :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:),dTaperMaxFreqWindowP(:)
    real(dp)                                     :: dLeftBand, dRightBand, dFFTFactor, dDOmega
    
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("Compute the Lagrangian Density",1)
    
    pcGrid=>cSpace%cGrid
    
	call InitSpace(cSpaceTemp, cSpace%iSpaceIdentifier, cSpace%bYSymm, &
					cSpace%iDimT,cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ, &
					cSpace%iStartT,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ, &
					0_i8b,0_i8b,0_i8b, &
					cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	call InitGrid(cSpaceTemp, pcGrid%iProcN, pcGrid%iProcID, (/ .false., .false., .false., .false./));
    call GridDistr0CreateEmpty(cSpaceTemp%cGrid);
    cSpaceTemp%cGrid%parD0  = pcGrid%parD0
    
	call InitSpace(phi, cSpace%iSpaceIdentifier, cSpace%bYSymm, &
					cSpace%iDimT,cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ, &
					cSpace%iStartT,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ, &
					0_i8b,0_i8b,0_i8b, &
					cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	call InitGrid(phi, pcGrid%iProcN, pcGrid%iProcID, (/ .false., .false., .false., .false./));
    call GridDistr1CreateEmpty(phi%cGrid);

    iStart	= 0;
    iDimT   = cSpace%iDimT      !time dimensions
    iDimX   = phi%cGrid%iD2XL   !x dimensions
    iDimY   = phi%cGrid%iD2YL   !y dimensions
    iDimZ   = phi%cGrid%iD2ZL   !z dimensions
    iDimW   = iDimT/2 + 1
    
    ! Here I use these variables in order to extend the domain for the spatial derivative
    ! Moreover I take into consideration the case where the Y symmetry is used in the domain
	wrap_error_padX = 0 ; wrap_error_padY = 0 ; wrap_error_padZ = 0; wrap_error_pad_MirrorY = 0;
	if (cSpace%bYSymm .EQV. .TRUE.)  wrap_error_pad_MirrorY = iDimY-1
    iDimX_Wrap = iDimX + wrap_error_padX   ! This is because zero padding is needed in order for the fft to be correct
    iDimY_Wrap = iDimY + wrap_error_padY  + wrap_error_pad_MirrorY
    iDimZ_Wrap = iDimZ + wrap_error_padZ 
    
	ALLOCATE(  dKvectorX(iDimX_Wrap), dKvectorRealX(iDimX_Wrap))
	ALLOCATE(  dKvectorY(iDimY_Wrap), dKvectorRealY(iDimY_Wrap))
	ALLOCATE(  dKvectorZ(iDimZ_Wrap), dKvectorRealZ(iDimZ_Wrap))
    
    dKvectorRealX = 0.0D0;
    dKvectorRealY = 0.0D0;
    dKvectorRealZ = 0.0D0;
	
    call PrintToLog("Compute the K vector",2)
	! First half is positive, the other half is negative
    dKvectorRealX(1:iDimX_Wrap) = (/(iLx, iLx = 0, iDimX_Wrap/2), (iLx - iDimX_Wrap, iLx = iDimX_Wrap/2+1,iDimX_Wrap - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimX_Wrap,dp))
    dKvectorX = im * dKvectorRealX ;
	
    ! Create the K vector in Ky space
	dKvectorRealY(1:iDimY_Wrap) = (/(iLy, iLy = 0, iDimY_Wrap/2), (iLy - iDimY_Wrap, iLy = iDimY_Wrap/2+1,iDimY_Wrap - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimY_Wrap,dp))
    dKvectorY = im * dKvectorRealY  ;           ! Normalization factor λ in [mm];
	
    ! Create the K vector in Kz space 
	dKvectorRealZ(1:iDimZ_Wrap) = (/(iLz, iLz = 0, iDimZ_Wrap/2), (iLz - iDimZ_Wrap, iLz = iDimZ_Wrap/2+1,iDimZ_Wrap - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimZ_Wrap,dp))
    dKvectorZ = im * dKvectorRealZ  ;           ! Normalization factor λ in [mm];
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(iDimT), dTaperMaxFreqWindow(iDimW),dMultFactor(iDimW),dTaperMaxFreqWindowP(iDimW))
	
	write(*,*) "Pressure,",MAXVAL(REAL(pcGrid%parD0)),  cSpace%cGrid%iProcID
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       ! use tapering of the field at start and end to prevent wraparound leakage
       ! build tapering window, the first two periods and the last two periods are tapered
       ! in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       ! multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+iDimT) = &
          pcGrid%parD0(iStart+1:iStart+iDimT) * dTaperSupportWindow
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
	
    ! --------------------------
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)  
	phi%cGrid%pacD1 = pcGrid%pacD1
    
    call TransformT_sml(pcGrid)  
    call TransformT_sml(phi%cGrid)  
	! Do this for both phi and cSpace in order to have the same dimensions 
	!============================================================================================================================================
    call PrintToLog("Compute the velocity potential",2)
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       ! Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! and is blown up by the double derivative...
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0.0_dp
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
	
	   dFFTFactor = 1.0_dp/real(pcGrid%iD0TL,dp)
	   dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(iDimT,dp) 
	   dMultFactor = - (/ (i,i=0,iDimW -1 ) /) * dDOmega * im
	   
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) = &
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) /dMultFactor /cMediumParams%rho0 * dTaperMaxFreqWindow * dFFTFactor 
		   phi%cGrid%pacD1(iStart+1) =0.0D0
           
		   pcGrid%pacD1(iStart+1:iStart+iDimW) = pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * dFFTFactor 
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
	  
    call TransformTInv(phi%cGrid)
    call TransformTInv(pcGrid)   
	pcGrid%pacD1 = pcGrid%pacD1/real(cMediumParams%rho0*cMediumParams%c0,dp)

    call ReorderDistr1ToDistr2(phi%cGrid)
    call ReorderDistr1ToDistr2(pcGrid)
    
	!============================================================================================================================================
    call PrintToLog("Multiply the velocity potential with the kx",2)

    ! here we perform the fft transform, multply with Kx and go back to time domain. 
	! Normalization has to be performed with respect to the forward fft
	! Remember that pacD2 has the form of [XYZ1T1 XYZ2T1 XYZ3T1 ... XYZNT1 XYZ1T2 ... XYZNT2 ... XYZ1TN XYZ2TN ... XYZNTN]
	! And you can find XYZ by = x + y * iDimX + z * iDimX * iDimY
    
    call DerivLookupInit(cModelParams%FDXOrder, iFDMinOrder, cDerivLookup, .false.)
	iStart = 0 ;
	ALLOCATE( arBuffer1square(iDimX_Wrap), arBuffer2square(iDimX_Wrap), arBuffer(iDimX_Wrap) )
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D     , iDimX_Wrap, arBuffer, arBuffer , fftw_forward , fftw_estimate)
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D_inv , iDimX_Wrap, arBuffer, arBuffer , fftw_backward, fftw_estimate)

    ! SpatialDerivative_w_wo_Ali functions gives the option to choose between FFT and FD Scheme
    ! Initially validation was tested and in order to be accurate, an agreement between those
    ! methods should have been checked.
    do iIndex = 0,pcGrid%iD2LocSize/pcGrid%iD2XL-1
       arBuffer = phi%cGrid%pacD2(iStart+1:iStart+iDimX)
       
        ! Each complex value stores two values of t, so we have to differentiate them both
        ! As differentiation is a simple scalar multiplication and addition,
        ! this happens automatically (real -> real, imag -> imag)
       call SpatialDerivative_w_wo_Ali(phi, cDerivLookup, arBuffer, dKvectorX, .false., 1)
	   
	   arBuffer1square = real(pcGrid%pacD2(iStart+1:iStart+iDimX),dp)
	   arBuffer2square = dimag(pcGrid%pacD2(iStart+1:iStart+iDimX))
	   
       pcGrid%pacD2(iStart+1:iStart+iDimX) = -(arBuffer1square**2 + im * arBuffer2square**2) + (real(arBuffer,dp)**2 + im * dimag(arBuffer)**2)
	   
       iStart=iStart + phi%cGrid%iD2YS
    end do
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D )
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D_inv )
	DEALLOCATE( arBuffer, arBuffer1square, arBuffer2square)
    
    call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "X,",MAXVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "X_MIN,",MINVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	
    ! ! ! ============================= COMPUTE THE Y COMPONENT ============================================================
	
    call PrintToLog("Multiply the velocity potential with the ky",2)
    
	ALLOCATE(   arBuffer(iDimY_Wrap))
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D     , iDimY_Wrap, arBuffer, arBuffer , fftw_forward , fftw_estimate)
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D_inv , iDimY_Wrap, arBuffer, arBuffer , fftw_backward, fftw_estimate)
    arBuffer = 0;
	
    do iLt = 0, phi%cGrid%iD2LocN-1 ! This is the loop for time instants stored locally
	
       ! This is the index that takes into account for the different starting point corresponding to a given time index
       Timestart  = 0 + iDimX * iDimY * iDimZ * iLt
	   
       do iLx = 0, iDimX-1 ! This is the loop for x points
          iStart  = 0 + iLx + Timestart  ! This is the index
          
          ! Now we calculate dxk p^2 so that afterwards we can apply the derivative to this term
          do iIndex = 0,iDimZ-1 ! This is the loop for z points
			 arBuffer(1:iDimY)   = phi%cGrid%pacD2( iStart + 1 : iStart + phi%cGrid%iD2ZS : phi%cGrid%iD2YS)
             if (cSpace%bYSymm == .true.) then;  arBuffer(iDimY+1:iDimY_Wrap) = arBuffer(2:iDimY); arBuffer(1:iDimY) = arBuffer(iDimY_Wrap:iDimY:-1); endif
             call SpatialDerivative_w_wo_Ali(phi, cDerivLookup, arBuffer, dKvectorY, .false., 1)
            
             if (cSpace%bYSymm == .true.) arBuffer(1:iDimY) = arBuffer(iDimY:iDimY_Wrap) 
             pcGrid%pacD2(iStart + 1 : iStart + phi%cGrid%iD2ZS : phi%cGrid%iD2YS) = &
			 pcGrid%pacD2(iStart + 1 : iStart + phi%cGrid%iD2ZS : phi%cGrid%iD2YS) + (real(arBuffer(1:iDimY),dp)**2 + im * dimag(arBuffer(1:iDimY))**2) 
														
             iStart		= iStart 	+	phi%cGrid%iD2ZS
             
          end do
          
          
       end do
	enddo
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D )
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D_inv)
	DEALLOCATE( arBuffer )
	
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "Y,",MAXVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "Y_MIN,",MINVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	  
    
    ! ! ! ============================= COMPUTE THE Z COMPONENT ============================================================
    call PrintToLog("Multiply the velocity potential with the kz",2) 	
      
	ALLOCATE(   arBuffer(iDimZ_Wrap)  )
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D     , iDimZ_Wrap, arBuffer, arBuffer , fftw_forward , fftw_estimate + fftw_unaligned)
    call dfftw_plan_dft_1d(phi%cGrid%cTransforms%iPlanTransform1D_inv , iDimZ_Wrap, arBuffer, arBuffer , fftw_backward, fftw_estimate + fftw_unaligned)
    	
	do iLt = 0, phi%cGrid%iD2LocN - 1! This is the loop for time instants stored locally
       
	   ! This is the index that takes into account for the different starting point corresponding to a ginve time index   
	   Timestart = 0 + phi%cGrid%iD2IS * iLt
	   
        do iLx = 0, iDimX * iDimY - 1
            
          iStart  = 0 + iLx + Timestart  ! This is the index for all the possible x,y for the same z
          arBuffer = phi%cGrid%pacD2(iStart+1 : iStart + phi%cGrid%iD2IS : phi%cGrid%iD2ZS)
          ! Compute the 1st spatial derivative  in the z direction
          call SpatialDerivative_w_wo_Ali(phi, cDerivLookup, arBuffer, dKvectorZ, .false., 1)
            
          pcGrid%pacD2(iStart + 1 : iStart + phi%cGrid%iD2IS : phi%cGrid%iD2ZS)  = &
		  pcGrid%pacD2(iStart + 1 : iStart + phi%cGrid%iD2IS : phi%cGrid%iD2ZS) + (real(arBuffer,dp)**2 + im * dimag(arBuffer)**2)
            
        end do   
	enddo
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D )
    call dfftw_destroy_plan(phi%cGrid%cTransforms%iPlanTransform1D_inv )
	DEALLOCATE( arBuffer )
    
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "Lagrangian,",MAXVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	write(*,*) "Lagrangian_MIN,",MINVAL(REAL(pcGrid%pacD2,dp)),  cSpace%cGrid%iProcID
    
	pcGrid%pacD2 = pcGrid%pacD2 *cMediumParams%rho0/2.0D0 
	phi%cGrid%pacD2 = pcGrid%pacD2 ! This is the Lagrangian Density
    
	call ReorderDistr2ToDistr1(pcGrid)
	
	call TransformT(pcGrid)   
	
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       ! Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       ! build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0.0_dp
       dRightBand = 0.1_dp
	   
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
       dFFTFactor = 1.0_dp/(2.0D0 * real(pcGrid%iD0TL,dp) )
       ! multiply the field with it, normalization with respect to the forward tranformation is performed 
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           pcGrid%pacD1(iStart+1:iStart+iDimW) = &
           pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * dFFTFactor
		   
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
	
    call TransformTInv_sml(pcGrid) 
	
	call ReorderDistr1ToDistr0(pcGrid)

	do i=1,cModelParams%numslices
        if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            filename = trim(trim(sOutputDir) // trim('Lagrangian') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                int2str(i)//int2str(0)
                
            if (cModelParams%xyzslicedim(i)=='t') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            elseif (cModelParams%xyzslicedim(i)=='x') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            elseif (cModelParams%xyzslicedim(i)=='y') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            elseif (cModelParams%xyzslicedim(i)=='z') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            end if
        end if
    end do
	
	call ReorderDistr0ToDistr1(pcGrid)
    
    call TransformT_sml(pcGrid)   
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       ! Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! and is blown up by the double derivative...
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0.0_dp
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
	
	   dFFTFactor = 1.0_dp/real(pcGrid%iD0TL,dp)
	   
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
		   pcGrid%pacD1(iStart+1:iStart+iDimW) = pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * dFFTFactor 
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
	  
    call TransformTInv(pcGrid) 
	 
	!============================================================================================================================================
    
    call PrintToLog("Allocate Slice Space", 3)
	call InitSpace(phiSliceMirror, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
						cSpace%cGrid%iProcN-1, iDimX_Wrap, iDimY_Wrap, iDimZ_Wrap, &
						0_i8b,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ,0_i8b,0_i8b,0_i8b, &
						cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	!iDimT is such that each processor gets one complex, i.e. two real slices..
    call InitGrid(phiSliceMirror, cSpace%cGrid%iProcN, cSpace%cGrid%iProcID, (/ .false., .false., .false., .false./));
	call GridDistr2CreateEmpty(phiSliceMirror%cGrid);

	call PrintToLog("Allocated", 3)

    call PrintToLog("Compute the Laplacian Operator of Lagrangian Density",2)
    ALLOCATE(dKVectorRealXYZ(phiSliceMirror%cGrid%iD1GlobN))
	do iLz = 0, iDimZ_Wrap-1
		do iLy = 0, iDimY_Wrap-1
			do iLx = 0, iDimX_Wrap-1
				dKVectorRealXYZ(iLx + iLy  * iDimX_Wrap + iLz * iDimX_Wrap * iDimY_Wrap + 1) =  (dKvectorRealX(iLx+1)**2 + dKvectorRealY(iLy+1)**2 + dKvectorRealZ(iLz+1)**2 ) 
			enddo
		enddo
	enddo 
    
    call DerivLookupDestroy(cDerivLookup);  
    call DerivLookupInit(cModelParams%FDXOrder, iFDMinOrder, cDerivLookup, .true.)
    
	do iLt = 0, phi%cGrid%iD2LocN-1 ! This is the loop for time instants stored locally
		
		write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iLt,  phi%cGrid%iD2LocN-1
		call PrintToLog(acTemp, 3);
		call PrintToLog("Obtain Mirrored Field slice", 4)
        ! If symmetry is used, the domain is only the half of the total, so a filling of the rest of the array should take place
        if (cSpace%bYSymm .EQV. .TRUE.) then
            call Distr2FillYMirroredXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        else
            ! Implement this in order to fill the array with the values of the current domain
            call Distr2ObtainXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        endif
        
        ! Compute the 2nd spatial derivative either using FFT or with FD
        call SpatialDerivative_w_wo_Ali(phiSliceMirror, cDerivLookup, arBuffer, dKVectorRealXYZ, .true., 2)
         
        call PrintToLog("Put Field slice", 4)
        ! Return the part of the array that is of interest to the initial phi space.
         if (cSpace%bYSymm .EQV. .TRUE.) then
            call Distr2PutYMirroredXYZBlock(phiSliceMirror%cGrid, phi%cGrid,  iLt)
        else
            ! Implement this in order to fill the array with the values of the current domain
            call Distr2PutXYZBlock(phiSliceMirror%cGrid, phi%cGrid,  iLt)
        endif
	enddo
    call DerivLookupDestroy(cDerivLookup);
    call DestructSpace(phiSliceMirror)
    
    DEALLOCATE(dKVectorRealXYZ)
    
	call ReorderDistr2ToDistr1(phi%cGrid) 
	
    call PrintToLog("Compute the time derivative of Lagrangian Density",2)

    call TransformT(pcGrid)    
    call TransformT(phi%cGrid)    
	
	dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(iDimT,dp)
	dMultFactor = -( (/ (i,i=0,iDimW-1) /) * dDOmega )**2
    dFFTFactor = 1.0_dp/real(2.0D0*cSpace%cGrid%iD0TL,dp)
	
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           pcGrid%pacD1(iStart+1:iStart+iDimW) = &
           pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor * dTaperMaxFreqWindow * dFFTFactor
		   
		   phi%cGrid%pacD1(iStart+1:iStart+iDimW) = &
           phi%cGrid%pacD1(iStart+1:iStart+iDimW)  * dTaperMaxFreqWindow * dFFTFactor 
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
     
    call TransformTInv_sml(pcGrid) 
    call TransformTInv_sml(phi%cGrid) 

    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(phi%cGrid)
	
	write(*,*) "time, end",MAXVAL(REAL(pcGrid%parD0,dp))/real(cMediumParams%c0**2,dp), pcGrid%iProcID
	write(*,*) "space, end",MAXVAL(REAL(phi%cGrid%parD0,dp)), pcGrid%iProcID
	
	pcGrid%parD0 =  pcGrid%parD0/real(cMediumParams%c0**2,dp) + phi%cGrid%parD0 ! Addition of the two terms in the lagrangian ( Kinetic and potential) 
    cSpace%cGrid%parD0 =  pcGrid%parD0 * real(cMediumParams%c0**2,dp) ! Normalization factor
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	
	call DestructSpace(phi) 
	DEALLOCATE( dKvectorX, dKvectorRealX)
	DEALLOCATE( dKvectorY, dKvectorRealY)
	DEALLOCATE( dKvectorZ, dKvectorRealZ)
    DEALLOCATE(dTaperSupportWindow,dTaperMaxFreqWindow,dMultFactor)
    
	write(*,*) "pc, end",MAXVAL(pcGrid%parD0), cSpace%cGrid%iProcID
    
	do i=1,cModelParams%numslices
        if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            filename = trim(trim(sOutputDir) // trim('LagrangianPressure') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                int2str(i)//int2str(0)
                
            if (cModelParams%xyzslicedim(i)=='t') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            elseif (cModelParams%xyzslicedim(i)=='x') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            elseif (cModelParams%xyzslicedim(i)=='y') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            elseif (cModelParams%xyzslicedim(i)=='z') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            end if
        end if
    end do
	!============================================================================================================================================
    ! call NonlinContrastOperator_Ali(cSpaceTemp);  
    ! pcGrid%parD0 = pcGrid%parD0 + cSpaceTemp%cGrid%parD0
    call DestructSpace(cSpaceTemp)
    
  END SUBROUTINE LagrangianDensity_Ali

  SUBROUTINE SpatialDerivative_w_wo_Ali(cSpace, cDerivLookup, arBuffer, dKVectorXYZ, Alias, Order)
  
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. The spatial 
    !   derivative is implemented in kx space.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)          :: cSpace 
    type(DerivLookup), intent(in)               :: cDerivLookup
    complex(dpc), intent(inout)                 :: arBuffer(:) , dKVectorXYZ(:)
    logical(lgt), intent(in)                    :: Alias
    integer(i8b), intent(in)                    :: Order
    
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
		
    real(dp)                                     :: dFFTFactor
    integer(i8b)                                 :: iLenXYZ
    
    iLenXYZ = size(arBuffer) ; 
    dFFTFactor	= 1.0D0/real(cSpace%cGrid%iD1GlobN,dp);	
    if (Alias) then
        if(Order == 1) then        
            call dfftw_execute(cSpace%cGrid%cTransforms%iPlanTransform1D,arBuffer,arBuffer)
            arBuffer = arBuffer * dKVectorXYZ / real(iLenXYZ,dp)
            call dfftw_execute(cSpace%cGrid%cTransforms%iPlanTransform1D_inv ,arBuffer,arBuffer)
        else
            call PrintToLog("Transform Field slice in XYZ", 4)
            call TransformXYZ(cSpace%cGrid, .true.)
            call PrintToLog("Multiply Result with k vector", 4)
            cSpace%cGrid%pacD2 = -dKVectorXYZ * cSpace%cGrid%pacD2 * dFFTFactor
            call PrintToLog("Inverse transform Field slice in XYZ", 4)
            call TransformXYZInv(cSpace%cGrid, .true.)
        endif
    else
        if(Order == 1) then
            call DerivativeComplex(arBuffer/cMediumParams%c0,arBuffer,iLenXYZ, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
        else
            call DerivativeComplex(cSpace%cGrid%pacD2/cMediumParams%c0**2,cSpace%cGrid%pacD2,cSpace%cGrid%iD1GlobN, cSpace.dDx**2, cDerivLookup.arWeights, cDerivLookup.aiPoints)
        endif
    endif
    
    END SUBROUTINE SpatialDerivative_w_wo_Ali
  
  SUBROUTINE LagrangianDensity_Simpl_Ali(cSpace)
    
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
    !Commented for Linux  !DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
    
    ! =============================================================================
    !
    !   Programmer: Libertario Demi
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.0     16022012  Original code (LD)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine NonlinContrastOperator computes the nonlinearity contrast
    !   source S^nlkappa(p) = dx(dx(kappa)p^2) for the given space. The 
    !   data in cSpace should be in Distribution 0. The spatial 
    !   derivative is implemented in kx space.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace   io   type(space)  space for which the contrast source
    !                              is determined
    !
    type(Space), target, intent(inout)::	cSpace
	type(space), target :: phiSliceMirror, phiSlice,cSpaceTemp
    
    ! *****************************************************************************
    ! 
    !   LOCAL PARAMETERS      
    !
    !   (x,y,z,t)index      i8b   loop counter over the positions in the space(time)
    !   iLen                i8b   length of a space(x,y,z) trace
    !   iStart              i8b   start of each space trace
    !   iindex              i8b   loop counter
    !   cDerivLookup  type(DerivLookup)  structure that contains the weights of the
    !                       Finite Difference stencil
    !   arBuffer1           dp    Temporary buffer containing a spatial(x,y,z) trace
    !   arBuffer2           dp    Temporary buffer containing a spatial(x,y,z) trace
    !
	type(Grid), pointer                          :: pcGrid
	type(Space), target 						 ::	phi
	
    INTEGER(8)                                   :: plan1D, plan1D_inv, plan3D, plan3D_INV
    integer(i4b)                                 :: iErr;
    character(len=1024)                          :: acTemp, filename
    type(DerivLookup)                            :: cDerivLookup
    
    
    integer(i8b)                                 :: iDimX, iDimY, iDimZ, iDimT, iDimW, iDimX_Wrap, iDimY_Wrap, iDimZ_Wrap
	integer(i8b)								 :: iIndex, iLx, iLy, iLz, iLt, iIndD, i, iStart, iStart2, Timestart, iIndex2
    integer(i8b)                      			 :: wrap_error_padX , wrap_error_padY, wrap_error_padZ, wrap_error_padFD, wrap_error_pad_MirrorY
	
	complex(dpc),allocatable					 :: arBuffer(:),dMultFactor(:)
    complex(dpc), allocatable   		      	 :: arBuffer1XYZ(:), arBuffer1XYZC(:), arBuffer2XYZC(:) , dKVectorRealXYZ(:)
	
	real(dp), dimension(cSpace.iDimT)            :: arBuffer1square,arBuffer2square,arBufferIn1square,arBufferIn2square
	
	complex(dpc)							     :: dKvectorX(cSpace.iDimX)     , dKvectorY(cSpace.iDimY)     , dKvectorZ(cSpace.iDimZ) 
	real(dp)									 :: dKvectorRealX(cSpace.iDimX) , dKvectorRealY(cSpace.iDimY) , dKvectorRealZ(cSpace.iDimZ) 
	real(dp)									 :: dTaperX(cSpace.iDimX) , dTaperY(cSpace.iDimY) , dTaperZ(cSpace.iDimZ) 
    
    real(dp), allocatable                        :: dTaperSupportWindow(:),dTaperMaxFreqWindow(:)
    real(dp)                                     :: dLeftBand, dRightBand, dFFTFactor, dDOmega, dLambdaMM
	complex(dpc) :: cFFTWtestarray(1),cFFTWtestarray2(1)
	
    
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
    !   DerivLookupInit
    !   DerivativeReal
    !   DerivLookupDestroy  
    ! *****************************************************************************
    
    
    call PrintToLog("Compute the Lagrangian Density using third and higher order simplifications",1)
    
    pcGrid=>cSpace%cGrid
    
	call InitSpace(phi, cSpace%iSpaceIdentifier, cSpace%bYSymm, &
					cSpace%iDimT,cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ, &
					cSpace%iStartT,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ, &
					0_i8b,0_i8b,0_i8b, &
					cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	call InitGrid(phi, pcGrid%iProcN, pcGrid%iProcID, (/ .false., .false., .false., .false./));
    call GridDistr1CreateEmpty(phi%cGrid);

    iStart	= 0;
    iDimT   = cSpace%iDimT   !time dimensions
    iDimZ   = cSpace%iDimZ   !z dimensions
    iDimX   = cSpace%iDimX   !x dimensions
    iDimY   = cSpace%iDimY   !y dimensions
    iDimW   = iDimT/2 + 1
    
	wrap_error_padX = 0 ; wrap_error_padY = 0 ; wrap_error_padZ = 0; wrap_error_padFD = 0*cModelParams%FDXOrder/2; wrap_error_pad_MirrorY = 0;
	if (cSpace%bYSymm .EQV. .TRUE.)  wrap_error_pad_MirrorY = iDimY
    iDimX_Wrap = iDimX + wrap_error_padX + 2*wrap_error_padFD  ! This is because zero padding is needed in order for the fft to be correct
    iDimY_Wrap = iDimY + wrap_error_padY + 2*wrap_error_padFD + wrap_error_pad_MirrorY
    iDimZ_Wrap = iDimZ + wrap_error_padZ + 2*wrap_error_padFD
    
	dKvectorRealX = 0.0D0
	dKvectorRealY = 0.0D0
	dKvectorRealZ = 0.0D0
	
	dLambdaMM = (cMediumParams%c0*1.0D3)/cModelParams%freq0
    call PrintToLog("Compute the K vector",2)
    ! Create the K vector in Kx space
    dKvectorRealX=(/(iLx, iLx = 0, iDimX/2), (iLx - iDimX, iLx = iDimX/2 + 1, iDimX - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimX,dp))
    dKvectorX = im * dKvectorRealX ;
	
    ! Create the K vector in Ky space
	dKvectorRealY=(/(iLy, iLy = 0, iDimY/2), (iLy - iDimy, iLy = iDimY/2 + 1, iDimY - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimY,dp))
    dKvectorY = im * dKvectorRealY ;
	
    ! Create the K vector in Kz space 
	dKvectorRealZ=(/(iLz, iLz = 0, iDimZ/2), (iLz - iDimZ, iLz = iDimZ/2 + 1, iDimZ - 1)/)*two_pi*2.0_dp*cSpace%dFnyq/(cMediumParams%c0 * real(iDimZ,dp))
    dKvectorZ = im * dKvectorRealZ  ;           ! Normalization factor λ in [mm];
    ! Calculating p^2
    !-----------------------------------------------------------------------------------------
    allocate(dTaperSupportWindow(iDimT), dTaperMaxFreqWindow(iDimW),dMultFactor(iDimW))
	
	write(*,*) "Pressure,",MAXVAL(REAL(pcGrid%parD0)),  cSpace%cGrid%iProcID
    if (cModelParams.UseSupportTapering .EQV. .true.) then
       ! use tapering of the field at start and end to prevent wraparound leakage
       
       ! build tapering window, the first two periods and the last two periods are tapered
       ! in the contrast source
       dLeftBand = 2.0_dp
       dRightBand = 2.0_dp
       dTaperSupportWindow=dTaperingWindow(iDimT,cSpace%dDt,dLeftBand,dRightBand)
       
       ! multiply the field with it
       iStart = 0
       do iIndex=0,pcGrid%iD0LocN-1
          pcGrid%parD0(iStart+1:iStart+iDimT) = &
          pcGrid%parD0(iStart+1:iStart+iDimT) * dTaperSupportWindow
          iStart=iStart+pcGrid%iD0IS
       end do
       
    end if
    
    ! --------------------------
    !First, transform to W-domain (use small transform, no wraparound regions)
    call ReorderDistr0ToDistr1(pcGrid)  
    phi%cGrid%pacD1  = pcGrid%pacD1 
    call TransformT_sml(phi%cGrid)  
	
	!============================================================================================================================================
    call PrintToLog("Compute the velocity potential",2)
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! and is blown up by the double derivative...
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0.0_dp
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
	
	   dFFTFactor = 1.0_dp/real(pcGrid%iD0TL,dp)
	   dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(iDimT,dp) 
	   dMultFactor = - (/ (i,i=0,iDimW -1 ) /) * dDOmega * im
	   
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) = &
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) /dMultFactor /cMediumParams%rho0 * dTaperMaxFreqWindow * dFFTFactor 
		   phi%cGrid%pacD1(iStart+1) =0.0D0
           iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
	  
    call TransformTInv(phi%cGrid)
    
	iStart=0
    do iIndex = 0, phi.cGrid.iD1LocN-1
       
       arBufferIn1square=real(phi%cGrid%pacD1(iStart+1:iStart+phi.cGrid.iD1IS),dp);
	   arBufferIn2square=dimag(phi%cGrid%pacD1(iStart+1:iStart+phi.cGrid.iD1IS));
	   
       phi%cGrid%pacD1(iStart+1:iStart+phi.cGrid.iD1IS) = arBufferIn1square**2+im*arBufferIn2square**2 ;
       
       iStart		= iStart + phi.cGrid.iD1IS;
       
    end do 
    
    ! Now  phi%cGrid% contains the phi^2
	pcGrid%pacD1 = phi%cGrid%pacD1
    
	call PrintToLog("Compute the Laplacian Operator of Lagrangian Density",2)

	call ReorderDistr1ToDistr2(phi%cGrid)
      
    call PrintToLog("Allocate Slice Space", 3)
	call InitSpace(phiSliceMirror, iSI_XYZSLICE, cModelParams.UseYSymmetry, &
						cSpace%cGrid%iProcN-1, iDimX_Wrap, iDimY_Wrap, iDimZ_Wrap, &
						0_i8b,cSpace%iStartX,cSpace%iStartY,cSpace%iStartZ,0_i8b,0_i8b,0_i8b, &
						cSpace%dFnyq, cSpace%dTanX, cSpace%dTanY, cSpace%dTanT)
	!iDimT is such that each processor gets one complex, i.e. two real slices..
    call InitGrid(phiSliceMirror, cSpace%cGrid%iProcN, cSpace%cGrid%iProcID, (/ .false., .false., .false., .false./));
	call GridDistr2CreateEmpty(phiSliceMirror%cGrid);

	call PrintToLog("Allocated", 3)	
    ALLOCATE(dKVectorRealXYZ(phiSliceMirror%cGrid%iD1GlobN))
    dFFTFactor	= 1.0D0/real(phiSliceMirror%cGrid%iD1GlobN,dp);	
	do iLz = 0, iDimZ_Wrap-1
		do iLy = 0, iDimY_Wrap-1
			do iLx = 0, iDimX_Wrap-1
				dKVectorRealXYZ(iLx + iLy  * iDimX_Wrap + iLz * iDimX_Wrap * iDimY_Wrap + 1) =  (dKvectorRealX(iLx+1)**2 + dKvectorRealY(iLy+1)**2 + dKvectorRealZ(iLz+1)**2 ) 
			enddo
		enddo
	enddo 
    call DerivLookupInit(cModelParams%FDXOrder, iFDMinOrder, cDerivLookup, .true.)
	do iLt = 0, phi%cGrid%iD2LocN-1 ! This is the loop for time instants stored locally
		
		write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iLt,  phi%cGrid%iD2LocN-1
		call PrintToLog(acTemp, 3);
		call PrintToLog("Obtain Mirrored Field slice", 4)
		call Distr2ObtainXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        if (cSpace%bYSymm .EQV. .TRUE.) call Distr2FillYMirroredXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        
        call DerivativeComplex(phiSliceMirror%cGrid%pacD2/cMediumParams%c0**2,phiSliceMirror%cGrid%pacD2,phiSliceMirror%cGrid%iD1GlobN, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
        ! call DerivativeComplex(phiSliceMirror%cGrid%pacD2/cMediumParams%c0,phiSliceMirror%cGrid%pacD2,phiSliceMirror%cGrid%iD1GlobN, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
		! call PrintToLog("Transform Field slice in XYZ", 4)
		! call TransformXYZ(phiSliceMirror%cGrid, .true.)
		! call PrintToLog("Multiply Result with k vector", 4)
		! phiSliceMirror%cGrid%pacD2 = -dKVectorRealXYZ * phiSliceMirror%cGrid%pacD2 * dFFTFactor
		! call PrintToLog("Inverse transform Field slice in XYZ", 4)
		! call TransformXYZInv(phiSliceMirror%cGrid, .true.)
        
        call PrintToLog("Put Field slice", 4)
        call Distr2PutXYZBlock(phiSliceMirror%cGrid, phi%cGrid,  iLt)
		
	enddo
	call ReorderDistr2ToDistr1(phi%cGrid)
	write(*,*) "Space deriv, MIN,",MINVAL(REAL(phi%cGrid%pacD1,dp)),  cSpace%cGrid%iProcID 
	write(*,*) "Space deriv, MAX",MAXVAL(REAL(phi%cGrid%pacD1,dp)),  cSpace%cGrid%iProcID 
	call MPI_BARRIER(MPI_COMM_WORLD,iErr)
	
    call PrintToLog("Compute the time derivative of Lagrangian Density",2)
    call TransformT(pcGrid)
    call TransformT(phi%cGrid)
	
	dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(iDimT,dp)
	dMultFactor = -( (/ (i,i=0,iDimW-1) /) * dDOmega )**2
	dFFTFactor = 1.0D0/(2.0D0*real(cSpace%cGrid%iD0TL,dp) ) ! This is because the 2nd derivative is two times 1st derivative so 1/L * 1/L
	
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0.0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           pcGrid%pacD1(iStart+1:iStart+iDimW) = &
           pcGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow * dMultFactor * dFFTFactor
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) = &
           phi%cGrid%pacD1(iStart+1:iStart+iDimW) * dTaperMaxFreqWindow  * dFFTFactor
		   
          iStart=iStart+pcGrid%iD1IS 
          
       end do
        
    end if
    call TransformTInv_sml(pcGrid)  
    call TransformTInv_sml(phi%cGrid) 
	
	pcGrid%pacD1=pcGrid%pacD1*1.0D0/real(cMediumParams%c0**2,dp)
	! pc now contains the Lagrangian
	write(*,*) "Time deriv, MAX",MAXVAL(REAL(pcGrid%pacD1,dp)),  cSpace%cGrid%iProcID  
	write(*,*) "Time deriv, MIN,",MINVAL(REAL(pcGrid%pacD1,dp)),  cSpace%cGrid%iProcID  
	call MPI_BARRIER(MPI_COMM_WORLD,iErr)
	
	! ! call ReorderDistr1ToDistr0(pcGrid)
	! ! call ReorderDistr1ToDistr0(phi%cGrid)
	
	! ! do i=1,cModelParams%numslices
        ! ! if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            ! ! filename = trim(trim(sOutputDir) // trim('Time_Deriv') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                ! ! int2str(i)//int2str(0)
                
            ! ! if (cModelParams%xyzslicedim(i)=='t') then
                ! ! call ExportSlice(trim(filename),"p",cSpace, &
                    ! ! (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    ! ! (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='x') then
                ! ! call ExportSlice(trim(filename),"p",cSpace, &
                    ! ! (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    ! ! (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='y') then
                ! ! call ExportSlice(trim(filename),"p",cSpace, &
                    ! ! (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    ! ! (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='z') then
                ! ! call ExportSlice(trim(filename),"p",cSpace, &
                    ! ! (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    ! ! (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            ! ! end if
        ! ! end if
    ! ! end do
	
	! ! do i=1,cModelParams%numslices
        ! ! if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            ! ! filename = trim(trim(sOutputDir) // trim('Z_Component') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                ! ! int2str(i)//int2str(0)
                
            ! ! if (cModelParams%xyzslicedim(i)=='t') then
                ! ! call ExportSlice(trim(filename),"p",phi, &
                    ! ! (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    ! ! (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='x') then
                ! ! call ExportSlice(trim(filename),"p",phi, &
                    ! ! (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    ! ! (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='y') then
                ! ! call ExportSlice(trim(filename),"p",phi, &
                    ! ! (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    ! ! (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            ! ! elseif (cModelParams%xyzslicedim(i)=='z') then
                ! ! call ExportSlice(trim(filename),"p",phi, &
                    ! ! (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    ! ! (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    ! ! cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            ! ! end if
        ! ! end if
    ! ! end do 
	
	! ! call ReorderDistr0ToDistr1(pcGrid)
	! ! call ReorderDistr0ToDistr1(phi%cGrid)
    
    pcGrid%pacD1 = phi%cGrid%pacD1 - pcGrid%pacD1
	pcGrid%pacD1 = pcGrid%pacD1*cMediumParams%rho0/4
    
	! pc now contains the Lagrangian
	write(*,*) "Lagrangian,",MAXVAL(REAL(pcGrid%pacD1,dp)),  cSpace%cGrid%iProcID
	write(*,*) "Lagrangian_MIN,",MINVAL(REAL(pcGrid%pacD1,dp)),  cSpace%cGrid%iProcID
	
	call ReorderDistr1ToDistr0(pcGrid)
	do i=1,cModelParams%numslices
        if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            filename = trim(trim(sOutputDir) // trim('Lagrangian') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                int2str(i)//int2str(0)
                
            if (cModelParams%xyzslicedim(i)=='t') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            elseif (cModelParams%xyzslicedim(i)=='x') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            elseif (cModelParams%xyzslicedim(i)=='y') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            elseif (cModelParams%xyzslicedim(i)=='z') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            end if
        end if
    end do
	call ReorderDistr0ToDistr1(pcGrid)
	!============================================================================================================================================
	phi%cGrid%pacD1 = pcGrid%pacD1 ! This is the Lagrangian Density
	call ReorderDistr1ToDistr2(phi%cGrid)
    dFFTFactor	= 1.0D0/real(phiSliceMirror%cGrid%iD1GlobN,dp);	
  
    call PrintToLog("Compute the Laplacian Operator of Lagrangian Density",2)
	do iLt = 0, phi%cGrid%iD2LocN-1 ! This is the loop for time instants stored locally
				
		write (acTemp, '("Start iOmega ", I5 ," out of ", I5)') iLt,  phi%cGrid%iD2LocN-1
		call PrintToLog(acTemp, 3);
		call PrintToLog("Obtain Mirrored Field slice", 4)
		call Distr2ObtainXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        if (cSpace%bYSymm .EQV. .TRUE.) call Distr2FillYMirroredXYZBlock(phi%cGrid, phiSliceMirror%cGrid,  iLt)
        
        ! call DerivativeComplex(phiSliceMirror%cGrid%pacD2/cMediumParams%c0,phiSliceMirror%cGrid%pacD2,phiSliceMirror%cGrid%iD1GlobN, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
        ! call DerivativeComplex(phiSliceMirror%cGrid%pacD2/cMediumParams%c0,phiSliceMirror%cGrid%pacD2,phiSliceMirror%cGrid%iD1GlobN, cSpace.dDx, cDerivLookup.arWeights, cDerivLookup.aiPoints)
		call PrintToLog("Transform Field slice in XYZ", 4)
		call TransformXYZ(phiSliceMirror%cGrid, .true.)
		call PrintToLog("Multiply Result with k vector", 4)
		phiSliceMirror%cGrid%pacD2 = -dKVectorRealXYZ * phiSliceMirror%cGrid%pacD2 * dFFTFactor
		call PrintToLog("Inverse transform Field slice in XYZ", 4)
		call TransformXYZInv(phiSliceMirror%cGrid, .true.)
        
        call PrintToLog("Put Field slice", 4)
        call Distr2PutXYZBlock(phiSliceMirror%cGrid, phi%cGrid,  iLt)
		
	enddo
    call DestructSpace(phiSliceMirror)
    DEALLOCATE(dKVectorRealXYZ)
    
	call ReorderDistr2ToDistr1(phi%cGrid) 
	
    call PrintToLog("Compute the time derivative of Lagrangian Density",2)

    call TransformT_sml(pcGrid)    
    call TransformT_sml(phi%cGrid)    
	
	dDOmega = two_pi * 2.0_dp * cSpace%dFnyq / real(iDimT,dp)
	dMultFactor = -( (/ (i,i=0,iDimW-1) /) * dDOmega )**2
    dFFTFactor = 1.0_dp/real(cSpace%cGrid%iD0TL,dp)
	
    if (cModelParams.UseFreqTapering .EQV. .true.) then
       !Tapering of the highest frequency part; otherwise the chopoff noise around the 
       ! Nyquist frequency is being distributed over the entire frequency axis by the p**2
       ! and is blown up by the double derivative...
       
       !build tapering window, the highest .1*f0 are tapered in the contrast source
       dLeftBand = 0
       dRightBand = 0.1_dp
       dTaperMaxFreqWindow=dTaperingWindow(iDimW,1.0_dp/(iDimT*cSpace%dDt),dLeftBand,dRightBand)
       
       !multiply the field with it, normalization with respect to the forward tranformation is performed
       iStart = 0
       do iIndex=0,pcGrid%iD1LocN-1
          
           pcGrid%pacD1(iStart+1:iStart+iDimW) = &
           pcGrid%pacD1(iStart+1:iStart+iDimW) * dMultFactor * dTaperMaxFreqWindow * dFFTFactor
		   
		   phi%cGrid%pacD1(iStart+1:iStart+iDimW) = &
           phi%cGrid%pacD1(iStart+1:iStart+iDimW)  * dTaperMaxFreqWindow * dFFTFactor 
          iStart=iStart+pcGrid%iD1IS
          
       end do
       
    end if
    
    !--------------------------------
    !Now, transform back to T-domain as though it was a grid with wraparound regions
    !however, the wraparound regions are now used as anti-aliasing regions...
     
    call TransformTInv_sml(pcGrid) 
    call TransformTInv_sml(phi%cGrid) 

    call ReorderDistr1ToDistr0(pcGrid)
    call ReorderDistr1ToDistr0(phi%cGrid)
	
	write(*,*) "time, end",MAXVAL(REAL(pcGrid%parD0,dp))/real(cMediumParams%c0**2,dp), pcGrid%iProcID
	write(*,*) "space, end",MAXVAL(REAL(phi%cGrid%parD0,dp)), pcGrid%iProcID
	
	pcGrid%parD0 =  pcGrid%parD0/real(cMediumParams%c0**2,dp) + phi%cGrid%parD0 ! Addition of the two terms in the lagrangian ( Kinetic and potential) 
    cSpace%cGrid%parD0 =  pcGrid%parD0 * real(cMediumParams%c0**2,dp) ! Normalization factor
	call MPI_BARRIER(MPI_COMM_WORLD, iErr)
	
	call DestructSpace(phi) 
    DEALLOCATE(dTaperSupportWindow,dTaperMaxFreqWindow,dMultFactor)
    
	write(*,*) "pc, end",MAXVAL(pcGrid%parD0), cSpace%cGrid%iProcID
    
	do i=1,cModelParams%numslices
        if ((cModelParams%xyzslicebeam(i)==0).or.(cModelParams%xyzslicebeam(i)==-1)) then
            filename = trim(trim(sOutputDir) // trim('LagrangianPressure') // int2str(cModelParams%iIter))//'_'//cModelParams%xyzslicedim(i)//&
                int2str(i)//int2str(0)
                
            if (cModelParams%xyzslicedim(i)=='t') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b, 0_i8b /), &
                    (/ 1_i8b, cSpace%iDimX, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimT, .true.);
            elseif (cModelParams%xyzslicedim(i)=='x') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b, 0_i8b /), &
                    (/ cSpace%iDimT, 1_i8b, cSpace%iDimY, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimX, .true.);
            elseif (cModelParams%xyzslicedim(i)=='y') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i), 0_i8b /), &
                    (/ cSpace%iDimT, cSpace%iDimX, 1_i8b, cSpace%iDimZ /), &
                    cModelParams.xyzsliceindex, cSpace%iDimY, .true.);
            elseif (cModelParams%xyzslicedim(i)=='z') then
                call ExportSlice(trim(filename),"p",cSpace, &
                    (/ 0_i8b, 0_i8b, 0_i8b, cModelParams%xyzsliceindex(i) /), &
                    (/ cSpace%iDimT, cSpace%iDimX, cSpace%iDimY, 1_i8b /), &
                    cModelParams.xyzsliceindex, cSpace%iDimZ, .true.);
            end if
        end if
    end do
	!============================================================================================================================================
    ! call NonlinContrastOperator_Ali(cSpaceTemp);  
    ! pcGrid%parD0 = pcGrid%parD0 + cSpaceTemp%cGrid%parD0
    ! call DestructSpace(cSpaceTemp)
      
  END SUBROUTINE LagrangianDensity_Simpl_Ali
  
  END MODULE ParnacContrastFunctions
