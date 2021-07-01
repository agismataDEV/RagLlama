! **************************************************************************

! L.D. 03-11-2009 Contrast Generation and Filtering

! **************************************************************************

MODULE ContrastGenFilt

    USE Types
	USE FFTW_cons
	USE ParnacGeneral


IMPLICIT NONE

CONTAINS

! This subroutine generates a function wich is 1 where there is a contrast
! present and 0 where no contrast is present. This function is oversampled
! in order to avoid an aliasing error. After the generation the function
! is filtered in the spatialfrequency-domain and then  trasformed back in
! the space domain and undersampled

SUBROUTINE ContrastGenerationandFiltering(Xdim,Ydim,Zdim,ContrastFiltered)

! DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_plan_dft_1d_"        :: dfftw_plan_dft_1d
! DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_execute_"            :: dfftw_execute   
! DEC$ ATTRIBUTES DLLIMPORT, ALIAS: "dfftw_destroy_plan_"       :: dfftw_destroy_plan
                                        
! Xdim represents the dimension in the X direction
! Zdim represents the dimension in the Z direction
! ContrastFiltered represents the 2-D matrix where the filtered contrast
! function will be written
! ContrastFilteredoversampled is the oversampled version of ContrastFiltered
! ContrastFilteredoversampled2 is the fourier transformed oversampled version 
! of ContrastFiltered
! indexx and indexz are used as an index in the do loop
! Size is an integer to allocate the size of the problem for the FFT and IFFT
! plan1D is an integer needed for the FFT routine
! MaxX and MaxZ are the dimensions of the domain in mm

  INTEGER(8), INTENT(in)                               :: Xdim,Zdim,Ydim
  REAL(8)                                              :: MaxX,MaxY,MaxZ
  INTEGER(4)                                           :: indexx,indexy,indexz, Size   
  INTEGER(8)                                           :: plan1D
  INTEGER(8)                                           :: a,a1,a2,a22,b,b1,b2,b22,a3,a33,b3,b33,b333,M,N,SHIFT,SHIFTORIZZONTAL
  REAL(8)   ,                            intent(inout) :: ContrastFiltered(:,:,:)
  REAL(8)   ,   ALLOCATABLE                            :: ContrastFilteredFromData(:,:,:)
  REAL(8)                                              :: lambda
  COMPLEX(dp),  ALLOCATABLE                            :: ContrastFilteredoversampled(:,:,:)
  COMPLEX(dp),  ALLOCATABLE                            :: ContrastFilteredoversampled2(:,:,:)
  COMPLEX(dp),  ALLOCATABLE                            :: ColumsArray(:)
  REAL(8)    ,  ALLOCATABLE                            :: GAUSSIAN(:)
  character(LEN = 1024)                                :: acTemp;
  
    
  !ALLOCATE (ContrastFilteredoversampled(2*Xdim,2*Zdim))
  ALLOCATE (ContrastFilteredoversampled(Xdim,Ydim,Zdim))
  ALLOCATE (ContrastFilteredoversampled2(2*Xdim,2*Ydim,2*Zdim))
  ContrastFilteredoversampled =0
  ContrastFilteredoversampled2=0
  
  ! Contrast Generation in the Space domain
  
  ! Oversampling
  !  DO indexz=1,2*Zdim
  !      DO indexx=1,2*Xdim
  
            !ContrastFilteredoversampled(indexx,indexz)=0
            
    DO indexz=1,Zdim
        DO indexy=1,Ydim
            DO indexx=1,Xdim
            
            ContrastFilteredoversampled(indexx,indexy,indexz)=0
           
            
              ENDDO
        ENDDO
  ENDDO
  
  !OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASsT.bin',STATUS='unknown',FORM='binary')
  !WRITE(3) ContrastFilteredoversampled
  !CLOSE(3, STATUS = 'keep')  
  !print *,'ContrastDimensions'
  !print *, M
  !print *, N
  
    
!!!!!  !oversampling only the double !!  
!!!!!   size=2*xdim
!!!!!  !perform the fft along x dimension
!!!!!    do indexx=1,2*zdim
!!!!!    
!!!!!        call dfftw_plan_dft_1d(plan1d,size,contrastfilteredoversampled(:,indexx),contrastfilteredoversampled2(:,indexx),fftw_forward,fftw_estimate+fftw_unaligned)
!!!!!        call dfftw_execute(plan1d)
!!!!!        call dfftw_destroy_plan(plan1d)
!!!!!        
!!!!!    enddo
!!!!!    
!!!!!  !perform the filtering along spatial frequency x
!!!!!    contrastfilteredoversampled2(xdim-int(xdim/2):xdim+int(xdim/2),:)=0
!!!!!    
!!!!!    size=2*zdim
!!!!!    allocate (columsarray(2*zdim))
!!!!!  !perform the fft along z dimension
!!!!!    do indexz=1,2*xdim
!!!!!        do indexx=1,2*zdim
!!!!!        columsarray(indexx)=contrastfilteredoversampled2(indexz,indexx)
!!!!!        enddo
!!!!!        call dfftw_plan_dft_1d(plan1d,size,columsarray,columsarray,fftw_forward,fftw_estimate+fftw_unaligned)
!!!!!        call dfftw_execute(plan1d)
!!!!!        call dfftw_destroy_plan(plan1d)
!!!!!        contrastfilteredoversampled2(indexz,:)=columsarray
!!!!!        do indexx=1,2*zdim
!!!!!        contrastfilteredoversampled2(indexz,indexx)=columsarray(indexx)
!!!!!        enddo
!!!!!    enddo
!!!!!    
!!!!!  !perform the filtering along spatial frequency z
!!!!!     contrastfilteredoversampled2(:,zdim-int(zdim/2):zdim+int(zdim/2))=0
!!!!!    
!!!!!    size=2*xdim
!!!!!  !perform the ifft along x dimension
!!!!!    do indexx=1,2*zdim
!!!!!    
!!!!!        call dfftw_plan_dft_1d(plan1d,size,contrastfilteredoversampled2(:,indexx),contrastfilteredoversampled2(:,indexx),fftw_backward,fftw_estimate+fftw_unaligned)
!!!!!        call dfftw_execute(plan1d)
!!!!!        call dfftw_destroy_plan(plan1d)
!!!!!        
!!!!!    enddo
!!!!!   
!!!!!    size=2*zdim
!!!!! 
!!!!!   !perform the ifft along z dimension
!!!!!    do indexz=1,2*xdim
!!!!!        do indexx=1,2*zdim
!!!!!        columsarray(indexx)=contrastfilteredoversampled2(indexz,indexx)
!!!!!        enddo
!!!!!        call dfftw_plan_dft_1d(plan1d,size,columsarray,columsarray,fftw_backward,fftw_estimate+fftw_unaligned)
!!!!!        call dfftw_execute(plan1d)
!!!!!        call dfftw_destroy_plan(plan1d)
!!!!!        contrastfilteredoversampled2(indexz,:)=columsarray
!!!!!        do indexx=1,2*zdim
!!!!!        contrastfilteredoversampled(indexz,indexx)=columsarray(indexx)
!!!!!        enddo
!!!!!    enddo 
!!!!!
!!!!!  
!!!!!   !allocate (contrastfiltered(xdim,zdim))
!!!!!  
!!!!!  do indexz=1,zdim
!!!!!        do indexx=1,xdim
!!!!!            
!!!!!            contrastfiltered(indexx,indexz)=abs(contrastfilteredoversampled(2*indexx,2*indexz))
!!!!!            
!!!!!        enddo
!!!!!  enddo
!!!!!    
!!!!!  !normalizzation due to the ifft
!!!!!    contrastfiltered=contrastfiltered/(4*zdim*xdim)
!!!!!  
!  OPEN(UNIT=3,FILE='\\srv247\hpc_data\Libertario\CONTRASTFILTERED.bin',STATUS='unknown',FORM='binary')
!  WRITE(3) ContrastFiltered
!  CLOSE(3, STATUS = 'keep')  

contrastfiltered=ContrastFilteredoversampled


! READ THE CONTRAST FROM FILE L.D. 15-04-2011
!print *, '----------------------------------------------------------------------------------------'
!print *, '--------------Reading Matrix From Data--------------------------------------------------'
!print *, '----------------------------------------------------------------------------------------'

ALLOCATE (ContrastFilteredFromData(XDim,YDim,Zdim))


! Edited by Agisilaos Matalliotakis 26/01/2020 , initially FILE='.\matrixcontrastinput.bin'
OPEN(UNIT=3,FILE=trim(sInputDir)//'matrixcontrastinput.bin',STATUS='unknown',FORM='binary')
READ(3) ContrastFilteredFromData
CLOSE(3, STATUS = 'keep')

contrastfiltered=ContrastFilteredFromData

if (size(contrastfiltered) /= size(ContrastFilteredFromData)) then
		write (acTemp, '("Error occurred during the execution of ContrastGenerationandFiltering, aborting program")');
		call PrintToLog(acTemp, -1);
		write (acTemp, '("the 3D matrix provided as contrast input has not the same dimensions as the 3D numerical domain in use")');
		call PrintToLog(acTemp, -1);
		stop
		else
		!write (acTemp, '("***********************************************************")');
		!call PrintToLog(acTemp, -1);
		!write (acTemp, '("*********** Reading cotrast matrix from file **************")');
		!call PrintToLog(acTemp, -1);
		!write (acTemp, '("***********************************************************")');
		!!call PrintToLog(acTemp, -1);
end if


END SUBROUTINE ContrastGenerationandFiltering

END MODULE ContrastGenFilt
