MODULE ParnacParamDef

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
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The module ParnacParamDef contains global declarations of structures that
!   contain all kinds of parameters for the program, many of which are input
!   through the configuration file. Most other parameters are derived from the
!   input parameters.
!
! *****************************************************************************
!
!   MODULES USED
!
	USE Types
!
! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   MediumInput    type   Definition of structure containing medium parameters
!   ModelInput     type   Definition of structure containing model parameters
!   SourceInput    type   Definition of structure containing source parameters
!   cMediumParams  struc  Structure containing medium parameters                
!   cModelParams   struc  Structure containing model parameters                
!   cSourceParams  struc  Structure containing source parameters                
!
! =============================================================================
!
	IMPLICIT NONE
	SAVE

	!*************************************************************************************
	!
	!	Structure containing parameters with physical information
	!
	!*************************************************************************************
	
	type MediumInput

		!medium constants
		real(dp) :: rho0          !Ambient density [kg m^{-3}]
		real(dp) :: c0            !Small-signal sound velocity [ m s^{-1}]
		real(dp) :: fc0           !Frequency at which sound velocity is given [Hz] 
    		real(dp) :: mu		  !Medium viscosity is given [Pa*sec]
    		real(dp) :: P0		  !Ambient Pressure [Pa]
		real(dp) :: a_att         !Frequency-dependent attenuation parameter [Np/cm MHz^b]
		real(dp) :: b_att         !Power in frequency power attenuation law []
		real(dp) :: alpha0_att	  !Derived frequency-dependent medium parameter
		real(dp) :: kappa0        !Ambient Compressibility
		real(dp) :: beta          !Nonlinearity parameter []
		real(dp) :: Z0            !Acoustic impedance

	end type MediumInput

	type(MediumInput):: cMediumParams

	!*************************************************************************************
	!
	!	Structure containing parameters for the algorithm
	!
	!*************************************************************************************
	
	type ModelInput

		!run parameters
		real(dp) :: Lx, Ly, Lt			! Width of interest [mm]
		real(dp) :: Sz, Lz				! Width of interest [mm]
		real(dp) :: Thetax, Thetay, Thetat	! The angles of the different planes [rad]
		logical(lgt) :: UseYSymmetry	   !Use Y symmetry in the beams (determined in specialchecks subroutine)
		
		real(dp) :: Xmin,Xmax    !Width of interest [mm]
		real(dp) :: Ymin,Ymax    !Height of interest  [mm]
		real(dp) :: Zmax		 !Distance of interest  [mm]
		real(dp) :: freq0        !Fundamental frequency [Hz]
		real(dp) :: Fnyq         !Maximum frequency (Nyquist) [factor relative to freq0]

		!FD parameters
		integer(i8b) :: FDTorder		!Max. order of finite difference schemes (even!)
		integer(i8b) :: FDXorder		!Max. order of finite difference schemes (even!)	

		!Parameters for the primary and contrast source operators
		integer(i8b) :: PrimarySourceType  !Use (volume) Q-source, (force) F-source, source from measured signature (sort of dQ/dt) or Plane wave
		integer(i8b) :: ContrastSourceType !Determine contrast source type	
		logical(lgt) :: UseAntiAliasing	   !Employ an anti-aliasing procedure in the contrast sources
		logical(lgt) :: UseSupportTapering !Taper the support at beginning and end 
		logical(lgt) :: UseFreqTapering !Taper the frequency at the end
		character(LEN=128) :: FieldFilename
		real(dp) :: xyzfielddim(3,2)	!Which dimension: x, y, z or t

		!Beams, iterations, save slices
		integer(i8b) :: numbeams						!Number of sub-beams in which to split the z-dimension
		integer(i8b),allocatable :: numiterations(:)	!Number of iterations for each beam
		integer(i8b) :: iIter							!Number of iterations for each beam
		integer(i8b) :: numslices						!Number of slices
		integer(i8b) :: slicesavespecifier				!Specifies what to save: first and last iterate (0), last only (1), or all iterates (2)
		character(LEN=1), allocatable :: xyzslicedim(:)	!Which dimension: x, y, z or t
		real(dp), allocatable :: xyzslicepos(:)			!Position in [mm]
		integer(i8b), allocatable :: xyzslicebeam(:)	!Which beam the slice is located; -1 if it may be all
		integer(i8b), allocatable :: xyzsliceindex(:)	!Index within the beam
		
		
		integer(i8b) :: PPW			! Points Per Wavelength, Default = 2
	end type ModelInput

	type(ModelInput) :: cModelParams

	!*************************************************************************************
	!
	!	Structure containing the source excitation parameters
	!
	!*************************************************************************************
	
	type SourceInput

		!source filename
		character(LEN=128) :: srcfilename !name of TXY source plane file

		!source signature params
		integer(i8b) :: srcsigntype !Type of source signature, or name of signature file
		real(dp) :: Pstart		  !Source pressure amplitude
		real(dp) :: Tpulse        !Total pulse length
		real(dp) :: Tdelay        !Delay of gaussian window or blackman pulse

		!source signature: gaussian pulse
		real(dp) :: Twidth        !Width of gaussian window 
		integer(i8b) :: power	  !Power of gaussian window

		!source signature: file
		character(LEN=128) :: srcsignfilename

		!source shape params
		integer(i8b) :: srcshapetype !Type of source shape

		!flat cylindrical source
		real(dp) :: radius		  !Source radius 
		real(dp) :: radfocus	  !Source focus (in z direction)

		!flat rectangular source
		real(dp) :: rectwidth	  !Width (in x direction)  
		real(dp) :: rectheight	  !Height (in y direction

		!flat triangular source
		real(dp) :: triangedge	  !Edge size
		real(dp) :: triangxcenter !Center in the x dimension

		!phased array source
		integer(i8b) :: numel	  !Number of phased array elements
		real(dp) :: elwidth		  !Width of each element (x)
		real(dp) :: elheight	  !Height of each element (y)
		real(dp) :: kerf		  !Distance between elements (x)
		real(dp) :: focusx		  !x-coordinate of focus
		real(dp) :: focusz		  !z-coordinate of focus
		real(dp) :: elevationfocusz !z-coordinate of focus in elevation (y) - direction ==0 if no el. focus
		character(LEN=128) :: phaseapodfilename 	!File in which element phases and apodizations are stored
		character(LEN=128) :: td_filename		 	!File in which element phases and apodizations are stored

!KH     !matrix array source
!KH     integer(i8b) :: matnumelx
!KH     integer(i8b) :: matnumely
!KH     real(dp) :: matelwidth
!KH     real(dp) :: matelheight
!KH     real(dp) :: matkerfx
!KH     real(dp) :: matkerfy
!KH     real(dp) :: matfocusx
!KH     real(dp) :: matfocusy
!KH     real(dp) :: matfocusz

        integer(i8b) :: NKC
        integer(i8b) :: CG
        integer(i8b) :: LCSM
        integer(i8b) :: BiCGSTAB
        integer(i8b) :: LINEARCASEIMPLEMENTATION
        integer(i8b) :: DFM
        integer(i8b) :: SHIFT
        
        real(i8b) :: alpha1
        real(i8b) :: b1
        real(i8b) :: SOS1
        real(i8b) :: KAPPA1
        real(i8b) :: BETA1
        
        real(i8b) :: alpha2
        real(i8b) :: b2
        real(i8b) :: SOS2
        real(i8b) :: KAPPA2
        real(i8b) :: BETA2

        real(i8b) :: alpha3
        real(i8b) :: b3
        real(i8b) :: SOS3
        real(i8b) :: KAPPA3
        real(i8b) :: BETA3
        
        real(i8b) :: alpha4
        real(i8b) :: b4
        real(i8b) :: SOS4
        real(i8b) :: KAPPA4
        real(i8b) :: BETA4
        
        real(i8b) :: alpha5
        real(i8b) :: b5
        real(i8b) :: SOS5
        real(i8b) :: KAPPA5
        real(i8b) :: BETA5
        
        !
	character(LEN=256) :: residualdirectory
        
	real(dp) :: arraywidth,maxphasedelay !Phase delay due to focusing (derived in normalize_constants)

	end type SourceInput

	type(SourceInput) :: cSourceParams
    
    type BubbleInput

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
    !   The subroutine remove_CR replaces any Carriage Return symbols (ASCII code
    !   13) by a space. This is necessary to prevent a mismatch between Fortran
    !   and C (or Fortran and Windows) line ending types.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   string   io   char   string to be processed
    !
    

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   lenstring   i4b   length of the string
    !   i           i4b   temporary variable
    
    character(len = 1024)	::  sBubbleDir 
    character(len = 1024)	::  GridPointsPressure 
    !--------------Cluster Parameters--------------
    integer(i8b)   		  ::  	N
    real(dp)  			  ::  	ClusterDimsRatio(2,3)
    character(LEN=128)	  ::	Distribution  
    real(dp)  			  ::  	MinInBetweenDist
    integer(i8b)		  :: 	ClusterSlicesN
    character(LEN=1), allocatable :: 	ClusterSliceDim(:)
    integer(i8b),allocatable  :: OutputIndexes(:)
    
    !------------------ Bubble parameters
    
    real(dp),allocatable  ::  R0(:) 
    real(dp)  ::  PDRange(2)
    real(dp)  ::  kappa_s

    !------------------ Medium parameters (water, Room temperature =20° and 1 atm ambient pressure)
    real(dp)  ::  sigma_w 
    real(dp)  ::  sigma_R0
    real(dp)  ::  sigma_R
    real(dp)  ::  gama
    real(dp)  ::  P_g0

    !---------------Marmottant model parameters
    real(dp)  ::  chi 
    real(dp)  ::  R_b
    real(dp)  ::  R_r
	real(dp)	:: time_norm
	real(dp)	:: rad_norm
    real(dp),allocatable  :: P_driv(:), T_driv(:)
	character(LEN=2) :: Solver_Method
	character(LEN=24) :: Solver_Normalize
   
    !------------- Experimental
    real(qp)  :: coeff_fit(11)
    real(qp)  :: A0c
    
    end type BubbleInput
    
    type(BubbleInput) :: BubbleParams
	
	type PointSourceCloudInput

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
    !   The subroutine remove_CR replaces any Carriage Return symbols (ASCII code
    !   13) by a space. This is necessary to prevent a mismatch between Fortran
    !   and C (or Fortran and Windows) line ending types.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   string   io   char   string to be processed
    !
    

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !
    !   lenstring   i4b   length of the string
    !   i           i4b   temporary variable
    
    !--------------Proton Cluster Parameters--------------
    integer(i8b)   		  ::  	N
    real(dp)  			  ::  	ClusterDimsRatio(2,3)
    real(dp)  			  ::  	MinInBetweenDist
    real(dp)  			  ::  	PointSourceAmplitude 
    real(dp)  			  ::  	Dist_Amplitude   
    real(dp),allocatable  ::    R0(:) 

    end type PointSourceCloudInput
    
    type(PointSourceCloudInput) :: PointSourceCloudParams


END MODULE ParnacParamDef
