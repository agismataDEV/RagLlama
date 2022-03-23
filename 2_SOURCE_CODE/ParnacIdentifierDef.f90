MODULE ParnacIdentifierDef

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
!   The module ParnacIdentifierDef contains the definition of a number of
!   identifiers to be used for various goals within the program.
!
! *****************************************************************************
!
!   MODULES USED
!
    USE Types

! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   All variables are of type integer(i8b)
!
!   Possibilities for the Space Identifier:
!   iSI_PRIMARYSOURCE  - Contains the primary source slice (iDimZ=1 (or more, 
!                        depending on numprocs))
!   iSI_FIELD          - Contains a field space
!   iSI_CONTRASTSOURCE - Contains a contrast source space
!   iSI_FIELDSLICE     - Contains a single plane in the Z-dimension of the 
!                        field  (iDimZ=1 (or more, depending on numprocs))
!   iSI_XYZSLICE       - Contains an XYZ slice (iDimT=iProcN-1, dist=2)
!   iSI_INHOMCONTRAST  - Contains the inhomogeneity contrast space 
!                        (iDimT=iProcN-1, dist=0)
!
!   Possibilities for the Source excitation type cModelParams%PrimarySourceType:
!   iEI_QSOURCE   - use a volume source excitation, i.e. a velocity jump source
!   iEI_FSOURCE   - use a force source excitation, i.e. a pressure jump source
!   iEI_DTQSOURCE - use a volume source excitation but define it with the 
!                   temporal derivative of the volume source; this is convenient
!                   for measured signatures and it is equivalent to the FieldII 
!                   input excitation (apart from a factor 1/c0)
!   iEI_PLANEWAVE - use a plane wave as the incident field (so not as a source!)
!
!   Possibilities for the Contrast Identifier cModelParams%ContrastSourceType:
!   iCI_NONLIN          - Nonlinear contrast operator
!   iCI_GAUSSIAN        - Inhomogeneity contrast: Gaussian blob shape
!   iCI_LUNEBERG        - Inhomogeneity contrast: Luneberg lens
!   iCI_COMPLEXCONTRAST - Inhomogeneity contrast: phase screen, cylinder and blob
!   iCI_SPHERE          - Inhomogeneity contrast: sphere
!   iCI_LINCONTRA        - For the Linearized Contrast Source Method
!   iCI_CGCONTRA        - For the CG method
!
!   Possibilities for the Slicesave Identifier cModelParams%slicesavespecifier:
!   iAI_FIRSTANDLAST - Save slices of the first and the last iteration
!   iAI_LAST         - Save slices of the last iteration
!   iAI_ALL          - Save slices of all iterations
!
!   Possibilities for the Source signature Identifier cSourceParams%srcsigntype:
!   iGI_GAUSSIAN - Harmonic signal modulated with a Gaussian envelope
!   iGI_BLACKMAN - Blackman pulse
!   iGI_FILE     - Signature definition stored as data in a .dat or a .bin file
!
!   Possibilities for the Source shape Identifier cSourceParams%srcshapetype
!   iHI_CYLINDRICAL - Cylindrical aperture, with or without focusing
!   iHI_RECTANGULAR - Rectangular aperture
!   iHI_TRIANGULAR  - Triangular aperture
!   iHI_POINTSOURCE - Point source aperture
!   iHI_PHASEDARRAY - Phased array aperture
!   iHI_MATRIXARRAY - Matrix array aperture
!
! =============================================================================
!
    IMPLICIT NONE

	!identifiers for the Space Identifier
	integer(i8b), parameter :: iSI_PRIMARYSOURCE=1
	integer(i8b), parameter :: iSI_FIELD=2
	integer(i8b), parameter :: iSI_CONTRASTSOURCE=3
	integer(i8b), parameter :: iSI_FIELDSLICE=4
	integer(i8b), parameter :: iSI_XYZSLICE=5
	integer(i8b), parameter :: iSI_INHOMCONTRAST=6

	!identifiers for the Source excitation type
	integer(i8b), parameter :: iEI_QSOURCE=1
	integer(i8b), parameter :: iEI_FSOURCE=2
	integer(i8b), parameter :: iEI_DTQSOURCE=3
	integer(i8b), parameter :: iEI_PLANEWAVE=4
	integer(i8b), parameter :: iEI_POINTSOURCECLOUD=5
	integer(i8b), parameter :: iEI_LOADFIELD=6

	!identifiers for the Contrast Identifier
	integer(i8b), parameter :: iCI_NONLIN=1
	integer(i8b), parameter :: iCI_COMPLEXCONTRAST=2
	integer(i8b), parameter :: iCI_SPHERE=3
	integer(i8b), parameter :: iCI_LUNEBERG=4
	integer(i8b), parameter :: iCI_BLOB=5
	integer(i8b), parameter :: iCI_LINCONTRA=6
	integer(i8b), parameter :: iCI_CGCONTRA=7
	integer(i8b), parameter :: iCI_CGCONVERSION=8 
	integer(i8b), parameter :: iCI_CGTHETATWO=9   
	integer(i8b), parameter :: iCI_NONLINKAPPADX=10
	integer(i8b), parameter :: iCI_NONLINKAPPADY=11  
	integer(i8b), parameter :: iCI_NONLINKAPPADZ=12  
	integer(i8b), parameter :: iCI_NONLINKAPPATHETADX=13
	integer(i8b), parameter :: iCI_NONLINKAPPATHETADY=14
	integer(i8b), parameter :: iCI_NONLINKAPPATHETADZ=15
	integer(i8b), parameter :: iCI_NONLINKAPPADIRECTIONDX=16
	integer(i8b), parameter :: iCI_NONLINKAPPADIRECTIONDY=17
	integer(i8b), parameter :: iCI_NONLINKAPPADIRECTIONDZ=18
	integer(i8b), parameter :: iCI_SCATTERER=19
	    

	!identifiers for the Slicesave Identifier
	integer(i8b), parameter :: iAI_FIRSTANDLAST=1
	integer(i8b), parameter :: iAI_LAST=2
	integer(i8b), parameter :: iAI_ALL=3

	!identifiers for the Source signature Identifier
	integer(i8b), parameter :: iGI_GAUSSIAN=1
	integer(i8b), parameter :: iGI_BLACKMAN=2
	integer(i8b), parameter :: iGI_FILE=3

	!identifiers for the Source shape Identifier
	integer(i8b), parameter :: iHI_CYLINDRICAL=1
	integer(i8b), parameter :: iHI_RECTANGULAR=2
	integer(i8b), parameter :: iHI_TRIANGULAR=3
	integer(i8b), parameter :: iHI_POINTSOURCE=4
	integer(i8b), parameter :: iHI_PHASEDARRAY=5
!KH	integer(i8b), parameter :: iHI_MATRIXARRAY=6

END MODULE ParnacIdentifierDef
