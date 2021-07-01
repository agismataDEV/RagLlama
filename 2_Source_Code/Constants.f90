MODULE Constants

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
!   The module Constants contains definitions of a number of often used
!   constants.
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
!   pi      dp  gives the number pi
!   two_pi  dp  gives 2*pi
!   half_pi dp  gives pi/2
!   euler   dp  gives Euler's constant \gamma
!   e       dp  gives the real number e = exp(1)
!   im      dpc gives the imaginary number sqrt(-1)
!
! =============================================================================

    IMPLICIT NONE

    SAVE

    real(dp), parameter ::                                                              &
        pi      = 3.141592653589793238462643383279502884197_dp,                         &
        two_pi  = 6.283185307179586476925286766559005768394_dp,                         &
        half_pi = 1.570796326794896619231321691639751442099_dp,							&
		euler   = 0.5772156649015328606065120900824024310422_dp,						&
		e		= 2.7182818284590452353602874713526624977572_dp		

     complex(dpc), parameter	::                                                          &
         im      = cmplx(0.0_dp,1.0_dp, dp)

END MODULE Constants
