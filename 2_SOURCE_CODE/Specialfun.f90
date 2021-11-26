
    MODULE SpecialFun

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen / Martin Verweij
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
    !   The module SpecialFun contains evaluations and faster approximations for
    !   a number of special functions used in the Parnac program.
    !
    ! *****************************************************************************
    !
    !   MODULES USED
    !
    USE Types
    USE Constants
    USE FFTW_cons
    include "fftw_wrappers.inc"

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
    !   cisi          sub   computes the Cosine and Sine Integral Ci(x) and
    !                       Si(x) for a real value x
    !   cisi4         sub   computes the Cosine and Sine Integral Ci(x) and
    !                       Si(x) for a real value x - much faster than cisi
    !   Expint_prime  sub   computes the Exponential Integral Ei(z)
    !                       multiplied with exp(z) for a complex value z
    !   dSinc         fun   computes the cardinal sine sin(x)/x for a real
    !                       value x
    !
    ! =============================================================================
    !
    CONTAINS

    SUBROUTINE cisi(x,ci,si)

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
    !   Source: W.H. Press, S.A. Teukolsky, W.T. Vetterling and B.P. Flannery,
    !           'Numerical Recipes in Fortran 77'. Cambridge: Cambridge University
    !           Press, 1996.
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine CiSi computes the cosine and sine integrals Ci(x) and Si(x).
    !   Ci(0) is returned as a large negative number and no error message is
    !   generated. For x<0 the routine returns Ci(-x) and you must supply the -i*pi
    !   yourself.
    !
    !   Method: see Press et al., Numerical Recipes, pp. 250vv. and 1125vv.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   x   i   dp    real value x at which Ci(x) and Si(x) are to be determined
    !   ci  o   dp    value of Ci(x)
    !   si  o   dp    value of Si(x)
    !
    REAL(dp), INTENT(in) :: x
    REAL(dp), INTENT(out) ::ci,si

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   MAXIT   i4b   maximum number of iterations allowed
    !   EPS     dp    relative error, or absolute error near a zero of Ci(x)
    !   FPMIN   dp    number near the smallest representable floating point number
    !   BIG     dp    number near the machine overflow limit
    !   TMIN    dp    dividing line between using the series and continued fraction
    !   i       i4b   loop counter in continued fraction
    !   k       i4b   loop counter in series expansion
    !   t       dp    abs(x)
    !   a       dp    temporary variable in continued fraction
    !   h       dpc   temporary variable in continued fraction
    !   b       dpc   temporary variable in continued fraction
    !   c       dpc   temporary variable in continued fraction
    !   d       dpc   temporary variable in continued fraction
    !   del     dpc   temporary variable in continued fraction
    !   fact    dp    factorial in series expansion
    !   term    dp    single term in series expansion
    !   sign    dp    sign of the term in series expansion
    !   err     dp    relative error of term in series expansion
    !   sum     dp    total sum for ci(x) and si(x) in series expansion
    !   sumc    dp    partial sum for ci(x) in series expansion
    !   sums    dp    partial sum for si(x) in series expansion
    !   odd     lgt   odd or even term in series expansion
    !
    INTEGER(i4b), PARAMETER :: MAXIT=100
    REAL(dp),PARAMETER :: EPS=1000*epsilon(x),FPMIN=4.0_dp*tiny(x),BIG=huge(x)*EPS,TMIN=2.0
    INTEGER(i4b) :: i,k
    REAL(dp) :: a,err,fact,sign,sum,sumc,sums,t,term
    COMPLEX(dpc) :: h,b,c,d,del
    LOGICAL(lgt) :: odd

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

    t=abs(x)
    if (t==0.0_dp) then
        si=0.0_dp
        ci=-BIG
        return
    end if

    if (t>TMIN) then
        !use continued fraction

        b=cmplx(1.0_dp,t,kind=dpc)
        c=BIG
        d=1.0_dp/b
        h=d
        do i=2,MAXIT
            a=-(i-1)**2
            b=b+2.0_dp
            d=1.0_dp/(a*d+b)
            c=b+a/c
            del=c*d
            h=h*del
            if ((abs(real(del)-1.0_dp)+abs(aimag(del)))<=EPS) exit
        end do
        if (i>MAXIT) stop 'cisi error: continued fraction failed'
        h=cmplx(cos(t),-sin(t),kind=dpc)*h
        ci=-real(h)
        si=half_pi+aimag(h)
    else
        if (t<sqrt(FPMIN)) then
            sumc=0.0_dp
            sums=t
        else
            !use series expansion
            sum=0.0_dp
            sums=0.0_dp
            sumc=0.0_dp
            sign=1.0_dp
            fact=1.0_dp
            odd=.true.
            do k=1,MAXIT
                fact=fact*t/k
                term=fact/k
                sum=sum+sign*term
                err=term/abs(sum)
                if (odd) then
                    sign=-sign
                    sums=sum
                    sum=sumc
                else
                    sumc=sum
                    sum=sums
                end if
                if (err < EPS) exit
                odd=.not. odd
            end do
            if (k>MAXIT) stop 'cisi error: MAXIT exceeded in cisi'
        end if
        si=sums
        ci=sumc+log(t)+EULER
    end if
    if (x <0.0_dp) si=-si

    END SUBROUTINE cisi

    subroutine cisi4(x,ci,si)

    ! =============================================================================
    !
    !   Programmer: M.D. Verweij
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   1.2     070101  Original code (MDV)
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
    !   The subroutine CiSi4 computes the cosine and sine integrals Ci(x) and Si(x).
    !   The routine is much faster than CiSi. It is only for x>0; for x<=0, it may
    !   result in an erroneous value or a runtime error.
    !   Relative error: 1e-3, 1e-5, 1e-7, depending on the number of terms
    !   included in the approximations.
    !
    !   Method: Taylor approximation or rational function approximation with
    !           the auxiliary functions f(x) and g(x). See M. Abramowitz and I.A.
    !           Stegun, 'Handbook of mathematical functions'. New York: Dover
    !           publications, 1968.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   x   i   dp    real positive value x at which Ci(x) and Si(x) are to be
    !                 determined
    !   ci  o   dp    value of Ci(x) with a relative error 1e-3, 1e-5 or 1e-7
    !   si  o   dp    value of Si(x) with a relative error 1e-3, 1e-5 or 1e-7
    !
    real(dp), intent(in) :: x
    real(dp), intent(out) :: si,ci

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   sq   dp    square of x
    !   f    dp    auxiliary function f(x)
    !   g    dp    auxiliary function g(x)
    !   s    dp    sin(x)
    !   c    dp    cos(x)
    !
    real(dp) :: sq,f,g,s,c

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

    !*** switch
    if (x <= 1e+0_dp) then
        !*** calculate Si and Ci for x<=1 using Taylor approximations
        sq=x*x
        !*** relative error 1e-3
        !03         si=((+1.666666666666667e-3_dp*sq&
        !03            -5.555555555555556e-2_dp)*sq &
        !03            +1e+0_dp)*x
        !03         ci=(+1.041666666666667e-2_dp*sq&
        !03            -2.5e-1_dp              )*sq&
        !03            +5.772156649015329e-1_dp+log(x)
        !*** relative error 1e-5
        !05         si=(((-2.834467120181406e-5_dp*sq&
        !05            +1.666666666666667e-3_dp)*sq  &
        !05            -5.555555555555556e-2_dp)*sq  &
        !05            +1e+0_dp)*x
        !05         ci=((-2.314814814814815e-4_dp*sq&
        !05            +1.041666666666667e-2_dp)*sq &
        !05            -2.5e-1_dp              )*sq &
        !05            +5.772156649015329e-1_dp+log(x)
        !*** relative error 1e-7
        si=((((+3.061924358220655e-7_dp*sq&
            -2.834467120181406e-5_dp)*sq   &
            +1.666666666666667e-3_dp)*sq   &
            -5.555555555555556e-2_dp)*sq   &
            +1e+0_dp)*x
        ci=(((+3.100198412698413e-6_dp*sq&
            -2.314814814814815e-4_dp)*sq  &
            +1.041666666666667e-2_dp)*sq  &
            -2.5e-1_dp              )*sq  &
            +5.772156649015329e-1_dp+log(x)
    else
        !*** calculate f and g for x>1 using rational function approximations
        !*** relative error 1e-3
        !03         f=((+2.784652088208877e-1_dp*x&
        !03           +1.584340895987836e+0_dp)*x &
        !03           +1.313137362705897e+0_dp)/  &
        !03           ((+2.786955097760614e-1_dp*x&
        !03           +1e+0_dp)*(x+1e+0_dp)**2)
        !03         g=(((+6.309775108591879e-1_dp*x&
        !03           +3.578080805383426e+0_dp)*x  &
        !03           +2.371756339729718e+0_dp)*x  &
        !03           +2.387834944980385e+0_dp)/   &
        !03           ((+6.311715544288667e-1_dp*x &
        !03           +1e+0_dp)*(x+1e+0_dp)**4)
        !*** relative error 1e-5
        !05         f=((((+6.099347983463999e-2_dp*x&
        !05           +6.514801175747663e-1_dp)*x   &
        !05           +2.475347340596883e+0_dp)*x   &
        !05           +2.993944684763584e+0_dp)*x   &
        !05           +1.422954929001066e+0_dp)/    &
        !05           (((+6.099311047241654e-2_dp*x &
        !05           +4.686499048628569e-1_dp)*x   &
        !05           +1e+0_dp)*(x+1e+0_dp)**3)
        !05         g=(((((+3.780264666622293e-2_dp*x&
        !05           +3.829201224596377e-1_dp)*x    &
        !05           +2.121730890478765e+0_dp)*x    &
        !05           +4.960992002047160e+0_dp)*x    &
        !05           +3.873071589020816e+0_dp)*x    &
        !05           +2.157558487954685e+0_dp)/     &
        !05           (((+3.780264176730856e-2_dp*x  &
        !05           +1.939109480682701e-1_dp)*x    &
        !05           +1e+0_dp)*(x+1e+0_dp)**5)
        !*** relative error 1e-7
        f=((((((+1.085051957486819e-1_dp*x&
            +1.668204866398108e+0_dp)*x     &
            +1.124991718241879e+1_dp)*x     &
            +3.560455310035718e+1_dp)*x     &
            +4.320834309430951e+1_dp)*x     &
            +2.094065499523842e+1_dp)*x     &
            +1.625127460448662e+0_dp)/      &
            (((((+1.085051974337423e-1_dp*x &
            +1.342687480828874e+0_dp)*x     &
            +7.113633189687764e+0_dp)*x     &
            +1.344695584637110e+1_dp)*x     &
            +1e+0_dp)*(x+1e+0_dp)**3)
        g=(((((((+2.128686646768257e-2_dp*x&
            +3.638317835431616e-1_dp)*x      &
            +3.004438448333865e+0_dp)*x      &
            +1.241230198911655e+1_dp)*x      &
            +2.164213978067445e+1_dp)*x      &
            +2.202987125965332e+1_dp)*x      &
            +9.863682327560993e+0_dp)*x      &
            +2.629409868517712e+0_dp)/       &
            (((((+2.128686670223448e-2_dp*x  &
            +2.573970478220391e-1_dp)*x      &
            +1.632402019056510e+0_dp)*x      &
            +3.638453979733940e+0_dp)*x      &
            +1e+0_dp)*(x+1e+0_dp)**5)
        !*** calculate sine and cosine for x>1
        s=sin(x)
        c=cos(x)
        !*** calculate values of sine and cosine integrals for x>1
        si=1.570796326794897e+0_dp-f*c-g*s
        ci=f*s-g*c
    end if
    !*** return
    end subroutine cisi4

    SUBROUTINE Expint_prime(z,E1p,err)

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
    !   The subroutine Expint_prime computes the alternative exponential integral
    !   E1(z) * exp(z).
    !
    !   Method: Series expansion or continued fraction. Generalized to a complex
    !           parameter from the methods described by W.H. Press, S.A. Teukolsky,
    !           W.T. Vetterling and B.P. Flannery, 'Numerical Recipes in Fortran
    !           77'. Cambridge: Cambridge University Press, 1996.See pp. 250vv. and
    !           1125vv.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   z    i   dpc   complex value z for which E1p(z) has to be determined.
    !   E1p  o   dpc   alternative exponential integral E1p(z)
    !   err  o   i8b   error number:
    !                   -2 for unsuccessful continued fraction
    !                   -1 for unsuccessful serie expansion
    !                    1 for successful serie expansion
    !                    2 for successful continued fraction
    !
    complex(dpc), INTENT(in) :: z
    complex(dpc), INTENT(out) :: E1p
    integer(i8b), INTENT(out) :: err

    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   MAXIT  i4b   maximum number of iterations allowed
    !   EPS    dp    relative error, or absolute error near a zero of E1p(x)
    !   FPMIN  dp    number near the smallest representable floating point number
    !   BIG    dp    number near the machine overflow limit
    !   i      i4b   loop counter
    !   selectSEregion   lgt   whether or not x falls in the region where the
    !                          series expansion is employed.
    !   rez    dp    real part of z
    !   imz    dp    imaginary part of z
    !   fact   dpc   factorial in the series expansion
    !   del    dpc   temporary variable in series expansion and continued fraction
    !   a      dp    temporary variable in continued fraction
    !   b      dpc   temporary variable in continued fraction
    !   c      dpc   temporary variable in continued fraction
    !   d      dpc   temporary variable in continued fraction
    !   h      dpc   temporary variable in continued fraction

    !   err     dp    relative error of term in series expansion


    INTEGER(i4b), PARAMETER :: MAXIT=600
    !	REAL(dp),PARAMETER :: EPS=epsilon(1D0)
    REAL(dp),PARAMETER :: EPS=1D-6
    REAL(dp),PARAMETER :: FPMIN=4.0_dp*tiny(1D0)
    REAL(dp),PARAMETER :: BIG=huge(1D0)
    INTEGER(i4b) :: i
    LOGICAL(lgt) :: selectSEregion
    COMPLEX(dpc) :: fact,del
    REAL(dp) :: a, rez, imz
    COMPLEX(dpc) :: b,c,d,h

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

    !first case: if |z|=0
    if (abs(z)==0.0_dp) then
        E1p=BIG
        err=0 ! err=0 denotes zero argument
        return
    end if

    !Select the region in which z is located
    !to determine the approximation to be used.
    ! The region consists of a triangle from (-18,0) to (-4,3.5) and (-4,-3.5)
    ! topped with an ellipse with center (-4.0) and axes (5.5,3.5)
    rez=real(z)
    imz=imag(z)

    selectSEregion=.false.
    if (rez.ge.-18.0_dp) then
        if (rez.le.-4.0_dp) then
            if (abs(imz).le.(3.5_dp/14.0_dp*rez+4.5_dp)) then
                selectSEregion=.true.
            endif
        elseif (rez.le.1.5_dp) then
            if (imz**2<(3.5_dp**2-(3.5_dp/5.5_dp*(rez+4.0_dp))**2)) then
                selectSEregion=.true.
            endif
        endif
    endif

    if (selectSEregion) then
        !Perform Series expansion
        err=-1 ! err=-1 denotes iteration overflow in Series expansion
        E1p = -EULER - log(z)
        fact = z
        del = z
        do i=2,MAXIT
            E1p = E1p + del
            fact = -z*fact/i
            del = fact/i
            if(abs(del).lt.EPS) then
                err=1 ! err=1 denotes successfull Series expansion
                exit
            endif
        enddo
        !		err=i ! override err: err is now number of iterations
        E1p=(E1p + del)*exp(z)

    else
        !Perform Continued fraction
        err=-2 ! err=-2 denotes iteration overflow in Continued fraction
        b=z+1.0_dp
        c=BIG
        d=1.0_dp/b
        h=d
        do i=1,MAXIT
            a=-i**2
            b=b+2.0_dp
            d=1.0_dp/(a*d+b)
            c=b+a/c
            del=c*d
            h=h*del
            if (abs(del-1.0_dp).lt.EPS) then
                err=2 ! err=2 denotes successfull Continued fraction
                exit
            end if
        end do
        !		err=i ! override err: err is now number of iterations

        !		If ((real(z).lt.0).and.(imag(z)==0)) then
        !			E1p = h - im*pi
        !		else
        E1p = h
        !		endif

    endif

    END SUBROUTINE Expint_prime

    FUNCTION dSinc(dArg)

    ! =============================================================================
    !
    !   Programmer: Koos Huijssen/ Agisilaos Matalliotakis (200614)
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
    !   The function dSinc computes the cardinal sine sin(x)/x of an input argument
    !   x. If x = 0, dSinc(x) = 1.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   dArg    i   dp   value for which dSinc(dArg) has to be determined.
    !   dSinc   o   dp   result
    !
    real(dp), intent(in) :: dArg(:)
    real(dp) :: dSinc(size(dArg))

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

    ! Logical acts give -1 if TRUE and 0 if FALSE
    dSinc = 0.0D0
    ! If the argument is 0 then change it to 1D-30/1D-30 in order to avoid 0/0 = 1
    dSinc = dsin(dArg+EPSILON(1.0D0))/real(dArg+EPSILON(1.0D0),dp)
   
    END FUNCTION dSinc

    SUBROUTINE LINSPACE(x,x_start,x_end,x_len)
	! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis // Date : 200902
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   3.0     200902  Original code (KAM)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine LINSPACE gemerates a vector x, when the starting value (x_start),
	!	the final value (x_end) and the length of interest (x_len) are known.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   Nominator             	i   i8b  		Nominator, Number we are looking of its divisors
    !   Denominator			  	i   i8b         Denominator, the divisor
    !   Result            		io  i8b         Closest Integer Divisor 
    !
	! *****************************************************************************
    implicit none
    real(dp), dimension(:), intent(out) :: x
    real(dp), intent(in)    			:: x_start, x_end
    integer(i8b),intent(in)  			:: x_len
	! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    ! *****************************************************************************
	real(dp)   							:: dx
    integer(i8b)		 				:: i
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
	
    dx = (x_end-x_start)*1.0D0/(x_len-1)

    x = [(x_start + dx*(i-1), i = 1,x_len)]
    END SUBROUTINE LINSPACE

    SUBROUTINE INTERP1D( xData, yData, xVal, yVal )
	! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis // Date : 200902
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   3.0     200902  Original code (KAM)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine FIND_CLOSEST_DIVISOR finds an integer divisor of a given
    !	number (Nominator) which is closest to a given second number (Denominator). 
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   xData             	i   i8b  		a vector of the x-values of the data to be interpolated
    !   yData			  	i   i8b         a vector of the y-values of the data to be interpolated
    !   xVal            	i   i8b         a vector of the x-values where interpolation should be performed
    !   yVal            	o   i8b         a vector of the resulting interpolated values
    !
    ! *****************************************************************************
    real(dp), intent(in) 	:: 	xData(:), yData(:), xVal(:)
    real(dp), intent(out) 	:: 	yVal(:)

	! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    ! *****************************************************************************
    integer(i8b)		    ::  inputIndex, i_start, i_end, i_min, i_max
	integer(i8b) 			::	xval_start, xVal_end, max_xData_pos, max_xVal_pos
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
    max_xData_pos = size(xData,1)
    max_xVal_pos  = size(xVal,1)
    
    i_min = minloc(xVal,1 , xVal > xData(1))
    xval_start = 1
        
    i_max = maxloc(xVal, 1, xVal< xData(max_xData_pos))+1
    xval_end = size(xVal,1)
    
    ! Check if there are values in the xVal higher than the max value of xData.
    ! If true, this means that the interpolation should take place using the last 2 values of the xData
    ! if there are more than 1 value, then do this for all the data and change the iteration range
    ! else  change on
    if (size(xVal,1)==1) then
        xval_end = 1
    elseif (i_max /= size(xVal,1)) then
        yVal(i_max : max_xVal_pos) = (yData(max_xData_pos)-yData(1))*1.0_dp/max_xData_pos*(xVal(i_max : max_xVal_pos)-xData(max_xData_pos))+yData(max_xData_pos)
        xval_end = i_max-1
    else
        yVal(i_max) = (yData(max_xData_pos)-yData(1))*1.0_dp/max_xData_pos*(xVal(i_max)-xData(max_xData_pos))+yData(max_xData_pos)
        xval_end = i_max - 1
    endif
    
    if (i_min /= 1) then
       xval_start=i_min+1
    endif
    
    do inputIndex = xval_start,  xval_end
        i_start = minloc(xData, 1, xData > xVal(inputIndex))-1
        i_end = minloc(xData, 1, xData > xVal(inputIndex))
 
        yVal(inputIndex) = (yData(i_end)-yData(i_start))*1.0_dp/(xData(i_end)-xData(i_start))*(xVal(inputIndex)-xData(i_start))+yData(i_start)
    end do
   
    END SUBROUTINE INTERP1D
	

    SUBROUTINE INTERP3D(LOCDATA,FDATA,LOCFINAL,FFINAL)
	! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis // Date : 200902
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   3.0     200902  Original code (KAM)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine FIND_CLOSEST_DIVISOR finds an integer divisor of a given
    !	number (Nominator) which is closest to a given second number (Denominator). 
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   xData             	i   i8b  		a vector of the x-values of the data to be interpolated
    !   yData			  	i   i8b         a vector of the y-values of the data to be interpolated
    !   xVal            	i   i8b         a vector of the x-values where interpolation should be performed
    !   yVal            	o   i8b         a vector of the resulting interpolated values
    !
    ! *****************************************************************************
    real(dp), intent(in) 		:: 		LOCDATA(:,:), FDATA(:,:), LOCFINAL(:)
    real(dp), intent(out) 		:: 		FFINAL(:)

	! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    ! *****************************************************************************
    integer(i8b)		 		:: 		inputIndex, dataIndex, i_start,i_end
    real(dp)    				::		x_d,xmin,xmax,y_d,ymin,ymax,z_d,zmin,zmax
    real(dp)    				:: 		c00(size(FFINAL,1)),c01(size(FFINAL,1)),c10(size(FFINAL,1)),c11(size(FFINAL,1)),c0(size(FFINAL,1)),c1(size(FFINAL,1))
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
  

    ! xmin is the minimum location of the point in x axis , found as the one smaller than the asked x location  
    ! xmax is the maximum location of the point in x axis , found as the one greater than the asked x location
    ! Same idea for the rest values
    xmin = LOCDATA(minloc(LOCDATA(:,1),1) ,1)
    xmax = LOCDATA(maxloc(LOCDATA(:,1),1) ,1)
    x_d  = (LOCFINAL(1) - xmin)/(xmax - xmin)

    ymin = LOCDATA(minloc(LOCDATA(:,2),1) ,2)
    ymax = LOCDATA(maxloc(LOCDATA(:,2),1) ,2)
    y_d  = (LOCFINAL(2) - ymin)/(ymax - ymin)

    zmin = LOCDATA(minloc(LOCDATA(:,3),1) ,3)
    zmax = LOCDATA(maxloc(LOCDATA(:,3),1) ,3)
    z_d  = (LOCFINAL(3) - zmin)/(zmax - zmin)
    
    c00=  FDATA(1,:)*(1-x_d) + FDATA(2,:)*x_d
    c01 = FDATA(3,:)*(1-x_d) + FDATA(4,:)*x_d
    c10 = FDATA(5,:)*(1-x_d) + FDATA(6,:)*x_d
    c11 = FDATA(7,:)*(1-x_d) + FDATA(8,:)*x_d

    c0 = c00*(1-y_d)+c10*y_d
    c1 = c01*(1-y_d)+c11*y_d

    FFINAL = c0*(1-z_d)+c1*z_d

    END SUBROUTINE INTERP3D
    
    SUBROUTINE FIND_CLOSEST_DIVISOR(Nominator, Denominator, Result)
    
    ! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis // Date : 200902
    !
    !   Language: Fortran 90
    !
    !   Version Date    Comment
    !   ------- -----   -------
    !   3.0     200902  Original code (KAM)
    !
    ! *****************************************************************************
    !
    !   DESCRIPTION
    !
    !   The subroutine FIND_CLOSEST_DIVISOR finds an integer divisor of a given
    !	number (Nominator) which is closest to a given second number (Denominator). 
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   Nominator             	i   i8b  		Nominator, Number we are looking of its divisors
    !   Denominator			  	i   i8b         Denominator, the divisor
    !   Result            		io  i8b         Closest Integer Divisor 
    !
    integer(i8b), intent(in)::          		Nominator, Denominator
    integer(i8b), intent(inout)::                      	Result
    ! *****************************************************************************
    !
    !   LOCAL PARAMETERS
    !
    !   i            i8b    Loop counter
    
    integer(i8b) :: i
    
    if (Denominator <=1) then
    	Result = 1 
    	return
    else if (mod(Nominator,Denominator)==0) then
    	Result = Denominator
    	return
    endif
    i=1
    do while(i<Denominator)
    	if (mod(Nominator,i)==0) Result = i
		i=i+1
    end do
    return
  END SUBROUTINE FIND_CLOSEST_DIVISOR
  
  FUNCTION Solve_Geometric_Series(a,n,S)
	  real(dp), intent(in) 		:: a
	  integer(i8b),intent(in) 	:: n,S
	  real(dp)					:: Solve_Geometric_Series
	  integer(i8b) 				:: i, Calc_S
	  real(dp)					:: r, startD, endD
	  Solve_Geometric_Series =1.0D0
	  if (S*1.0D0/n == a) return	  
	  i=1
	  Calc_S = 0
	  do while(Calc_S<1)
	  i=i+1
	  Calc_S = INT(a*(i**n-1)/(i-1)/S,8)
	  enddo
	  
	  r = real(i,dp)
	  startD = r-1.0D0
	  endD = r
	  do while(abs(S-Calc_S) > a*1E2)
		  r = (endD+startD)/2
		  Calc_S = INT(a*(r**n-1)/(r-1),8)		  
		  	endD = r*abs(Calc_S>S) + (1-abs(Calc_S>S))*endD
		  	startD = startD*abs(Calc_S>S) + (1-abs(Calc_S>S))*r
	  enddo
	  Solve_Geometric_Series = r
	  
	  return
  
  END FUNCTION Solve_Geometric_Series

END MODULE SpecialFun
