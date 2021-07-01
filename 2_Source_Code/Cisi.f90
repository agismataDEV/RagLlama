
MODULE NR_Cisi
	USE Types
	USE Constants


IMPLICIT NONE
CONTAINS

SUBROUTINE cisi(x,ci,si)
	REAL(dp), INTENT(in) :: x
	REAL(dp), INTENT(out) ::ci,si
	INTEGER(i4b), PARAMETER :: MAXIT=100
	REAL(dp),PARAMETER :: EPS=1000*epsilon(x),FPMIN=4.0_dp*tiny(x),BIG=huge(x)*EPS,TMIN=2.0

	! Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is returned as a large
	! negative number and no error message is generated. For x<0 the routine returns Ci(-x)
	! and you must supply the -i*pi yourself.
	! Parameters: MAXIT is the maximum number of iterations allowed; EPS is the relative error,
	! or absolute error near a zero of Ci(x); FPMIN is a number near the smallest representable
	! floating point number; BIG is a number naer the machine overflow limit; TMIN is the 
	! dividing line between using the series and continued fraction; EULER=gamma (in constants)
	! see Numerical Recipes, pp. 250vv. and 1125vv.

	INTEGER(i4b) :: i,k
	REAL(dp) :: a,err,fact,sign,sum,sumc,sums,t,term
	COMPLEX(dpc) :: h,b,c,d,del
	LOGICAL(lgt) :: odd

	t=abs(x)
	if (t==0.0_dp) then
		si=0.0_dp
		ci=-BIG
		return
	end if

	if (t>TMIN) then
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

	!*** CiSi4
	!***
	!*** Fast evaluation of Cosine and Sine integrals.
	!*** Relative error: 1e-3, 1e-5, 1e-7.
	!*** Author: M.D. Verweij - Laboratory of Electromagnetic Research -
	!***         Delft University of Technology.
	!*** Version 1.2 - January 3, 2007.
	!***

	!*** declarations
	use types
	implicit none
	real(dp), intent(in) :: x
	real(dp), intent(out) :: si,ci
	real(dp) :: sq,f,g,s,c

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
!03		         f=((+2.784652088208877e-1_dp*x&
!03		           +1.584340895987836e+0_dp)*x &
!03		           +1.313137362705897e+0_dp)/  &
!03		           ((+2.786955097760614e-1_dp*x& 
!03		           +1e+0_dp)*(x+1e+0_dp)**2)
!03		         g=(((+6.309775108591879e-1_dp*x&
!03		           +3.578080805383426e+0_dp)*x  &
!03		           +2.371756339729718e+0_dp)*x  &
!03		           +2.387834944980385e+0_dp)/   &
!03		           ((+6.311715544288667e-1_dp*x &
!03		           +1e+0_dp)*(x+1e+0_dp)**4)
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

END MODULE NR_Cisi
