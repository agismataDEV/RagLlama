    MODULE RP
	USE Types
	USE Constants
	USE ParnacGeneral
	USE ParnacDataDef    
	USE ParnacParamDef
	USE ParnacIdentifierDef  
	USE SpecialFun
    USE ParnacParamInit, ONLY : &
        ScattererInit
    USE Constants  
	USE ParnacTransformFilter, ONLY : &
       dTaperingWindow, INTERP1DFREQ
	!   *****************************************************************************
	!
	!   GLOBAL DECLARATIONS
	!
	!   none
	!
	IMPLICIT NONE

    CONTAINS 

    SUBROUTINE RP_SOLVER(RealPressIn, RealTimeIn, n_samples, iBubble, V_dd_pad, RealTimeOut)
    
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
    !   The subroutine RP_SOLVER_ODEPACK solves the Rayleigh - Plesset equation
    !   using the ODEPACK solver (ODEPACK.f90).
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   cSpace             io   type(space)  Space from which the data array needs
    !                                        to be stored.
    !   MinDistance  i    i8b          Minimum distance between microbubbles
    !
	
    integer(i8b), intent(in) ::  iBubble, n_samples
    real(dp), intent(in)     ::  RealPressIn(n_samples), RealTimeIn(n_samples)
    real(dp), intent(inout)    ::  V_dd_pad(n_samples), RealTimeOut(n_samples)
	
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
	
    integer(i8b), parameter ::  NEQ_ODEPACK(1) = 2 , NEQ_RK = 5		!   Number of first order ode's
    integer(i4b), parameter ::  lrw = 697			!   LRW   :IN     Declared length of RWORK (in user's DIMENSION statement).
    integer(i4b), parameter ::  liw = 45            !   LIW   :IN     Declared length of IWORK (in user's DIMENSION statement).
    real(dp)    			::  atol(NEQ_ODEPACK(1)), rtol, dtout, t, tout, t_ref, y(NEQ_RK), yp(NEQ_RK) 
	real(dp)				::	TOLSF, ATOL_UP(NEQ_ODEPACK(1)), RTOL_UP
    real(dp)				::  rwork(lrw)
    integer(i4b) 			::  iwork(liw)
    integer(i4b) 			::  i, iopar, iopt, iout, istate, itask, itol, mf, jt, flag
	real(dp) 				::  norm_factor(3), P_interp(1)
	real(dp)				::  R_bub(n_samples,3), dTaperSupportWindowN(n_samples), P_bub(n_samples)
 
	character(len = 1024)	::  actemp  
	
	logical					::  log_isnan
    
    call ScattererInit()
    call R_exp(iBubble)  
      
	! If ScattererParams(cSourceParams%iScCloud)%rad_norm = ScattererParams(cSourceParams%iScCloud)%R0 , The result of the RP Solver is normalized to x = R/ScattererParams(cSourceParams%iScCloud)%rad_norm = R/R0 , so the solved value is multiplied by ScattererParams(cSourceParams%iScCloud)%R0 at the end of the code
    ! atol and rtol shoudld be 1E-5
    ! If ScattererParams(cSourceParams%iScCloud)%rad_norm =  1.0D0, the result of the RP solver ir R, so the solved value should not be multiplied with any factor
    ! atol and rtol shoudld be 1E-9. DLSODE input should change to MARMOTTANT_R

	if (trim(ScattererParams(cSourceParams%iScCloud)%Solver_Normalize) =='time') ScattererParams(cSourceParams%iScCloud)%time_norm = cModelParams%freq0
	if (trim(ScattererParams(cSourceParams%iScCloud)%Solver_Normalize) =='radius') ScattererParams(cSourceParams%iScCloud)%rad_norm  = ScattererParams(cSourceParams%iScCloud)%R0(iBubble)
	if (trim(ScattererParams(cSourceParams%iScCloud)%Solver_Normalize) =='freq0_and_radius') then
		 ScattererParams(cSourceParams%iScCloud)%time_norm = cModelParams%freq0
		 ScattererParams(cSourceParams%iScCloud)%rad_norm  = ScattererParams(cSourceParams%iScCloud)%R0(iBubble)
	elseif (trim(ScattererParams(cSourceParams%iScCloud)%Solver_Normalize) =='Minnaert_and_radius') then
		 ScattererParams(cSourceParams%iScCloud)%time_norm = SQRT(cMediumParams%P0/cMediumParams%rho0)/ScattererParams(cSourceParams%iScCloud)%R0(iBubble)
		 ScattererParams(cSourceParams%iScCloud)%rad_norm  = ScattererParams(cSourceParams%iScCloud)%R0(iBubble)
	endif 

    itol = 2            !  if atol scalar, itol = 1 and  if atol array, itol = 2, if atol and rtol array, itol = 4
    rtol = 10.0**(floor(log10(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)/ScattererParams(cSourceParams%iScCloud)%rad_norm))-8.0D0)
    atol = 10.0**(floor(log10(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)/ScattererParams(cSourceParams%iScCloud)%rad_norm))-8.0D0)
	 
	dTaperSupportWindowN = dTaperingWindow(n_samples,(RealTimeIn(2)-RealTimeIn(1))* cModelParams%freq0,2.0_dp,2.0_dp)

    ALLOCATE(ScattererParams(cSourceParams%iScCloud)%P_driv(n_samples ),ScattererParams(cSourceParams%iScCloud)%T_driv(n_samples ))
	ScattererParams(cSourceParams%iScCloud)%P_driv = RealPressIn * dTaperSupportWindowN
	ScattererParams(cSourceParams%iScCloud)%T_driv = RealTimeIn * ScattererParams(cSourceParams%iScCloud)%time_norm
	
	! Medium parameters (water, Room temperature =20° and 1 atm ambient pressure) 
    ScattererParams(cSourceParams%iScCloud)%P_g0 = cMediumParams%P0+2.0D0*ScattererParams(cSourceParams%iScCloud)%sigma_R0/ScattererParams(cSourceParams%iScCloud)%R0(iBubble)

    ! Marmottant model parameters
    ScattererParams(cSourceParams%iScCloud)%R_b =  ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0D0/sqrt(ScattererParams(cSourceParams%iScCloud)%sigma_R0*1.0D0/ScattererParams(cSourceParams%iScCloud)%chi+1.0D0)
    ScattererParams(cSourceParams%iScCloud)%R_r =  ScattererParams(cSourceParams%iScCloud)%R_b*sqrt(ScattererParams(cSourceParams%iScCloud)%sigma_w*1.0D0/ScattererParams(cSourceParams%iScCloud)%chi+1.0D0)  
     
    y(1) = ScattererParams(cSourceParams%iScCloud)%R0(iBubble)/ScattererParams(cSourceParams%iScCloud)%rad_norm - 1.0D0
    y(2) = 0.0D0 ; y(3) = 0.0D0 ; y(4) = iBubble*1.0D0; 
	R_bub = 0.0D0 ;     R_bub(1,:) = y(1:3) 
    
    tout 	=	0.0D0
    t 		= 	0.0D0
	RealTimeOut  = ScattererParams(cSourceParams%iScCloud)%T_driv
	
	if (ScattererParams(cSourceParams%iScCloud)%Solver_Method=='RK') then !For the case of Runge Kutta Solver
		flag = 1 			!  This is for RK
		
		do iout = 2,n_samples
			tout = ScattererParams(cSourceParams%iScCloud)%T_driv(iout) ; 
			call R8_RKF45(EXP_SIGMA, NEQ_RK, y, yp, t, tout, rtol, atol, flag )
			R_bub(iout,:) = y(1:3) 
			RealTimeOut(iout) = t 
		enddo
		
	else 
		!For the case of ODEPACK solver	
		mf = 22             !   memory flag : 10 (Adams) Stiff - 22 (BDF) Non-stiff code
		jt  = 2 			!   Jacobian Type : 2 - For DLSODA instead of mf
		iopt = 1            !   0 for not optional inputs, 1 for optional inputs
		itask = 1
		istate = 1		
		CALL XSETF(0) ! 0 if  no messages should be printed in ODEPACK call, comment for normal messaging				y(5) = ScattererParams(cSourceParams%iScCloud)%P_driv(iout)
    	do iout = 2,n_samples -1 
		
				ATOL_UP = ATOL
				RTOL_UP = RTOL
				tout = ScattererParams(cSourceParams%iScCloud)%T_driv(iout) !tout + (RealTimeIn(2)-RealTimeIn(1)) * ScattererParams(cSourceParams%iScCloud)%time_norm ; 
				y(5) = ScattererParams(cSourceParams%iScCloud)%P_driv(iout)
				call dlsode(MARMOTTANT,NEQ_ODEPACK,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,JACDUM,mf)

				do while (ISTATE <0)
				
					if (ISTATE == -2) then
						TOLSF = RWORK(14) ; if (isnan(TOLSF)) TOLSF = 1.0D0
						write(*,*) "Initial tolerance : ", ATOL_UP,RTOL_UP
						write(*,*) "Tolerance factor used ", TOLSF
						ATOL_UP = ATOL_UP * TOLSF * 2.0D0
						RTOL_UP = RTOL_UP * TOLSF * 2.0D0
						ISTATE = 3
						if (ATOL_UP(1) .GE. atol(1)*8.0D0) then; ITASK = 4; ISTATE = 2;RWORK(1) = tout; end if
					elseif (ISTATE == -1) then
						IWORK(6)=0 
						IWORK(6) = INT(IWORK(11)*1.5, i4b);
						! tout = ScattererParams(cSourceParams%iScCloud)%T_driv(iout) + 1D-15;
						! write(*,*) "MAXSTEPS, MXSTEPS USED BEFORE " , IWORK(6), IWORK(11)
						! write(*,*) "Updated tolerance : ", ATOL_UP,RTOL_UP
						ISTATE = 3
					elseif(ISTATE == -3) then 
						RWORK(5:7) = 0 ; IWORK(5:7) = 0;
						tout = RWORK(13)-RWORK(11)/2.0D0
						ATOL_UP = ATOL_UP  * 2.0D0
						RTOL_UP = RTOL_UP  * 2.0D0
						IWORK(6) = INT(IWORK(11) * 2, i4b);
						write(*,*) "TOUT CHANGED TO = ", tout,"ATOL = ", ATOL_UP
						write(*,*) "ISTATE = " ,ISTATE
						ISTATE = 3
					elseif(ISTATE == -4) then
						RWORK(5:7) = 0 ; IWORK(5:7) = 0;
						write(*,*) "RESET DEFAULT VALUES"
						write(*,*) "ISTATE = " ,ISTATE
						ISTATE = 3
					else
						write(*,*) "ISTATE = " ,ISTATE
						ISTATE = 3
					endif  
					 
					call dlsoda(MARMOTTANT,NEQ_ODEPACK,y,t,tout,itol,RTOL_UP,ATOL_UP,itask,istate,iopt,rwork,lrw,iwork,liw,JACDUM,jt)
					
					IF (ISTATE>0) RWORK((/5,7,8,9/)) = 0 
					IF (ISTATE>0) IWORK(5:6) = 0 
					
				end do
				R_bub(iout,:) = y(1:3)
				RealTimeOut(iout) = t
				
		enddo
	endif 
    norm_factor = ScattererParams(cSourceParams%iScCloud)%time_norm**(/0.0D0,1.0D0,2.0D0/)*ScattererParams(cSourceParams%iScCloud)%rad_norm
    R_bub(:,1) = (R_bub(:,1) +1.0D0)*norm_factor(1) ; R_bub(:,2) = R_bub(:,2)*norm_factor(2) ;R_bub(:,3) = R_bub(:,3)*norm_factor(3)
    RealTimeOut = RealTimeOut*1.0D0/ScattererParams(cSourceParams%iScCloud)%time_norm

	! Find volume acceleration d^2V/dt^2 [m^3/s^2] 
    ! This way the temporal derivative is calculated analytically so a spectral difference method is not needed.
    V_dd_pad = REAL(4.0D0*pi*R_bub(:,1)*(R_bub(:,1)*R_bub(:,3)+2.0D0*R_bub(:,2)**2),dp)!*dTaperSupportWindowN
	! P_bub  = REAL( ScattererParams(cSourceParams%iScCloud)%P_g0*( R_bub(:,1) /ScattererParams(cSourceParams%iScCloud)%R0(iBubble) )**(-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama) * (1.0D0-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama/cMediumParams%c0 * R_bub(:,2)),dp)
	
    DEALLOCATE(ScattererParams(cSourceParams%iScCloud)%P_driv,ScattererParams(cSourceParams%iScCloud)%T_driv)
	
    END SUBROUTINE RP_SOLVER

    SUBROUTINE MARMOTTANT(neq, t, R, Rdot)
    
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
    !   The subroutine MARMOTTANT_X solves the Marmottant model of the RP Equation
    !   using a simplification of R=x*R0 and t=τ/f0 . This way the radius is normalized and the
    !	time is normalized so the equation can handle the more unstable cases in a nicer way.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   neq		  i    i8b          number of equations (R_dot , R_dotdot)
    !   t		  r    dp           time stamp for solving the equation
    !   R		  r    dp           solution of the equation ( same length with y)
    !   Rdot		  r    dp           Rdot variable for the solution of the equation
    !   sigma_R		  r    dp           surface tension as a function of Radius
    !   P_elas		  r    dp           damping due to surface tension
    !   P_vis		  r    dp           damping due to the viscosity of the fluid 

    integer neq, iBubble
    real(dp) t, R(*), Rdot(*), sigma_R , R_norm(3)
	real(dp) P_interp(1)  , R_VanderWaals
	real(dp) P_elas, P_vis, P_gas, Damp_ac, Damp_visc, P_total 

    iBubble= NINT(R(4))
	ScattererParams(cSourceParams%iScCloud)%kappa_s  = (1.5D-9)*EXP(8.0D5*ScattererParams(cSourceParams%iScCloud)%R0(iBubble))
    ! call INTERP1D(ScattererParams(cSourceParams%iScCloud)%T_driv,ScattererParams(cSourceParams%iScCloud)%P_driv,real((/t/),dp), P_interp);
	 P_interp = NINT(R(5))
    ! In this method , the solver solves for x = R/R0 which is easier because it does not have to deal with really low numbers
    ! Accuracy meaning atol and rtol should be increased in this case ( Basically it is the division of the atol and rtol of the other method over R0)
 
	R_norm(1) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * real(R(1) + 1.0D0,dp)   ! Normalized R = r *  (1+x)
	R_norm(2) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm * R(2) ! Normalized R = r * tau*  (1+x)
	
	sigma_R = ScattererParams(cSourceParams%iScCloud)%chi*( R_norm(1)**2 /ScattererParams(cSourceParams%iScCloud)%R_b**2-1.0D0)
	
    if ( R_norm(1) .lt. ScattererParams(cSourceParams%iScCloud)%R_b) then                        	!rupture state
       sigma_R = 0.0D0
    else  if ( R_norm(1) .ge. ScattererParams(cSourceParams%iScCloud)%R_r) then 					!buckled state
        sigma_R = ScattererParams(cSourceParams%iScCloud)%sigma_w
    end if

	! call VanDerWaalsHardCoreRadius(R,R_VanderWaals)
	! P_gas     =  ScattererParams(cSourceParams%iScCloud)%P_g0*R_VanderWaals
    P_gas     =  ScattererParams(cSourceParams%iScCloud)%P_g0*( R_norm(1) /ScattererParams(cSourceParams%iScCloud)%R0(iBubble) )**(-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama)
	Damp_ac   =  1.0D0-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama/cMediumParams%c0 * R_norm(2)
	Damp_visc =  4.0D0*cMediumParams%mu * R_norm(2)/R_norm(1) 
    P_elas    =  2.0D0*sigma_R/R_norm(1)
    P_vis     =  4.0D0*ScattererParams(cSourceParams%iScCloud)%kappa_s * R_norm(2)/R_norm(1)**2
	P_total   = (P_gas*Damp_ac - cMediumParams%P0 - P_interp(1) - Damp_visc - P_elas - P_vis)
	
    Rdot(1) = R(2)  ! Rdot , radius velocity 
    R_norm(3) = (P_total/cMediumParams%rho0 - 3.0D0/2.0D0*R_norm(2)**2)/R_norm(1)  
	Rdot(2) = R_norm(3)/(ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm**2 )! Rddot , radius acceleration
    R(3) = Rdot(2)
    END SUBROUTINE MARMOTTANT
	
    !===================================== Experimental Solver ==========================================
    
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
    !   The subroutine R_EXP and Marmottant_EXP solves the ODE for the microbubble 
    !	response. The main difference is that in this case , the curve used for the
    !	surface tension is generated by experimental results taken from TU Twente.
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   neq		  i    i8b          number of equations (R_dot , R_dotdot)
    !   t		  r    dp           time stamp for solving the equation
    !   R		  r    dp           solution of the equation ( same length with y)
    !   Rdot		  r    dp           Rdot variable for the solution of the equation
    !   sigma_R		  r    dp           surface tension as a function of Radius
    !   P_elas		  r    dp           damping due to surface tension
    !   P_vis		  r    dp           damping due to the viscosity of the fluid 
    SUBROUTINE R_EXP(iBubble)
    
    ! ==================== Experimental values
    integer(i8b)		::  iBubble
    real(dp) 		 	::	Rexp(2*int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0D+9/2.0))
    real(qp)			::	sigma_RExp(2*int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0D+9/2.0))
    real(qp)			::	A(2*int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0D+9/2.0)), A0, A_m(2*int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0D+9/2.0))
    integer(i8b)		::	i, j, b
     
    ! ======================  Proceedure for A0 corrected
    Rexp = [( i*1.0D0/1000.0+ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0E+6 , i = -int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1E+9/2.0),int(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1E+9/2.0) )] 

    A0      = 4.0D0*pi*(ScattererParams(cSourceParams%iScCloud)%R0(iBubble)*1.0E+6)**2;
    A  	    = 4.0D0*pi*(Rexp)**2
    A_m     = A*1.0D0/A0
 
    ! First case is created by doing the multiplication for every R_test
    sigma_RExp = 0.0D0
    do i = 1,size(ScattererParams(cSourceParams%iScCloud)%coeff_fit)
    	sigma_RExp = sigma_RExp + ScattererParams(cSourceParams%iScCloud)%coeff_fit(i)*A_m**(11 - i)
    end do
    sigma_RExp = sigma_RExp*abs(A_m >= 0.9216D0 .AND. A_m <= 1.116D0) + 0.0D0 * abs(A_m < 0.9216D0) + ScattererParams(cSourceParams%iScCloud)%sigma_w * abs(A_m > 1.116D0)

	b = minloc(  abs(sigma_RExp - ScattererParams(cSourceParams%iScCloud)%sigma_R0) ,1)
    ScattererParams(cSourceParams%iScCloud)%A0c = 2.0D0*A0 - 4.0D0*pi*(Rexp(b))**2;
    
    END SUBROUTINE R_EXP
	
    SUBROUTINE EXP_SIGMA(neq, t, R, Rdot)
    ! =============================================================================
    !
    !   Programmer: Agisilaos Matalliotakis // Date : 211106
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
    !   The subroutine EXP_SIGMA solves the ODE for the microbubble 
    !	response. The main difference is that in this case , the curve used for the
    !	surface tension is generated by experimental results taken from TU Twente.
	!	Moreover, a normalization factor is assumed for R = R0*(1+x) and t = tau / y
    !
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   neq		  i    i8b          number of equations (R_dot , R_dotdot)
    !   t		  r    dp           time stamp for solving the equation
    !   R		  r    dp           solution of the equation ( same length with y)
    !   Rdot		  r    dp           Rdot variable for the solution of the equation
    !   sigma_R		  r    dp           surface tension as a function of Radius
    !   P_elas		  r    dp           damping due to surface tension
    !   P_vis		  r    dp           damping due to the viscosity of the fluid 
	integer 		:: 	neq
	integer(i8b) 	::  iBubble, i 
    real(qp)		::	A_R,A_m,sigma_Rcoeff(11), sigma_R
    real(dp) 		::	t, R(*), Rdot(*) , R_norm(3)
	real(qp) 		::	P_elas, P_vis, P_gas, Damp_ac, Damp_visc, P_total
	real(dp)		::	P_interp(1) , R_VanderWaals

    iBubble= NINT(R(4))
	ScattererParams(cSourceParams%iScCloud)%kappa_s  = (1.5D-9)*EXP(8.0D5*ScattererParams(cSourceParams%iScCloud)%R0(iBubble))
	call interp1D(ScattererParams(cSourceParams%iScCloud)%T_driv,ScattererParams(cSourceParams%iScCloud)%P_driv,real((/t/),dp), P_interp);

	R_norm(1) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * real(R(1) + 1.0D0,dp)   ! Normalized R = r *  (1+x)
	R_norm(2) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm * R(2) ! Normalized R = r * tau*  (1+x)
	
	! R = R0*(1+x)
    A_R = 4.0D0*pi*(R_norm(1) *1.0D+6)**2.0D0;
    A_m = A_R*1.0D0/ScattererParams(cSourceParams%iScCloud)%A0c;
    
    ! First case is created by doing the multiplication for every R_test
    sigma_R = 0.0D0
    do i = 1,size(ScattererParams(cSourceParams%iScCloud)%coeff_fit)
    	sigma_Rcoeff(i)  =  ScattererParams(cSourceParams%iScCloud)%coeff_fit(i)*A_m**(11 - i)
    end do
    sigma_R = sum(sigma_Rcoeff)
    sigma_R = sigma_R*abs(A_m >= 0.9216D0 .AND. A_m <= 1.116D0) + 0.0D0 * abs(A_m < 0.9216D0) + ScattererParams(cSourceParams%iScCloud)%sigma_w * abs(A_m > 1.116D0)

	call VanDerWaalsHardCoreRadius(R,R_VanderWaals)
	! P_gas     =  ScattererParams(cSourceParams%iScCloud)%P_g0*R_VanderWaals
    P_gas     =  ScattererParams(cSourceParams%iScCloud)%P_g0*( R_norm(1) /ScattererParams(cSourceParams%iScCloud)%R0(iBubble) )**(-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama)
	Damp_ac   =  1.0D0-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama/cMediumParams%c0 * R_norm(2)
	Damp_visc =  4.0D0*cMediumParams%mu * R_norm(2)/R_norm(1)
    P_elas    =  2.0D0*sigma_R/R_norm(1)
    P_vis     =  4.0D0*ScattererParams(cSourceParams%iScCloud)%kappa_s * R_norm(2)/R_norm(1)**2
	P_total   = real(P_gas*Damp_ac - cMediumParams%P0 - R(5) - Damp_visc - P_elas - P_vis,dp)
	
    Rdot(1) 	= R(2)  ! Rdot , radius velocity 
    R_norm(3) 	= (P_total/cMediumParams%rho0 - 3.0D0/2.0D0*R_norm(2)**2)/R_norm(1)  
	Rdot(2) 	= R_norm(3)/(ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm**2 )! Rddot , radius acceleration
    R(3) 		= Rdot(2)
    END SUBROUTINE EXP_SIGMA
	
	 SUBROUTINE MARMOTTANT_EXP_NNORM(neq, t, R, Rdot)
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
    !   The subroutine MARMOTTANT_EXP_NNORM solves the ODE for the microbubble 
    !	response. The marmottant model is used in order to describe the behaviour
	!	of 1 microbubble. Moreover, a normalization factor is assumed for R = R0*x 
	! 	and t = tau / y
	!	See "A model for large amplitude oscillations of coated bubbles accounting for buckling and rupture"
    !	Marmottant et. al.
	!
    ! *****************************************************************************
    !
    !   INPUT/OUTPUT PARAMETERS
    !
    !   neq		  i    i8b          number of equations (R_dot , R_dotdot)
    !   t		  r    dp           time stamp for solving the equation
    !   R		  r    dp           solution of the equation ( same length with y)
    !   Rdot		  r    dp           Rdot variable for the solution of the equation
    !   sigma_R		  r    dp           surface tension as a function of Radius
    !   P_elas		  r    dp           damping due to surface tension
    !   P_vis		  r    dp           damping due to the viscosity of the fluid 
    integer neq, iBubble
    real(dp) t, R(*), Rdot(*), sigma_R , R_norm(3)
	real(dp) P_interp(1) 
	real(dp) P_elas, P_vis, P_gas, Damp_ac, Damp_visc, P_total 

    iBubble= NINT(R(4))
	ScattererParams(cSourceParams%iScCloud)%kappa_s  = (1.5D-9)*EXP(8.0D5*ScattererParams(cSourceParams%iScCloud)%R0(iBubble))
    ! call interp1D(ScattererParams(cSourceParams%iScCloud)%T_driv,ScattererParams(cSourceParams%iScCloud)%P_driv,real((/t/),dp), P_interp);
	
    ! In this method , the solver solves for x = R/R0 which is easier because it does not have to deal with really low numbers
    ! Accuracy meaning atol and rtol should be increased in this case ( Basically it is the division of the atol and rtol of the other method over R0)
 
	R_norm(1) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * real(R(1),dp)   ! Normalized R = r *  (1+x)
	R_norm(2) =  ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm * R(2) ! Normalized R = r * tau*  (1+x)
	
	sigma_R = ScattererParams(cSourceParams%iScCloud)%chi*( R_norm(1)**2 /ScattererParams(cSourceParams%iScCloud)%R_b**2-1.0D0)
	
    if ( R_norm(1) .lt. ScattererParams(cSourceParams%iScCloud)%R_b) then                        	!rupture state
       sigma_R = 0.0D0
    else  if ( R_norm(1) .ge. ScattererParams(cSourceParams%iScCloud)%R_r) then 																	!buckled state
        sigma_R = ScattererParams(cSourceParams%iScCloud)%sigma_w
    end if
	

    P_gas     =  ScattererParams(cSourceParams%iScCloud)%P_g0*( R_norm(1) /ScattererParams(cSourceParams%iScCloud)%R0(iBubble) )**(-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama)
	Damp_ac   =  1.0D0-3.0D0*ScattererParams(cSourceParams%iScCloud)%gama/cMediumParams%c0 * R_norm(2)
	Damp_visc =  4.0D0*cMediumParams%mu * R_norm(2)/R_norm(1)
    P_elas    =  2.0D0*sigma_R/R_norm(1)
    P_vis     =  4.0D0*ScattererParams(cSourceParams%iScCloud)%kappa_s * R_norm(2)/R_norm(1)**2
	P_total   = (P_gas*Damp_ac - cMediumParams%P0 - R(5) - Damp_visc - P_elas - P_vis)
	
    Rdot(1) = R(2)  ! Rdot , radius velocity 
    R_norm(3) = (P_total/cMediumParams%rho0 - 3.0D0/2.0D0*R_norm(2)**2)/R_norm(1)  
	Rdot(2) = R_norm(3)/(ScattererParams(cSourceParams%iScCloud)%rad_norm * ScattererParams(cSourceParams%iScCloud)%time_norm**2 )! Rddot , radius acceleration
    R(3) = Rdot(2)
    END SUBROUTINE MARMOTTANT_EXP_NNORM
    
	
	SUBROUTINE VanDerWaalsHardCoreRadius(R,R_VanderWaals)
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
    !   The subroutine VanDerWaalsHardCoreRadius includes a change in the radius of 
	! 	the pressure inside the microbubble using the Van der Waals hardcore radius.
	!	see "Phase diagrams for sonoluminescing bubbles", Hilgenfeldt et. al.
    !
    ! *****************************************************************************
    
    real(dp),intent(out) :: R_VanderWaals
    real(dp),intent(in)  ::	R(*)
    real(dp) 			 :: h, Nominator, Denominator
	integer				 :: iBubble
	iBubble = NINT(R(4))
	
	h = ScattererParams(cSourceParams%iScCloud)%R0(iBubble)/8.86D0
	Nominator   = (ScattererParams(cSourceParams%iScCloud)%R0(iBubble)			  )**3 - h**3
	Denominator = (ScattererParams(cSourceParams%iScCloud)%rad_norm * (R(1) + 1) )**3 - h**3
	R_VanderWaals = (Nominator/Denominator) ** (ScattererParams(cSourceParams%iScCloud)%gama)
	
    END SUBROUTINE VanDerWaalsHardCoreRadius


!    !===================================== END OF Experimental Solver ==========================================
!    
   ! SUBROUTINE JAC1 (neq, t, R, ml, mu, pd, nrowpd)
   ! integer 			::	neq, ml, mu, nrowpd, iBubble
   ! real(dp) 		::	t, pd(nrowpd,2) , R(4), sigma_R
   ! real(dp) 		::	P_elas, P_vis, P_gas, Damp_ac, Damp_visc, P_interp(1) 
   
   ! iBubble = NINT(R(4))

   ! if (R(1) .lt. ScattererParams(cSourceParams%iScCloud)%R_b) then                               !buckled state
       ! sigma_R = 0.0
   ! else if (R(1) .gt. ScattererParams(cSourceParams%iScCloud)%R_r) then                        !rupture state
       ! sigma_R = ScattererParams(cSourceParams%iScCloud)%sigma_w
   ! else                                            !Elastic state
       ! sigma_R = ScattererParams(cSourceParams%iScCloud)%chi*((R(1)**2)/(ScattererParams(cSourceParams%iScCloud)%R_b**2)-1.0)
   ! end if

   ! P_elas =  2*sigma_R*1.0_dp/R(1)
   ! P_vis  =  4*ScattererParams(cSourceParams%iScCloud)%kappa_s*R(2)*1.0_dp/(R(1)**2)

   ! pd(1,1) = 0.0d0
   ! pd(1,2) = 1.0d0
   ! pd(2,1) = (-(3*ScattererParams(cSourceParams%iScCloud).gama)*(1 - 3*ScattererParams(cSourceParams%iScCloud).gama*R(2)/cMediumParams%c0)*(ScattererParams(cSourceParams%iScCloud).P_g0*(ScattererParams(cSourceParams%iScCloud).R0(iBubble)/R(1)) &
       ! **((3*ScattererParams(cSourceParams%iScCloud).gama)*(1 - 3*ScattererParams(cSourceParams%iScCloud).gama*R(2)/cMediumParams%c0)-1)  + 4*cMediumParams%mu*R(2)/(R(1)**2)+ &
       ! 2*sigma_R/(R(1) **2) + 8*ScattererParams(cSourceParams%iScCloud).kappa_s*R(2)/(R(1)**3))/cMediumParams%rho0- 3.0/2.0*R(2)**2)/R(1) &
       ! -((ScattererParams(cSourceParams%iScCloud).P_g0*(ScattererParams(cSourceParams%iScCloud).R0(iBubble)/R(1))**(3*ScattererParams(cSourceParams%iScCloud).gama)*(1-3*ScattererParams(cSourceParams%iScCloud).gama*R(2)/cMediumParams%c0) - &
       ! cMediumParams%P0- R(3) -4*cMediumParams%mu*R(2)/R(1)-P_elas-P_vis)/cMediumParams%rho0-3/2*R(2)**2)/(R(1)**2)
   ! pd(2,2) = (ScattererParams(cSourceParams%iScCloud).P_g0*(ScattererParams(cSourceParams%iScCloud).R0(iBubble)/R(1))**(3*ScattererParams(cSourceParams%iScCloud).gama)*( -3*ScattererParams(cSourceParams%iScCloud).gama/cMediumParams%c0)-  &
       ! 4*cMediumParams%mu/R(1) - 4*ScattererParams(cSourceParams%iScCloud).kappa_s/(R(1)**2))/cMediumParams%rho0 - 3*R(2)
   ! return
   ! END SUBROUTINE JAC1
   SUBROUTINE JACDUM
   implicit none 
   return
   end
    END MODULE RP



