    MODULE RP
	USE ParnacGeneral
    USE Specialfun
    USE ParnacDataDef
    USE ParnacParamDef
    USE ParnacParamInit, ONLY : &
        ScattererInit
    USE Constants
	 USE ParnacTransformFilter, ONLY : &
       dTaperingWindow, INTERP1DFREQ
    
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
    real(dp), intent(out)    ::  V_dd_pad(n_samples), RealTimeOut(n_samples)
	
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
	
    integer(i8b), parameter ::  NEQ_ODEPACK(1) = 2 , NEQN_RK = 4			!   Number of first order ode's
    integer(i4b), parameter ::  lrw = 697			!   LRW   :IN     Declared length of RWORK (in user's DIMENSION statement).
    integer(i4b), parameter ::  liw = 45            !   LIW   :IN     Declared length of IWORK (in user's DIMENSION statement).
    real(sp)    			::  atol(1), rtol(1), dtout, t, tout, y(NEQN_RK), yp(NEQN_RK) , TOLSF, ATOL_UP(1), RTOL_UP(1)
    real(sp)				::  rwork(lrw)
    integer(i4b) 			::  iwork(liw)
    integer(i4b) 			::  i, iopar, iopt, iout, istate, itask, itol, mf, jt, flag
	real(dp) 				::  norm_factor(3), P_interp(1)
	real(dp)				::  R_bub(n_samples,3), dTaperSupportWindowN(n_samples)
	
	character(len = 1024)	::  	actemp
	
	logical					:: log_isnan
    
    call ScattererInit()
    call R_exp(iBubble)
    
	flag = 1 			!  This is for RK
    itol = 1            !  if atol scalar, itol = 1 and  if atol array, itol = 2	
    mf = 12             !   memory flag :10 -(Adams) - Non-stiff code
    jt  = 2 			!   Jacobian Type : 2 - For DLSODA instead of mf
    iopt = 1            !   0 for not optional inputs, 1 for optional inputs
    itask = 1
    istate = 1
    ! IWORK(6)=0
	! IWORK(6) = 2000; 
	IWORK(7) = 0;
	IWORK(7) = 1;
	! CALL XSETF(0) ! 0 if  no messages should be printed in ODEPACK call
    
	! If ScattererParams%rad_norm = ScattererParams%R0 , The result of the RP Solver is normalized to x = R/ScattererParams%rad_norm = R/R0 , so the solved value is multiplied by ScattererParams%R0 at the end of the code
    ! atol and rtol shoudld be 1E-5
    ! If ScattererParams%rad_norm =  1.0D0, the result of the RP solver ir R, so the solved value should not be multiplied with any factor
    ! atol and rtol shoudld be 1E-9. DLSODE input should change to MARMOTTANT_R

	if (trim(ScattererParams%Solver_Normalize) =='time') ScattererParams%time_norm = cModelParams%freq0
	if (trim(ScattererParams%Solver_Normalize) =='radius') ScattererParams%rad_norm  = ScattererParams%R0(iBubble)
	if (trim(ScattererParams%Solver_Normalize) =='freq0_and_radius') then
		 ScattererParams%time_norm = cModelParams%freq0
		 ScattererParams%rad_norm  = ScattererParams%R0(iBubble)
	elseif (trim(ScattererParams%Solver_Normalize) =='Minnaert_and_radius') then
		 ScattererParams%time_norm = SQRT(cMediumParams%P0/cMediumParams%rho0/ScattererParams%R0(iBubble)**2)
		 ScattererParams%rad_norm  = ScattererParams%R0(iBubble)
	endif
	 
    rtol = 10.0**(floor(log10(ScattererParams%R0(iBubble)/ScattererParams%rad_norm))-6.0D0)
    atol = 10.0**(floor(log10(ScattererParams%R0(iBubble)/ScattererParams%rad_norm))-6.0D0)
	
	dTaperSupportWindowN = dTaperingWindow(n_samples,(RealTimeIn(2)-RealTimeIn(1))* cModelParams%freq0,2.0_dp,2.0_dp)
	
	ScattererParams%P_driv = RealPressIn * dTaperSupportWindowN
	ScattererParams%T_driv = RealTimeIn * ScattererParams%time_norm;
		  
	! Medium parameters (water, Room temperature =20° and 1 atm ambient pressure) 
    ScattererParams.P_g0 = cMediumParams%P0+2.0D0*ScattererParams.sigma_R0/ScattererParams.R0(iBubble)

    ! Marmottant model parameters
    ScattererParams.R_b =  ScattererParams.R0(iBubble)*1.0D0/sqrt(ScattererParams.sigma_R0*1.0D0/ScattererParams.chi+1.0D0)
    ScattererParams.R_r =  ScattererParams.R_b*sqrt(ScattererParams.sigma_w*1.0D0/ScattererParams.chi+1.0D0)  
     
    y(1) = ScattererParams%R0(iBubble)/ScattererParams%rad_norm - 1.0D0                    
    y(2) = 0.0d0 ; y(3) = 0.0d0 ; y(4) = iBubble*1.0D0;
    R_bub(1,:) = y(1:3)
    
    tout = ScattererParams%T_driv(1)
    t = ScattererParams%T_driv(1)
	RealTimeOut = 0;
    RealTimeOut(1) = ScattererParams%T_driv(1)
	if (ScattererParams%Solver_Method=='RK') then !For the case of Runge Kutta Solver
		do iout = 2,n_samples
			tout = RealTimeIn(iout) ; 
			call R8_RKF45(MARMOTTANT_NORM, NEQN_RK, y, yp, t, tout, rtol, atol, flag )
			R_bub(iout,:) = y(1:3) 
			RealTimeOut(iout) = t 
		enddo
	else		!For the case of ODEPACK solver
    	do iout = 2,n_samples
				ATOL_UP = ATOL
				RTOL_UP = RTOL
			    ! t  = RealTimeIn(iout-1) 
				tout = ScattererParams%T_driv(iout)
				call slsoda(MARMOTTANT_NORM,NEQ_ODEPACK,y,t,tout,itol,rtol,atol,itask,istate,iopt,rwork,lrw,iwork,liw,JAC1,jt)
				do while (ISTATE <0)
					if (ISTATE == -2) then
						TOLSF = RWORK(14) 
						write(*,*) "Initial tolerance : ", ATOL_UP,RTOL_UP
						write(*,*) "Tolerance factor used ", TOLSF
						ATOL_UP = ATOL_UP * TOLSF * 2.0D0
						RTOL_UP = RTOL_UP * TOLSF * 2.0D0
						write(*,*) "Updated tolerance : ", ATOL_UP,RTOL_UP
					elseif (ISTATE == -1) then
						IWORK(6)=0
						IWORK(6) = INT(IWORK(11) * 2, i4b);
						ATOL_UP = ATOL_UP * 2.0D0
						RTOL_UP = RTOL_UP * 2.0D0
						write(*,*) "MAXSTEPS, MXSTEPS USED BEFORE " , IWORK(6), IWORK(11)
						write(*,*) "Updated tolerance : ", ATOL_UP,RTOL_UP
					elseif(ISTATE == -3) then
						! RWORK(6)  = ScattererParams%T_driv(2) - ScattererParams%T_driv(1)
						RWORK(5:7) = 0 ; IWORK(5:7) = 0;
						tout = MAX(ScattererParams%T_driv(iout), RealTimeOut(iout-1))
						write(*,*) "TOUT CHANGED TO = ", tout
						write(*,*) "ISTATE = " ,ISTATE
					elseif(ISTATE == -4) then
						RWORK(5:7) = 0 ; IWORK(5:7) = 0;
						ATOL_UP = ATOL
						RTOL_UP = RTOL
						write(*,*) "RESET DEFAULT VALUES"
						write(*,*) "ISTATE = " ,ISTATE
					else
						write(*,*) "ISTATE = " ,ISTATE
					endif
					ISTATE = 3
					call slsoda(MARMOTTANT_NORM,NEQ_ODEPACK,y,t,tout,itol,RTOL_UP,ATOL_UP,itask,istate,iopt,rwork,lrw,iwork,liw,JAC1,jt)
					IF (ISTATE>0) RWORK(5:7:2) = 0 
					IF (ISTATE>0) IWORK(5:6) = 0 
				end do
				R_bub(iout,:) = y(1:3)
				RealTimeOut(iout) = t
		enddo
	endif 
	
    norm_factor = ScattererParams%time_norm**(/0.0D0,1.0D0,2.0D0/)*ScattererParams%rad_norm
    R_bub(:,1) = (R_bub(:,1) + 1.0D0)*norm_factor(1) ; R_bub(:,2) = R_bub(:,2)*norm_factor(2) ;R_bub(:,3) = R_bub(:,3)*norm_factor(3)
    RealTimeOut = RealTimeOut*1.0D0/ScattererParams%time_norm
	
	! Check to see if results are realistic
	log_isnan = sum(isnan(R_bub(:,1))*-1) + sum(isnan(R_bub)*-1) + sum(isnan(R_bub(:,3))*-1)
	if ( log_isnan>0 .OR. REAL(maxval(abs(R_bub)),dp)>1.0D+32 .OR. maxval(abs(R_bub)) > HUGE(maxval(abs(R_bub))) ) then ! Check Nan , Large values or Infinity
	    write(acTemp,"('Error in ODE SOLVER , Bubble No. ', I7, ' .')") iBubble
        call PrintToLog(acTemp,1) 
        stop
	endif
	
	! Filter the edges (not absolutely needed)
	R_bub(:,1)= R_bub(:,1) * dTaperSupportWindowN
	R_bub(:,2)= R_bub(:,2) * dTaperSupportWindowN
	R_bub(:,3)= R_bub(:,3) * dTaperSupportWindowN
		  
	! Find volume acceleration d^2V/dt^2 [m^3/s^2] 
    ! This way the temporal derivative is calculated analytically so a spectral difference method is not needed.
    V_dd_pad = 4.0D0*pi*R_bub(:,1)*(R_bub(:,1)*R_bub(:,3)+2.0D0*R_bub(:,2)**2)

    END SUBROUTINE RP_SOLVER

    SUBROUTINE MARMOTTANT_NORM(neq, t, R, Rdot)
    
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
    real(sp) t, R(*), Rdot(*), sigma_R
	real(dp) P_interp(1) 
	real(dp) P_elas, P_vis, P_gas, Damp_ac, Damp_visc

    iBubble= NINT(R(4))
    call interp1D(ScattererParams%T_driv,ScattererParams%P_driv,real((/t/),dp), P_interp);
    ! In this method , the solver solves for x = R/R0 which is easier because it does not have to deal with really low numbers
    ! Accuracy meaning atol and rtol should be increased in this case ( Basically it is the division of the atol and rtol of the other method over R0)

    if ( (R(1)+1) *ScattererParams%rad_norm .lt. ScattererParams%R_b) then                        	!rupture state
       sigma_R = 0.0D0
    else if ( (R(1)+1) *ScattererParams%rad_norm .ge. ScattererParams%R_b .AND. &
    		  (R(1)+1) *ScattererParams%rad_norm .le. ScattererParams%R_r) then 					!Elastic state
        sigma_R = ScattererParams.sigma_R0 + ScattererParams%chi*((ScattererParams%rad_norm * (R(1)+1) /ScattererParams%R_b)**2-1.0D0)
    else  																				!buckled state
        sigma_R = ScattererParams%sigma_w
    end if
	
    P_gas     =  ScattererParams%P_g0*(ScattererParams%R0(iBubble)/ScattererParams%rad_norm * 1.0D0/(R(1) + 1) )**(3.0D0*ScattererParams%gama)
	Damp_ac   =  1.0D0-3.0D0*ScattererParams%gama/cMediumParams%c0*ScattererParams%rad_norm*ScattererParams%time_norm*R(2)
	Damp_visc =  4.0D0*cMediumParams%mu*R(2)*ScattererParams%time_norm/(1+R(1))
    P_elas    =  2.0D0*sigma_R/(ScattererParams%rad_norm * (R(1)+1) )
	ScattererParams%kappa_s  = (1.5D-9)*EXP(8.0D5*ScattererParams%R0(iBubble))
    P_vis     =  4.0D0*ScattererParams%kappa_s*R(2)*ScattererParams%time_norm/(ScattererParams%rad_norm * (R(1) + 1.0D0) ** 2)
	
    Rdot(1) = R(2)  
    Rdot(2) = ((P_gas*Damp_ac - cMediumParams%P0 - R(5) - Damp_visc - P_elas - P_vis)/(cMediumParams%rho0 * ScattererParams%rad_norm**2* ScattererParams%time_norm**2)-(3.0D0/2.0D0)*R(2)**2)*1.0D0/(R(1)+1)
    R(3) = Rdot(2)
    END SUBROUTINE MARMOTTANT_NORM
    
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
    real(dp) 		 	::	Rexp(2*int(ScattererParams%R0(iBubble)*1.0D+9/2.0))
    real(qp)			::	sigma_RExp(2*int(ScattererParams%R0(iBubble)*1.0D+9/2.0))
    real(qp)			::	A(2*int(ScattererParams%R0(iBubble)*1.0D+9/2.0)), A0, A_m(2*int(ScattererParams%R0(iBubble)*1.0D+9/2.0))
    integer(i8b)		::	i, j, b
    
    ! ======================  Proceedure for A0 corrected
    Rexp = [( i*1.0D0/1000.0+ScattererParams%R0(iBubble)*1.0E+6 , i = -int(ScattererParams%R0(iBubble)*1E+9/2.0),int(ScattererParams%R0(iBubble)*1E+9/2.0) )] 

    A0      = 4.0D0*pi*(ScattererParams%R0(iBubble)*1.0E+6)**2;
    A  	    = 4.0D0*pi*(Rexp)**2
    A_m     = A*1.0D0/A0
 
    ! First case is created by doing the multiplication for every R_test
    sigma_RExp = 0.0D0
    do i = 1,size(ScattererParams%coeff_fit)
    	sigma_RExp = sigma_RExp + ScattererParams%coeff_fit(i)*A_m**(11 - i)
    end do
    sigma_RExp = sigma_RExp*abs(A_m >= 0.9216D0 .AND. A_m <= 1.116D0) + 0.0D0 * abs(A_m < 0.9216D0) + ScattererParams%sigma_w * abs(A_m > 1.116D0)

	b = minloc(  abs(sigma_RExp - ScattererParams%sigma_R0) ,1)
    ScattererParams%A0c = 2.0D0*A0 - 4.0D0*pi*(Rexp(b))**2;
    
    END SUBROUTINE R_EXP
    
    SUBROUTINE MARMOTTANT_EXP(neq, t, R, Rdot)
	integer 		:: 	neq
	integer(i8b) 	::  iBubble, i 
    real(qp)		::	A_R,A_m,sigma_Rcoeff(11), sigma_R
    real(sp) 		::	t, R(*), Rdot(*)
	real(dp) 		::	P_elas, P_vis, P_gas, Damp_ac, Damp_visc, P_interp(1) 

    iBubble= NINT(R(4))
	call interp1D(ScattererParams%T_driv,ScattererParams%P_driv,real((/t/),dp), P_interp);
	
    A_R = 4.0D0*pi*(ScattererParams%R0(iBubble)*R(1)*1.0D+6)**2.0D0;
    A_m = A_R*1.0D0/ScattererParams%A0c;
    
    ! First case is created by doing the multiplication for every R_test
    sigma_R = 0.0D0
    do i = 1,size(ScattererParams%coeff_fit)
    	sigma_Rcoeff(i)  =  ScattererParams%coeff_fit(i)*A_m**(11 - i)
    end do
    sigma_R = sum(sigma_Rcoeff)
    sigma_R = sigma_R*abs(A_m >= 0.9216D0 .AND. A_m <= 1.116D0) + 0.0D0 * abs(A_m < 0.9216D0) + ScattererParams%sigma_w * abs(A_m > 1.116D0)

	! COmpute the pressure values for each component in RP Equation
	P_gas     =  ScattererParams%P_g0*(ScattererParams%R0(iBubble)/ScattererParams%rad_norm * 1.0D0/(R(1) + 1) )**(3.0D0*ScattererParams%gama)
	Damp_ac   =  1.0D0-3.0D0*ScattererParams%gama/cMediumParams%c0*ScattererParams%rad_norm*R(2)*ScattererParams%time_norm
	Damp_visc =  4.0D0*cMediumParams%mu*ScattererParams%time_norm * R(2)/(1+R(1))
    P_elas    =  2.0D0*sigma_R/(ScattererParams%rad_norm * (R(1)+1) )
	ScattererParams%kappa_s  = (1.5D-9)*EXP(8.0D5*ScattererParams%R0(iBubble))
    P_vis     =  4.0D0*ScattererParams%kappa_s*ScattererParams%time_norm/ScattererParams%rad_norm * R(2) / (R(1) + 1) ** 2
	
    Rdot(1) = R(2)  
    Rdot(2) = ((P_gas*Damp_ac - cMediumParams%P0 - P_interp(1) - Damp_visc - P_elas - P_vis)/(cMediumParams%rho0 * ScattererParams%rad_norm**2* ScattererParams%time_norm**2)-(3.0D0/2.0D0)*R(2)**2)*1.0D0/(R(1)+1)
    R(3) = Rdot(2)
    END SUBROUTINE MARMOTTANT_EXP
    
!    !===================================== END OF Experimental Solver ==========================================
!    
!    SUBROUTINE JAC1 (neq, t, R, ml, mu, pd, nrowpd)
!    integer neq, ml, mu, nrowpd
!    double precision t, pd(nrowpd,neq)

!    if (R(1) .lt. ScattererParams%R_b) then                               !buckled state
!        sigma_R = 0.0
!    else if (R(1) .gt. ScattererParams%R_r) then                        !rupture state
!        sigma_R = ScattererParams%sigma_w
!    else                                            !Elastic state
!        sigma_R = ScattererParams%chi*((R(1)**2)/(ScattererParams%R_b**2)-1.0)
!    end if

!    P_elas =  2*sigma_R*1.0_dp/R(1)
!    P_vis  =  4*ScattererParams%kappa_s*R(2)*1.0_dp/(R(1)**2)

!    pd(1,1) = 0.0d0
!    pd(1,2) = 1.0d0
!    pd(2,1) = (-(3*ScattererParams.gama)*(1 - 3*ScattererParams.gama*R(2)/cMediumParams%c0)*(ScattererParams.P_g0*(ScattererParams.R0/R(1)) &
!        **((3*ScattererParams.gama)*(1 - 3*ScattererParams.gama*R(2)/cMediumParams%c0)-1)  + 4*cMediumParams%mu*R(2)/(R(1)**2)+ &
!        2*sigma_R/(R(1) **2) + 8*ScattererParams.kappa_s*R(2)/(R(1)**3))/cMediumParams%rho0- 3.0/2.0*R(2)**2)/R(1) &
!        -((ScattererParams.P_g0*(ScattererParams.R0/R(1))**(3*ScattererParams.gama)*(1-3*ScattererParams.gama*R(2)/cMediumParams%c0) - &
!        cMediumParams%P0- R(3) -4*cMediumParams%mu*R(2)/R(1)-P_elas-P_vis)/cMediumParams%rho0-3/2*R(2)**2)/(R(1)**2)
!    pd(2,2) = (ScattererParams.P_g0*(ScattererParams.R0/R(1))**(3*ScattererParams.gama)*( -3*ScattererParams.gama/cMediumParams%c0)-  &
!        4*cMediumParams%mu/R(1) - 4*ScattererParams.kappa_s/(R(1)**2))/cMediumParams%rho0 - 3*R(2)
!    return
!    END SUBROUTINE JAC1
    END MODULE RP



