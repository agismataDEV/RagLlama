        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLHIN__genmod
          INTERFACE 
            SUBROUTINE DLHIN(NEQ,N,T0,Y0,YDOT,F,TOUT,UROUND,EWT,ITOL,   &
     &ATOL,Y,TEMP,H0,NITER,IER)
              INTEGER(KIND=4) :: NEQ(*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: T0
              REAL(KIND=8) :: Y0(*)
              REAL(KIND=8) :: YDOT(*)
              EXTERNAL F
              REAL(KIND=8) :: TOUT
              REAL(KIND=8) :: UROUND
              REAL(KIND=8) :: EWT(*)
              INTEGER(KIND=4) :: ITOL
              REAL(KIND=8) :: ATOL(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: TEMP(*)
              REAL(KIND=8) :: H0
              INTEGER(KIND=4) :: NITER
              INTEGER(KIND=4) :: IER
            END SUBROUTINE DLHIN
          END INTERFACE 
        END MODULE DLHIN__genmod
