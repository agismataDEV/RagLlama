        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLSODIS__genmod
          INTERFACE 
            SUBROUTINE DLSODIS(RES,ADDA,JAC,NEQ,Y,YDOTI,T,TOUT,ITOL,RTOL&
     &,ATOL,ITASK,ISTATE,IOPT,RWORK,LRW,IWORK,LIW,MF)
              INTEGER(KIND=4) :: LIW
              INTEGER(KIND=4) :: LRW
              EXTERNAL RES
              EXTERNAL ADDA
              EXTERNAL JAC
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YDOTI(*)
              REAL(KIND=8) :: T
              REAL(KIND=8) :: TOUT
              INTEGER(KIND=4) :: ITOL
              REAL(KIND=8) :: RTOL(*)
              REAL(KIND=8) :: ATOL(*)
              INTEGER(KIND=4) :: ITASK
              INTEGER(KIND=4) :: ISTATE
              INTEGER(KIND=4) :: IOPT
              REAL(KIND=8) :: RWORK(LRW)
              INTEGER(KIND=4) :: IWORK(LIW)
              INTEGER(KIND=4) :: MF
            END SUBROUTINE DLSODIS
          END INTERFACE 
        END MODULE DLSODIS__genmod
