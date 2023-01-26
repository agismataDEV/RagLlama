        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DLSODAR__genmod
          INTERFACE 
            SUBROUTINE DLSODAR(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,     &
     &ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,JT,G,NG,JROOT)
              INTEGER(KIND=4) :: NG
              INTEGER(KIND=4) :: LIW
              INTEGER(KIND=4) :: LRW
              EXTERNAL F
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
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
              EXTERNAL JAC
              INTEGER(KIND=4) :: JT
              EXTERNAL G
              INTEGER(KIND=4) :: JROOT(NG)
            END SUBROUTINE DLSODAR
          END INTERFACE 
        END MODULE DLSODAR__genmod
