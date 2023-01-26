        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DEWSET__genmod
          INTERFACE 
            SUBROUTINE DEWSET(N,ITOL,RTOL,ATOL,YCUR,EWT)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ITOL
              REAL(KIND=8) :: RTOL(*)
              REAL(KIND=8) :: ATOL(*)
              REAL(KIND=8) :: YCUR(N)
              REAL(KIND=8) :: EWT(N)
            END SUBROUTINE DEWSET
          END INTERFACE 
        END MODULE DEWSET__genmod
