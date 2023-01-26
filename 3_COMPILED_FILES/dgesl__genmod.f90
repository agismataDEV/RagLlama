        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGESL__genmod
          INTERFACE 
            SUBROUTINE DGESL(A,LDA,N,IPVT,B,JOB)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(*)
              REAL(KIND=8) :: B(*)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DGESL
          END INTERFACE 
        END MODULE DGESL__genmod
