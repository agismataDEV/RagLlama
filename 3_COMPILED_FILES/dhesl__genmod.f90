        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DHESL__genmod
          INTERFACE 
            SUBROUTINE DHESL(A,LDA,N,IPVT,B)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(*)
              REAL(KIND=8) :: B(*)
            END SUBROUTINE DHESL
          END INTERFACE 
        END MODULE DHESL__genmod
