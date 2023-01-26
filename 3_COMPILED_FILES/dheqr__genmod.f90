        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DHEQR__genmod
          INTERFACE 
            SUBROUTINE DHEQR(A,LDA,N,Q,INFO,IJOB)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: Q(*)
              INTEGER(KIND=4) :: INFO
              INTEGER(KIND=4) :: IJOB
            END SUBROUTINE DHEQR
          END INTERFACE 
        END MODULE DHEQR__genmod
