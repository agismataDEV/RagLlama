        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DHELS__genmod
          INTERFACE 
            SUBROUTINE DHELS(A,LDA,N,Q,B)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: Q(*)
              REAL(KIND=8) :: B(*)
            END SUBROUTINE DHELS
          END INTERFACE 
        END MODULE DHELS__genmod
