        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGEFA__genmod
          INTERFACE 
            SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGEFA
          END INTERFACE 
        END MODULE DGEFA__genmod
