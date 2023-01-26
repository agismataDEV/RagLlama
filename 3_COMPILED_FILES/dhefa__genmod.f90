        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DHEFA__genmod
          INTERFACE 
            SUBROUTINE DHEFA(A,LDA,N,IPVT,INFO,JOB)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: A(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: INFO
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DHEFA
          END INTERFACE 
        END MODULE DHEFA__genmod
