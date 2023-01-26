        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGBSL__genmod
          INTERFACE 
            SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: ABD(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ML
              INTEGER(KIND=4) :: MU
              INTEGER(KIND=4) :: IPVT(*)
              REAL(KIND=8) :: B(*)
              INTEGER(KIND=4) :: JOB
            END SUBROUTINE DGBSL
          END INTERFACE 
        END MODULE DGBSL__genmod
