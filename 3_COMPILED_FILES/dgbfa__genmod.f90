        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DGBFA__genmod
          INTERFACE 
            SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
              INTEGER(KIND=4) :: LDA
              REAL(KIND=8) :: ABD(LDA,*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: ML
              INTEGER(KIND=4) :: MU
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: INFO
            END SUBROUTINE DGBFA
          END INTERFACE 
        END MODULE DGBFA__genmod
