        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DATP__genmod
          INTERFACE 
            SUBROUTINE DATP(NEQ,Y,SAVF,P,WGHT,HL0,WK,F,W)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: P(*)
              REAL(KIND=8) :: WGHT(*)
              REAL(KIND=8) :: HL0
              REAL(KIND=8) :: WK(*)
              EXTERNAL F
              REAL(KIND=8) :: W(*)
            END SUBROUTINE DATP
          END INTERFACE 
        END MODULE DATP__genmod
