        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSOLPK__genmod
          INTERFACE 
            SUBROUTINE DSOLPK(NEQ,Y,SAVF,X,EWT,WM,IWM,F,PSOL)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL F
              EXTERNAL PSOL
            END SUBROUTINE DSOLPK
          END INTERFACE 
        END MODULE DSOLPK__genmod
