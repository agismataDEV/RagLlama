        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSETPK__genmod
          INTERFACE 
            SUBROUTINE DSETPK(NEQ,Y,YSV,EWT,FTEM,SAVF,JOK,WM,IWM,F,JAC)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YSV(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: FTEM(*)
              REAL(KIND=8) :: SAVF(*)
              INTEGER(KIND=4) :: JOK
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL F
              EXTERNAL JAC
            END SUBROUTINE DSETPK
          END INTERFACE 
        END MODULE DSETPK__genmod
