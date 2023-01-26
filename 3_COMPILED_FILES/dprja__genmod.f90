        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPRJA__genmod
          INTERFACE 
            SUBROUTINE DPRJA(NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WM,IWM,F,JAC)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: FTEM(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL F
              EXTERNAL JAC
            END SUBROUTINE DPRJA
          END INTERFACE 
        END MODULE DPRJA__genmod
