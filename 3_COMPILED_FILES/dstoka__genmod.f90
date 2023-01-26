        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSTOKA__genmod
          INTERFACE 
            SUBROUTINE DSTOKA(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVX,ACOR,WM,IWM&
     &,F,JAC,PSOL)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: YH1(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: SAVX(*)
              REAL(KIND=8) :: ACOR(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL F
              EXTERNAL JAC
              EXTERNAL PSOL
            END SUBROUTINE DSTOKA
          END INTERFACE 
        END MODULE DSTOKA__genmod
