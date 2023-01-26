        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSTODI__genmod
          INTERFACE 
            SUBROUTINE DSTODI(NEQ,Y,YH,NYH,YH1,EWT,SAVF,SAVR,ACOR,WM,IWM&
     &,RES,ADDA,JAC,PJAC,SLVS)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: YH1(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: SAVR(*)
              REAL(KIND=8) :: ACOR(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL RES
              EXTERNAL ADDA
              EXTERNAL JAC
              EXTERNAL PJAC
              EXTERNAL SLVS
            END SUBROUTINE DSTODI
          END INTERFACE 
        END MODULE DSTODI__genmod
