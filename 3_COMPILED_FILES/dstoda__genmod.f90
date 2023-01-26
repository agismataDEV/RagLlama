        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSTODA__genmod
          INTERFACE 
            SUBROUTINE DSTODA(NEQ,Y,YH,NYH,YH1,EWT,SAVF,ACOR,WM,IWM,F,  &
     &JAC,PJAC,SLVS)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: YH1(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: ACOR(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL F
              EXTERNAL JAC
              EXTERNAL PJAC
              EXTERNAL SLVS
            END SUBROUTINE DSTODA
          END INTERFACE 
        END MODULE DSTODA__genmod
