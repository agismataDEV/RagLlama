        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPREPJI__genmod
          INTERFACE 
            SUBROUTINE DPREPJI(NEQ,Y,YH,NYH,EWT,RTEM,SAVR,S,WM,IWM,RES, &
     &JAC,ADDA)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: RTEM(*)
              REAL(KIND=8) :: SAVR(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: WM(*)
              INTEGER(KIND=4) :: IWM(*)
              EXTERNAL RES
              EXTERNAL JAC
              EXTERNAL ADDA
            END SUBROUTINE DPREPJI
          END INTERFACE 
        END MODULE DPREPJI__genmod
