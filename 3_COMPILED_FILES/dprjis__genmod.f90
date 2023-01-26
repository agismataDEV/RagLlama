        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPRJIS__genmod
          INTERFACE 
            SUBROUTINE DPRJIS(NEQ,Y,YH,NYH,EWT,RTEM,SAVR,S,WK,IWK,RES,  &
     &JAC,ADDA)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: RTEM(*)
              REAL(KIND=8) :: SAVR(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IWK(*)
              EXTERNAL RES
              EXTERNAL JAC
              EXTERNAL ADDA
            END SUBROUTINE DPRJIS
          END INTERFACE 
        END MODULE DPRJIS__genmod
