        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPREPI__genmod
          INTERFACE 
            SUBROUTINE DPREPI(NEQ,Y,S,YH,SAVR,EWT,RTEM,IA,JA,IC,JC,WK,  &
     &IWK,IPPER,RES,JAC,ADDA)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: YH(*)
              REAL(KIND=8) :: SAVR(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: RTEM(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              INTEGER(KIND=4) :: IC(*)
              INTEGER(KIND=4) :: JC(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IWK(*)
              INTEGER(KIND=4) :: IPPER
              EXTERNAL RES
              EXTERNAL JAC
              EXTERNAL ADDA
            END SUBROUTINE DPREPI
          END INTERFACE 
        END MODULE DPREPI__genmod
