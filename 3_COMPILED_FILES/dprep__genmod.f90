        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPREP__genmod
          INTERFACE 
            SUBROUTINE DPREP(NEQ,Y,YH,SAVF,EWT,FTEM,IA,JA,WK,IWK,IPPER,F&
     &,JAC)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: FTEM(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IWK(*)
              INTEGER(KIND=4) :: IPPER
              EXTERNAL F
              EXTERNAL JAC
            END SUBROUTINE DPREP
          END INTERFACE 
        END MODULE DPREP__genmod
