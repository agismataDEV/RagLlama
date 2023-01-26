        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPRJS__genmod
          INTERFACE 
            SUBROUTINE DPRJS(NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WK,IWK,F,JAC)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: EWT(*)
              REAL(KIND=8) :: FTEM(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IWK(*)
              EXTERNAL F
              EXTERNAL JAC
            END SUBROUTINE DPRJS
          END INTERFACE 
        END MODULE DPRJS__genmod
