        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAINVGS__genmod
          INTERFACE 
            SUBROUTINE DAINVGS(NEQ,T,Y,WK,IWK,TEM,YDOT,IER,RES,ADDA)
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IWK(*)
              REAL(KIND=8) :: TEM(*)
              REAL(KIND=8) :: YDOT(*)
              INTEGER(KIND=4) :: IER
              EXTERNAL RES
              EXTERNAL ADDA
            END SUBROUTINE DAINVGS
          END INTERFACE 
        END MODULE DAINVGS__genmod
