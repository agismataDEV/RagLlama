        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAINVG__genmod
          INTERFACE 
            SUBROUTINE DAINVG(RES,ADDA,NEQ,T,Y,YDOT,MITER,ML,MU,PW,IPVT,&
     &IER)
              EXTERNAL RES
              EXTERNAL ADDA
              INTEGER(KIND=4) :: NEQ
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YDOT(*)
              INTEGER(KIND=4) :: MITER
              INTEGER(KIND=4) :: ML
              INTEGER(KIND=4) :: MU
              REAL(KIND=8) :: PW(*)
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: IER
            END SUBROUTINE DAINVG
          END INTERFACE 
        END MODULE DAINVG__genmod
