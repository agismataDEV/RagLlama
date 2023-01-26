        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DAIGBT__genmod
          INTERFACE 
            SUBROUTINE DAIGBT(RES,ADDA,NEQ,T,Y,YDOT,MB,NB,PW,IPVT,IER)
              EXTERNAL RES
              EXTERNAL ADDA
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: T
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YDOT(*)
              INTEGER(KIND=4) :: MB
              INTEGER(KIND=4) :: NB
              REAL(KIND=8) :: PW(*)
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: IER
            END SUBROUTINE DAIGBT
          END INTERFACE 
        END MODULE DAIGBT__genmod
