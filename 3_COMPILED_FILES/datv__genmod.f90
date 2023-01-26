        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DATV__genmod
          INTERFACE 
            SUBROUTINE DATV(NEQ,Y,SAVF,V,WGHT,FTEM,F,PSOL,Z,VTEM,WP,IWP,&
     &HL0,JPRE,IER,NPSL)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: V(*)
              REAL(KIND=8) :: WGHT(*)
              REAL(KIND=8) :: FTEM(*)
              EXTERNAL F
              EXTERNAL PSOL
              REAL(KIND=8) :: Z(*)
              REAL(KIND=8) :: VTEM(*)
              REAL(KIND=8) :: WP(*)
              INTEGER(KIND=4) :: IWP(*)
              REAL(KIND=8) :: HL0
              INTEGER(KIND=4) :: JPRE
              INTEGER(KIND=4) :: IER
              INTEGER(KIND=4) :: NPSL
            END SUBROUTINE DATV
          END INTERFACE 
        END MODULE DATV__genmod
