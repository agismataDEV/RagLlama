        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DPCGS__genmod
          INTERFACE 
            SUBROUTINE DPCGS(NEQ,TN,Y,SAVF,R,WGHT,N,MAXL,DELTA,HL0,JPRE,&
     &MNEWT,F,PSOL,NPSL,X,P,W,Z,LPCG,WP,IWP,WK,IFLAG)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: TN
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: R(*)
              REAL(KIND=8) :: WGHT(*)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: MAXL
              REAL(KIND=8) :: DELTA
              REAL(KIND=8) :: HL0
              INTEGER(KIND=4) :: JPRE
              INTEGER(KIND=4) :: MNEWT
              EXTERNAL F
              EXTERNAL PSOL
              INTEGER(KIND=4) :: NPSL
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: P(*)
              REAL(KIND=8) :: W(*)
              REAL(KIND=8) :: Z(*)
              INTEGER(KIND=4) :: LPCG
              REAL(KIND=8) :: WP(*)
              INTEGER(KIND=4) :: IWP(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DPCGS
          END INTERFACE 
        END MODULE DPCGS__genmod
