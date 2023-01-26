        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DUSOL__genmod
          INTERFACE 
            SUBROUTINE DUSOL(NEQ,TN,Y,SAVF,B,WGHT,N,DELTA,HL0,MNEWT,PSOL&
     &,NPSL,X,WP,IWP,WK,IFLAG)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: TN
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: B(*)
              REAL(KIND=8) :: WGHT(*)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DELTA
              REAL(KIND=8) :: HL0
              INTEGER(KIND=4) :: MNEWT
              EXTERNAL PSOL
              INTEGER(KIND=4) :: NPSL
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: WP(*)
              INTEGER(KIND=4) :: IWP(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DUSOL
          END INTERFACE 
        END MODULE DUSOL__genmod
