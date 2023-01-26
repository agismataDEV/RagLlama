        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSPIGMR__genmod
          INTERFACE 
            SUBROUTINE DSPIGMR(NEQ,TN,Y,SAVF,B,WGHT,N,MAXL,MAXLP1,KMP,  &
     &DELTA,HL0,JPRE,MNEWT,F,PSOL,NPSL,X,V,HES,Q,LGMR,WP,IWP,WK,DL,IFLAG&
     &)
              INTEGER(KIND=4) :: MAXLP1
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: TN
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: B(*)
              REAL(KIND=8) :: WGHT(*)
              INTEGER(KIND=4) :: MAXL
              INTEGER(KIND=4) :: KMP
              REAL(KIND=8) :: DELTA
              REAL(KIND=8) :: HL0
              INTEGER(KIND=4) :: JPRE
              INTEGER(KIND=4) :: MNEWT
              EXTERNAL F
              EXTERNAL PSOL
              INTEGER(KIND=4) :: NPSL
              REAL(KIND=8) :: X(*)
              REAL(KIND=8) :: V(N,*)
              REAL(KIND=8) :: HES(MAXLP1,*)
              REAL(KIND=8) :: Q(*)
              INTEGER(KIND=4) :: LGMR
              REAL(KIND=8) :: WP(*)
              INTEGER(KIND=4) :: IWP(*)
              REAL(KIND=8) :: WK(*)
              REAL(KIND=8) :: DL(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DSPIGMR
          END INTERFACE 
        END MODULE DSPIGMR__genmod
