        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSPIOM__genmod
          INTERFACE 
            SUBROUTINE DSPIOM(NEQ,TN,Y,SAVF,B,WGHT,N,MAXL,KMP,DELTA,HL0,&
     &JPRE,MNEWT,F,PSOL,NPSL,X,V,HES,IPVT,LIOM,WP,IWP,WK,IFLAG)
              INTEGER(KIND=4) :: MAXL
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: TN
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: SAVF(*)
              REAL(KIND=8) :: B(*)
              REAL(KIND=8) :: WGHT(*)
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
              REAL(KIND=8) :: HES(MAXL,MAXL)
              INTEGER(KIND=4) :: IPVT(*)
              INTEGER(KIND=4) :: LIOM
              REAL(KIND=8) :: WP(*)
              INTEGER(KIND=4) :: IWP(*)
              REAL(KIND=8) :: WK(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DSPIOM
          END INTERFACE 
        END MODULE DSPIOM__genmod
