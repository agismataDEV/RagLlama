        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DROOTS__genmod
          INTERFACE 
            SUBROUTINE DROOTS(NG,HMIN,JFLAG,X0,X1,G0,G1,GX,X,JROOT)
              INTEGER(KIND=4) :: NG
              REAL(KIND=8) :: HMIN
              INTEGER(KIND=4) :: JFLAG
              REAL(KIND=8) :: X0
              REAL(KIND=8) :: X1
              REAL(KIND=8) :: G0(NG)
              REAL(KIND=8) :: G1(NG)
              REAL(KIND=8) :: GX(NG)
              REAL(KIND=8) :: X
              INTEGER(KIND=4) :: JROOT(NG)
            END SUBROUTINE DROOTS
          END INTERFACE 
        END MODULE DROOTS__genmod
