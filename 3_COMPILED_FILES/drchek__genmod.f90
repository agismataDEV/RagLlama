        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DRCHEK__genmod
          INTERFACE 
            SUBROUTINE DRCHEK(JOB,G,NEQ,Y,YH,NYH,G0,G1,GX,JROOT,IRT)
              INTEGER(KIND=4) :: NYH
              INTEGER(KIND=4) :: JOB
              EXTERNAL G
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: G0(*)
              REAL(KIND=8) :: G1(*)
              REAL(KIND=8) :: GX(*)
              INTEGER(KIND=4) :: JROOT(*)
              INTEGER(KIND=4) :: IRT
            END SUBROUTINE DRCHEK
          END INTERFACE 
        END MODULE DRCHEK__genmod
