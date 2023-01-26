        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIPREP__genmod
          INTERFACE 
            SUBROUTINE DIPREP(NEQ,Y,RWORK,IA,JA,IPFLAG,F,JAC)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: RWORK(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              INTEGER(KIND=4) :: IPFLAG
              EXTERNAL F
              EXTERNAL JAC
            END SUBROUTINE DIPREP
          END INTERFACE 
        END MODULE DIPREP__genmod
