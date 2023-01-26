        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DIPREPI__genmod
          INTERFACE 
            SUBROUTINE DIPREPI(NEQ,Y,S,RWORK,IA,JA,IC,JC,IPFLAG,RES,JAC,&
     &ADDA)
              INTEGER(KIND=4) :: NEQ(*)
              REAL(KIND=8) :: Y(*)
              REAL(KIND=8) :: S(*)
              REAL(KIND=8) :: RWORK(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              INTEGER(KIND=4) :: IC(*)
              INTEGER(KIND=4) :: JC(*)
              INTEGER(KIND=4) :: IPFLAG
              EXTERNAL RES
              EXTERNAL JAC
              EXTERNAL ADDA
            END SUBROUTINE DIPREPI
          END INTERFACE 
        END MODULE DIPREPI__genmod
