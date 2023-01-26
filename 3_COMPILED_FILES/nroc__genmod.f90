        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NROC__genmod
          INTERFACE 
            SUBROUTINE NROC(N,IC,IA,JA,A,JAR,AR,P,FLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IC(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A(*)
              INTEGER(KIND=4) :: JAR(*)
              REAL(KIND=8) :: AR(*)
              INTEGER(KIND=4) :: P(*)
              INTEGER(KIND=4) :: FLAG
            END SUBROUTINE NROC
          END INTERFACE 
        END MODULE NROC__genmod
