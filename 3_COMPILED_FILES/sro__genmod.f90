        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SRO__genmod
          INTERFACE 
            SUBROUTINE SRO(N,IP,IA,JA,A,Q,R,DFLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IP(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A(*)
              INTEGER(KIND=4) :: Q(*)
              INTEGER(KIND=4) :: R(*)
              LOGICAL(KIND=4) :: DFLAG
            END SUBROUTINE SRO
          END INTERFACE 
        END MODULE SRO__genmod
