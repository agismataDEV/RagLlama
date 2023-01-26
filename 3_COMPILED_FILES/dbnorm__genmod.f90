        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DBNORM__genmod
          INTERFACE 
            FUNCTION DBNORM(N,A,NRA,ML,MU,W)
              INTEGER(KIND=4) :: NRA
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(NRA,N)
              INTEGER(KIND=4) :: ML
              INTEGER(KIND=4) :: MU
              REAL(KIND=8) :: W(N)
              REAL(KIND=8) :: DBNORM
            END FUNCTION DBNORM
          END INTERFACE 
        END MODULE DBNORM__genmod
