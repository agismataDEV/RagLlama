        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DFNORM__genmod
          INTERFACE 
            FUNCTION DFNORM(N,A,W)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(N,N)
              REAL(KIND=8) :: W(N)
              REAL(KIND=8) :: DFNORM
            END FUNCTION DFNORM
          END INTERFACE 
        END MODULE DFNORM__genmod
