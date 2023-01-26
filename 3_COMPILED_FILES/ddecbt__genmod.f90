        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DDECBT__genmod
          INTERFACE 
            SUBROUTINE DDECBT(M,N,A,B,C,IP,IER)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,M,N)
              REAL(KIND=8) :: B(M,M,N)
              REAL(KIND=8) :: C(M,M,N)
              INTEGER(KIND=4) :: IP(M,N)
              INTEGER(KIND=4) :: IER
            END SUBROUTINE DDECBT
          END INTERFACE 
        END MODULE DDECBT__genmod
