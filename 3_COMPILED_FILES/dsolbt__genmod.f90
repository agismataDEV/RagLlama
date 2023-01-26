        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DSOLBT__genmod
          INTERFACE 
            SUBROUTINE DSOLBT(M,N,A,B,C,Y,IP)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: M
              REAL(KIND=8) :: A(M,M,N)
              REAL(KIND=8) :: B(M,M,N)
              REAL(KIND=8) :: C(M,M,N)
              REAL(KIND=8) :: Y(M,N)
              INTEGER(KIND=4) :: IP(M,N)
            END SUBROUTINE DSOLBT
          END INTERFACE 
        END MODULE DSOLBT__genmod
