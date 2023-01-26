        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:41 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DDOT__genmod
          INTERFACE 
            FUNCTION DDOT(N,DX,INCX,DY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DY(*)
              INTEGER(KIND=4) :: INCY
              REAL(KIND=8) :: DDOT
            END FUNCTION DDOT
          END INTERFACE 
        END MODULE DDOT__genmod
