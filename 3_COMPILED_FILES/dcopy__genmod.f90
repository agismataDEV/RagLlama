        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:41 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DCOPY__genmod
          INTERFACE 
            SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: DX(*)
              INTEGER(KIND=4) :: INCX
              REAL(KIND=8) :: DY(*)
              INTEGER(KIND=4) :: INCY
            END SUBROUTINE DCOPY
          END INTERFACE 
        END MODULE DCOPY__genmod
