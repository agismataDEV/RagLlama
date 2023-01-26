        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DINTDY__genmod
          INTERFACE 
            SUBROUTINE DINTDY(T,K,YH,NYH,DKY,IFLAG)
              INTEGER(KIND=4) :: NYH
              REAL(KIND=8) :: T
              INTEGER(KIND=4) :: K
              REAL(KIND=8) :: YH(NYH,*)
              REAL(KIND=8) :: DKY(*)
              INTEGER(KIND=4) :: IFLAG
            END SUBROUTINE DINTDY
          END INTERFACE 
        END MODULE DINTDY__genmod
