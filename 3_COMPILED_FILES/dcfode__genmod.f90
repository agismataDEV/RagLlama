        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DCFODE__genmod
          INTERFACE 
            SUBROUTINE DCFODE(METH,ELCO,TESCO)
              INTEGER(KIND=4) :: METH
              REAL(KIND=8) :: ELCO(13,12)
              REAL(KIND=8) :: TESCO(3,12)
            END SUBROUTINE DCFODE
          END INTERFACE 
        END MODULE DCFODE__genmod
