        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ODRV__genmod
          INTERFACE 
            SUBROUTINE ODRV(N,IA,JA,A,P,IP,NSP,ISP,PATH,FLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A(*)
              INTEGER(KIND=4) :: P(*)
              INTEGER(KIND=4) :: IP(*)
              INTEGER(KIND=4) :: NSP
              INTEGER(KIND=4) :: ISP(*)
              INTEGER(KIND=4) :: PATH
              INTEGER(KIND=4) :: FLAG
            END SUBROUTINE ODRV
          END INTERFACE 
        END MODULE ODRV__genmod
