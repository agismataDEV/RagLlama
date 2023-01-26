        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:39 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CDRV__genmod
          INTERFACE 
            SUBROUTINE CDRV(N,R,C,IC,IA,JA,A,B,Z,NSP,ISP,RSP,ESP,PATH,  &
     &FLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: R(*)
              INTEGER(KIND=4) :: C(*)
              INTEGER(KIND=4) :: IC(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A(*)
              REAL(KIND=8) :: B(*)
              REAL(KIND=8) :: Z(*)
              INTEGER(KIND=4) :: NSP
              INTEGER(KIND=4) :: ISP(*)
              REAL(KIND=8) :: RSP(*)
              INTEGER(KIND=4) :: ESP
              INTEGER(KIND=4) :: PATH
              INTEGER(KIND=4) :: FLAG
            END SUBROUTINE CDRV
          END INTERFACE 
        END MODULE CDRV__genmod
