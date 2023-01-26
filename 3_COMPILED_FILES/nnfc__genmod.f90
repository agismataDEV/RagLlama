        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE NNFC__genmod
          INTERFACE 
            SUBROUTINE NNFC(N,R,C,IC,IA,JA,A,Z,B,LMAX,IL,JL,IJL,L,D,UMAX&
     &,IU,JU,IJU,U,ROW,TMP,IRL,JRL,FLAG)
              INTEGER(KIND=4) :: N
              INTEGER(KIND=4) :: R(*)
              INTEGER(KIND=4) :: C(*)
              INTEGER(KIND=4) :: IC(*)
              INTEGER(KIND=4) :: IA(*)
              INTEGER(KIND=4) :: JA(*)
              REAL(KIND=8) :: A(*)
              REAL(KIND=8) :: Z(*)
              REAL(KIND=8) :: B(*)
              INTEGER(KIND=4) :: LMAX
              INTEGER(KIND=4) :: IL(*)
              INTEGER(KIND=4) :: JL(*)
              INTEGER(KIND=4) :: IJL(*)
              REAL(KIND=8) :: L(*)
              REAL(KIND=8) :: D(*)
              INTEGER(KIND=4) :: UMAX
              INTEGER(KIND=4) :: IU(*)
              INTEGER(KIND=4) :: JU(*)
              INTEGER(KIND=4) :: IJU(*)
              REAL(KIND=8) :: U(*)
              REAL(KIND=8) :: ROW(*)
              REAL(KIND=8) :: TMP(*)
              INTEGER(KIND=4) :: IRL(*)
              INTEGER(KIND=4) :: JRL(*)
              INTEGER(KIND=4) :: FLAG
            END SUBROUTINE NNFC
          END INTERFACE 
        END MODULE NNFC__genmod
