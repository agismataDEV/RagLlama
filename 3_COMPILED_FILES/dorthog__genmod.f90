        !COMPILER-GENERATED INTERFACE MODULE: Mon May 23 17:39:40 2022
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DORTHOG__genmod
          INTERFACE 
            SUBROUTINE DORTHOG(VNEW,V,HES,N,LL,LDHES,KMP,SNORMW)
              INTEGER(KIND=4) :: LDHES
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: VNEW(*)
              REAL(KIND=8) :: V(N,*)
              REAL(KIND=8) :: HES(LDHES,*)
              INTEGER(KIND=4) :: LL
              INTEGER(KIND=4) :: KMP
              REAL(KIND=8) :: SNORMW
            END SUBROUTINE DORTHOG
          END INTERFACE 
        END MODULE DORTHOG__genmod
