MODULE Types
	
    USE MPI_F08
    USE OMP_LIB
    
    IMPLICIT NONE
    SAVE

!DEC$ NOFREEFORM
!IBM* SOURCEFORM (FIXED)
		  !#ifdef COMPILE_LOCAL
	    	  !include "mpi.f"             ! Uncomment this , if USE MPI Not working
		  !#else
		  !include "mpif_ibm.h"
		  !#endif
!IBM* SOURCEFORM (FREE(F90))
!DEC$ FREEFORM

                !********************************************************************************
                !
                !       different number representations, derived from NRTYPES.F90
                !
                !       i   = integer
                !       sp  = single precision
                !       dp  = double precision
                !       dpc = complex
                !       lgt = logical type
                !
                !       The parameters with an suffix S, contain the memory-size of that kind
                !
                !********************************************************************************

				integer, parameter :: i8b   			= selected_int_kind(10)
                integer, parameter :: i4b               = selected_int_kind(9)
                integer, parameter :: i2b               = selected_int_kind(4)
                integer, parameter :: i1b               = selected_int_kind(2)
                integer(i8b), parameter :: i8bS         = sizeof(int(1,i8b))
                integer(i8b), parameter :: i4bS         = sizeof(int(1,i4b))
                integer(i8b), parameter :: i2bS         = sizeof(int(1,i2b))
                integer(i8b), parameter :: i1bS         = sizeof(int(1,i1b))

                integer, parameter :: sp                = kind(1.0)
                integer, parameter :: dp                = selected_real_kind(15, 307)
                integer, parameter :: qp                = selected_real_kind(32, 4931)
                integer(i8b), parameter :: spS          = sizeof(real(1, sp));
                integer(i8b), parameter :: dpS          = sizeof(real(1, dp));
                integer(i8b), parameter :: qpS          = sizeof(real(1, qp));

                integer, parameter :: spc               = kind((1.0,1.0))
                integer, parameter :: dpc               = kind((1.0d0,1.0d0))
                integer, parameter :: qpc               = kind((1.0q0,1.0q0))
                integer(i8b), parameter :: spcS         = sizeof(cmplx(1,1, spc));
                integer(i8b), parameter :: dpcS         = sizeof(cmplx(1,1, dpc));

                integer, parameter :: lgt               = kind(.true.)
                integer(i8b), parameter :: lgtS         = sizeof(.true.)
END MODULE Types
