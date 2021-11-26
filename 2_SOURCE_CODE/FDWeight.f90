MODULE FDWEIGHT

! =============================================================================
!
!   Programmer: Koos Huijssen
!
!   Language: Fortran 90
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The module FDWEIGHT contains the subroutine FDWEIGHTS that computes the
!   weights in a finite difference template, for arbitrary number of points
!   and on a not necessarily equidistant grid.
!
! *****************************************************************************
!
!   MODULES USED
!
!   none
!
! *****************************************************************************
!
!   GLOBAL DECLARATIONS
!
!   none
!
      IMPLICIT NONE
!
! *****************************************************************************
!
!   CONTAINED SUBROUTINES AND FUNCTIONS
!
!   FDWEIGHTS   sub    computes the weights of an FD template.
!
! =============================================================================
!
      CONTAINS

      SUBROUTINE FDWEIGHTS(XI,X,N,M,C)

! =============================================================================
!
!   Programmer: Koos Huijssen
!
!   Language: Fortran 77
!
!   Version Date    Comment
!   ------- -----   -------
!   1.0     090505  Original code (KH)
!
!   Source: B. Fornberg, 'A practical guide to pseudospectral methods'.
!           New York: Cambridge University Press, 1996
!
! *****************************************************************************
!
!   DESCRIPTION
!
!   The subroutine FDWEIGHTS computes the weights in a finite difference 
!   template, for arbitrary differentiation order, arbitrary number of points
!   and on a not necessarily equidistant grid.
!
!   Method: see B. Fornberg, 'A practical guide to pseudospectral methods', 
!   Chapter 3 and Appendix C   
!
! *****************************************************************************
!
!   INPUT/OUTPUT PARAMETERS
!
!   XI  i   dp    point at which the FD approximations are to be accurate
!   X   i   dp    X-coordinates for the grid points, array dimensioned X(0:N)
!   N   i   int   maximum index in X -> total number of grid points is N+1
!   M   i   int   highest order of the derivative to be approximated
!   C   o   dp    weights, array dimensioned C(0:N,0:M).
!                 On return, the element C(J,K) contains the weight to be
!                 applied at point X(J) when the K-th derivative is 
!                 approximated by a stencil extending over X(0),X(1),...,X(N).
! 
      REAL*8  XI
      REAL*8  X
      INTEGER N
      INTEGER M
      REAL*8  C
      DIMENSION X(0:N),C(0:N,0:M)

! *****************************************************************************
! 
!   LOCAL PARAMETERS      
!
!   C1   dp   temporary variable in dp
!   C2   dp   temporary variable in dp
!   C3   dp   temporary variable in dp
!   C4   dp   temporary variable in dp
!   C5   dp   temporary variable in dp
!   I    int  loop counter in int
!   J    int  loop counter in int
!   K    int  loop counter in int
!   MN   int  temporary variable in int
!   
      REAL*8 C1
      REAL*8 C2
      REAL*8 C3
      REAL*8 C4
      REAL*8 C5

      INTEGER I 
      INTEGER J
      INTEGER K
      INTEGER MN

! *****************************************************************************
!
!   I/O
!
!   none
!   
! *****************************************************************************
!
!   SUBROUTINES/FUNCTIONS CALLED
!
!   none
!   
! =============================================================================

      C1     = 1.0D0
      C4     = X(0)-XI
      DO 10 K=0,M
        DO 10 J=0,N
10        C(J,K)=0.0D0
      C(0,0) = 1.0D0
      DO 50 I=1,N
        MN = MIN(I,M)
        C2    = 1.0D0
        C5    = C4
        C4    = X(I)-XI
        DO 40 J=0,I-1
          C3 = X(I)-X(J)
          C2 = C2*C3
	      DO 20 K=MN,1,-1
20          C(I,K) = C1*(K*C(I-1,K-1)-C5*C(I-1,K))/C2
          C(I,0)     = -C1*C5*C(I-1,0)/C2
          DO 30 K = MN,1,-1
30	        C(J,K) = (C4*C(J,K)-K*C(J,K-1))/C3
40        C(J,0) = C4*C(J,0)/C3
50      C1=C2

      RETURN

      END SUBROUTINE

      END MODULE
