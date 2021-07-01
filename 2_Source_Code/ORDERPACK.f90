MODULE ORDERPACK
USE Types

contains

Subroutine D_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! __________________________________________________________
      Real (kind=dp), Dimension (:), Intent (InOut) :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real (kind=dp), Dimension (Size(XDONT)) :: XINDT
      integer(i8b), Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      integer(i8b) :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call D_MRGRNK (XINDT, JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine D_ctrper

Subroutine R_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real, Dimension (Size(XDONT)) :: XINDT
      integer(i8b), Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      integer(i8b) :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call R_MRGRNK (XINDT, JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine R_ctrper
Subroutine I_ctrper (XDONT, PCLS)
!   Permute array XVALT randomly, but leaving elements close
!   to their initial locations (nearbyness is controled by PCLS).
! _________________________________________________________________
!   The routine takes the 1...size(XVALT) index array as real
!   values, takes a combination of these values and of random
!   values as a perturbation of the index array, and sorts the
!   initial set according to the ranks of these perturbated indices.
!   The relative proportion of initial order and random order
!   is 1-PCLS / PCLS, thus when PCLS = 0, there is no change in
!   the order whereas the new order is fully random when PCLS = 1.
!   Michel Olagnon - May 2000.
! _________________________________________________________________
! __________________________________________________________
      integer(i8b), Dimension (:), Intent (InOut)  :: XDONT
      Real, Intent (In) :: PCLS
! __________________________________________________________
!
      Real, Dimension (Size(XDONT)) :: XINDT
      integer(i8b), Dimension (Size(XDONT)) :: JWRKT
      Real :: PWRK
      integer(i8b) :: I
!
      Call Random_number (XINDT(:))
      PWRK = Min (Max (0.0, PCLS), 1.0)
      XINDT = Real(Size(XDONT)) * XINDT
      XINDT = PWRK*XINDT + (1.0-PWRK)*(/ (Real(I), I=1,size(XDONT)) /)
      Call I_MRGRNK (INT(XINDT,8), JWRKT)
      XDONT = XDONT (JWRKT)
!
End Subroutine I_ctrper
!================================================== UNISTA ==============================================

Subroutine D_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      Real (kind=dp), Dimension (:), Intent (InOut) :: XDONT
      integer(i8b), Intent (Out) :: NUNI
! __________________________________________________________
!
      integer(i8b), Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      integer(i8b) :: ICRS
! __________________________________________________________
      Call D_UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine D_unista

Subroutine R_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (InOut) :: XDONT
      integer(i8b), Intent (Out) :: NUNI
! __________________________________________________________
!
      integer(i8b), Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      integer(i8b) :: ICRS
! __________________________________________________________
      Call R_UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine R_unista
Subroutine I_unista (XDONT, NUNI)
!   UNISTA = (Stable unique) Removes duplicates from an array,
!            leaving unique entries in the order of their first
!            appearance in the initial set.
!  Michel Olagnon - Feb. 2000
! __________________________________________________________
! __________________________________________________________
      integer(i8b), Dimension (:), Intent (InOut)  :: XDONT
      integer(i8b), Intent (Out) :: NUNI
! __________________________________________________________
!
      integer(i8b), Dimension (Size(XDONT)) :: IWRKT
      Logical, Dimension (Size(XDONT)) :: IFMPTYT
      integer(i8b) :: ICRS
! __________________________________________________________
      Call I_UNIINV (XDONT, IWRKT)
      IFMPTYT = .True.
      NUNI = 0
      Do ICRS = 1, Size(XDONT)
         If (IFMPTYT(IWRKT(ICRS))) Then
            IFMPTYT(IWRKT(ICRS)) = .False.
            NUNI = NUNI + 1
            XDONT (NUNI) = XDONT (ICRS)
         End If
      End Do
      Return
!
End Subroutine I_unista

!================================================== MRGNK ==============================================


Subroutine D_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=dp), Dimension (:), Intent (In) :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real (kind=dp) :: XVALA, XVALB
!
      integer(i8b), Dimension (SIZE(IRNGT)) :: JWRKT
      integer(i8b) :: LMTNA, LMTNC, IRNG1, IRNG2
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine D_mrgrnk

Subroutine R_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      Real :: XVALA, XVALB
!
      integer(i8b), Dimension (SIZE(IRNGT)) :: JWRKT
      integer(i8b) :: LMTNA, LMTNC, IRNG1, IRNG2
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine R_mrgrnk
Subroutine I_mrgrnk (XDONT, IRNGT)
! __________________________________________________________
!   MRGRNK = Merge-sort ranking of an array
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      integer(i8b), Dimension (:), Intent (In)  :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IRNGT
! __________________________________________________________
      integer(i8b) :: XVALA, XVALB
!
      integer(i8b), Dimension (SIZE(IRNGT)) :: JWRKT
      integer(i8b) :: LMTNA, LMTNC, IRNG1, IRNG2
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IRNGT))
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IRNGT (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) <= XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo(NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 2) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!   Shortcut for the case when the max of A is smaller
!   than the min of B. This line may be activated when the
!   initial set is already close to sorted.
!
!          IF (XDONT(IRNGT(JINDA)) <= XDONT(IRNGT(IINDB))) CYCLE
!
!  One steps in the C subset, that we build in the final rank array
!
!  Make a copy of the rank array for the merge iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
!
            XVALA = XDONT (JWRKT(IINDA))
            XVALB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XVALA > XVALB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XVALB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XVALA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
      Return
!
End Subroutine I_mrgrnk
! =============================================================== UNIINV =====================================================

Subroutine D_uniinv (XDONT, IGOEST)
! __________________________________________________________
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      Real (kind=dp), Dimension (:), Intent (In) :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IGOEST
! __________________________________________________________
      Real (kind=dp) :: XTST, XDONA, XDONB
!
! __________________________________________________________
      integer(i8b), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
      integer(i8b) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IGOEST (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) < XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XDONA = XDONT (JWRKT(IINDA))
            XDONB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XDONA > XDONB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XDONB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XDONA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = D_NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
      Else
        XTST = D_NEARLESS (XDONT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XDONT(IRNG) > XTST) Then
            XTST = XDONT (IRNG)
            NUNI = NUNI + 1
         End If
         IGOEST (IRNG) = NUNI
!
      End Do
!
      Return
!
End Subroutine D_uniinv

Subroutine R_uniinv (XDONT, IGOEST)
! __________________________________________________________
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! _________________________________________________________
      Real, Dimension (:), Intent (In) :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IGOEST
! __________________________________________________________
      Real    :: XTST, XDONA, XDONB
!
! __________________________________________________________
      integer(i8b), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
      integer(i8b) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IGOEST (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) < XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XDONA = XDONT (JWRKT(IINDA))
            XDONB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XDONA > XDONB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XDONB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XDONA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = R_NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
      Else
        XTST = R_NEARLESS (XDONT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XDONT(IRNG) > XTST) Then
            XTST = XDONT (IRNG)
            NUNI = NUNI + 1
         End If
         IGOEST (IRNG) = NUNI
!
      End Do
!
      Return
!
End Subroutine R_uniinv
Subroutine I_uniinv (XDONT, IGOEST)
! __________________________________________________________
!   UNIINV = Merge-sort inverse ranking of an array, with removal of
!   duplicate entries.
!   The routine is similar to pure merge-sort ranking, but on
!   the last pass, it sets indices in IGOEST to the rank
!   of the value in the ordered set with duplicates removed.
!   For performance reasons, the first 2 passes are taken
!   out of the standard loop, and use dedicated coding.
! __________________________________________________________
! __________________________________________________________
      integer(i8b), Dimension (:), Intent (In)  :: XDONT
      integer(i8b), Dimension (:), Intent (Out) :: IGOEST
! __________________________________________________________
      integer(i8b) :: XTST, XDONA, XDONB
!
! __________________________________________________________
      integer(i8b), Dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
      integer(i8b) :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
      integer(i8b) :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB
!
      NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
!
      Select Case (NVAL)
      Case (:0)
         Return
      Case (1)
         IGOEST (1) = 1
         Return
      Case Default
         Continue
      End Select
!
!  Fill-in the index array, creating ordered couples
!
      Do IIND = 2, NVAL, 2
         If (XDONT(IIND-1) < XDONT(IIND)) Then
            IRNGT (IIND-1) = IIND - 1
            IRNGT (IIND) = IIND
         Else
            IRNGT (IIND-1) = IIND
            IRNGT (IIND) = IIND - 1
         End If
      End Do
      If (Modulo (NVAL, 2) /= 0) Then
         IRNGT (NVAL) = NVAL
      End If
!
!  We will now have ordered subsets A - B - A - B - ...
!  and merge A and B couples into     C   -   C   - ...
!
      LMTNA = 2
      LMTNC = 4
!
!  First iteration. The length of the ordered subsets goes from 2 to 4
!
      Do
         If (NVAL <= 4) Exit
!
!   Loop on merges of A and B into C
!
         Do IWRKD = 0, NVAL - 1, 4
            If ((IWRKD+4) > NVAL) Then
               If ((IWRKD+2) >= NVAL) Exit
!
!   1 2 3
!
               If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
!
!   1 3 2
!
               If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
                  IRNG2 = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNG2
!
!   3 1 2
!
               Else
                  IRNG1 = IRNGT (IWRKD+1)
                  IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
                  IRNGT (IWRKD+2) = IRNG1
               End If
               Exit
            End If
!
!   1 2 3 4
!
            If (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
!
!   1 3 x x
!
            If (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) Then
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
               If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   1 3 2 4
                  IRNGT (IWRKD+3) = IRNG2
               Else
!   1 3 4 2
                  IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+4) = IRNG2
               End If
!
!   3 x x x
!
            Else
               IRNG1 = IRNGT (IWRKD+1)
               IRNG2 = IRNGT (IWRKD+2)
               IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
               If (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) Then
                  IRNGT (IWRKD+2) = IRNG1
                  If (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) Then
!   3 1 2 4
                     IRNGT (IWRKD+3) = IRNG2
                  Else
!   3 1 4 2
                     IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
                     IRNGT (IWRKD+4) = IRNG2
                  End If
               Else
!   3 4 1 2
                  IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
                  IRNGT (IWRKD+3) = IRNG1
                  IRNGT (IWRKD+4) = IRNG2
               End If
            End If
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 4
         Exit
      End Do
!
!  Iteration loop. Each time, the length of the ordered subsets
!  is doubled.
!
      Do
         If (2*LMTNA >= NVAL) Exit
         IWRKF = 0
         LMTNC = 2 * LMTNC
!
!   Loop on merges of A and B into C
!
         Do
            IWRK = IWRKF
            IWRKD = IWRKF + 1
            JINDA = IWRKF + LMTNA
            IWRKF = IWRKF + LMTNC
            If (IWRKF >= NVAL) Then
               If (JINDA >= NVAL) Exit
               IWRKF = NVAL
            End If
            IINDA = 1
            IINDB = JINDA + 1
!
!  One steps in the C subset, that we create in the final rank array
!
!  Make a copy of the rank array for the iteration
!
            JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
            XDONA = XDONT (JWRKT(IINDA))
            XDONB = XDONT (IRNGT(IINDB))
!
            Do
               IWRK = IWRK + 1
!
!  We still have unprocessed values in both A and B
!
               If (XDONA > XDONB) Then
                  IRNGT (IWRK) = IRNGT (IINDB)
                  IINDB = IINDB + 1
                  If (IINDB > IWRKF) Then
!  Only A still with unprocessed values
                     IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
                     Exit
                  End If
                  XDONB = XDONT (IRNGT(IINDB))
               Else
                  IRNGT (IWRK) = JWRKT (IINDA)
                  IINDA = IINDA + 1
                  If (IINDA > LMTNA) Exit! Only B still with unprocessed values
                  XDONA = XDONT (JWRKT(IINDA))
               End If
!
            End Do
         End Do
!
!  The Cs become As and Bs
!
         LMTNA = 2 * LMTNA
      End Do
!
!   Last merge of A and B into C, with removal of duplicates.
!
      IINDA = 1
      IINDB = LMTNA + 1
      NUNI = 0
!
!  One steps in the C subset, that we create in the final rank array
!
      JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
      If (IINDB <= NVAL) Then
        XTST = I_NEARLESS (Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
      Else
        XTST = I_NEARLESS (XDONT(JWRKT(1)))
      Endif
      Do IWRK = 1, NVAL
!
!  We still have unprocessed values in both A and B
!
         If (IINDA <= LMTNA) Then
            If (IINDB <= NVAL) Then
               If (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) Then
                  IRNG = IRNGT (IINDB)
                  IINDB = IINDB + 1
               Else
                  IRNG = JWRKT (IINDA)
                  IINDA = IINDA + 1
               End If
            Else
!
!  Only A still with unprocessed values
!
               IRNG = JWRKT (IINDA)
               IINDA = IINDA + 1
            End If
         Else
!
!  Only B still with unprocessed values
!
            IRNG = IRNGT (IWRK)
         End If
         If (XDONT(IRNG) > XTST) Then
            XTST = XDONT (IRNG)
            NUNI = NUNI + 1
         End If
         IGOEST (IRNG) = NUNI
!
      End Do
!
      Return
!
End Subroutine I_uniinv

Function D_nearless (XVAL) result (D_nl)
!  Nearest value less than given value
! __________________________________________________________
      Real (kind=dp), Intent (In) :: XVAL
      Real (kind=dp) :: D_nl
! __________________________________________________________
      D_nl = nearest (XVAL, -1.0_dp)
      return
!
End Function D_nearless
Function R_nearless (XVAL) result (R_nl)
!  Nearest value less than given value
! __________________________________________________________
      Real, Intent (In) :: XVAL
      Real :: R_nl
! __________________________________________________________
      R_nl = nearest (XVAL, -1.0)
      return
!
End Function R_nearless
Function I_nearless (XVAL) result (I_nl)
!  Nearest value less than given value
! __________________________________________________________
      integer(i8b), Intent (In) :: XVAL
      integer(i8b) :: I_nl
! __________________________________________________________
      I_nl = XVAL - 1
      return
!
End Function I_nearless


RECURSIVE SUBROUTINE QUICKSORTHOARE(SortMe, lo, hi)
    implicit none
    INTEGER(i8b), intent(inout)               :: SortMe(*)
    INTEGER(i8b), intent(in)                  :: lo
    INTEGER(i8b), intent(in)                  :: hi
    INTEGER(i8b)                              :: Pivot

    IF(lo < hi) THEN
        Pivot = HOARE_PARTITION(SortMe, lo, hi)
        call QUICKSORTHOARE(SortMe, lo, Pivot)
        call QUICKSORTHOARE(SortMe, Pivot+1, hi)
    END IF

END SUBROUTINE QUICKSORTHOARE

FUNCTION HOARE_PARTITION(SortMe, lo, hi)
    INTEGER(i8b), dimension(*), intent(InOut) :: SortMe
    INTEGER(i8b), intent(in)                  :: lo
    INTEGER(i8b), intent(in)                  :: hi
    INTEGER(i8b)                              :: HOARE_PARTITION

    INTEGER(i8b) :: intTemp
    INTEGER(i8b) :: Pivot
    INTEGER(i8b) :: i, j

    Pivot = SortMe(lo+((hi-lo)/2))
    i = lo
    j = hi

    DO WHILE(.TRUE.)
        DO WHILE (SortMe(i) < Pivot)
            i = i + 1
        END DO

        DO WHILE (Pivot < SortMe(j))
            j = j - 1
        END DO

        IF (i >= j) THEN
            HOARE_PARTITION = j
            EXIT
        END IF

        intTemp   = SortMe(i)
        SortMe(i) = SortMe(j)
        SortMe(j) = intTemp

    END DO

END FUNCTION HOARE_PARTITION

END MODULE ORDERPACK
